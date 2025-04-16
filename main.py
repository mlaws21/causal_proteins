# controller -- just calls into all of my bs

# Todo this could be a chatbot lol
# Another idea could be LLM aided generation of mutation parsing/extracting code

# things we can do

# generate data
import pandas as pd
import numpy as np
import threading
import os
import pickle
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # 0 = all logs, 1 = filter INFO, 2 = filter INFO+WARNING, 3 = filter all but ERROR

from typing import List, Dict, Optional, Tuple
from data_generation.mutgen import generate_synth, split_and_unique
from alphafold_driver.driver import fold_all
from sequence_adjust.proteinbert import protein_bert_scores
from sequence_adjust.blosum import blosum_scores
from alphafold_driver.alt_align import align_all
from datetime import datetime
import argparse




log_lock = threading.Lock()

def normalize(nums):
    return (nums - np.mean(nums)) / np.std(nums)
# before this need to create this list
# note need to ensure ids are filename safe
# name = name.replace(":", "_").replace(">", "_").lower()
def generate_data(ref_sequence: str, mutated_sequences: str | Tuple[list[str], list[str], list[bool]], sequence_score_method, project_name: str, coeffs, protein_length=253, num_datapoints=1000):
    '''
    TODO

    Parameters:
    - ref_sequence (str): A reference sequence or a path to a file containing it.
    - mutated_sequences (str | Tuple[list[str], list[str], list[str]]): Either a csv file containing 
      mutated sequences, clinical significance, or a tuple of the list of mutation ids, mutated 
      sequences, and significance (disease=True) directly.
    - sequence_score_method (str): ["proteinbert", "blosum"] or default to none
    - project_name (str)
    - coeffs
    - protein_length (int)
    - num_datapoints (int, default: 1000): number of datapoints to generate

    Returns:
    TODO
    '''
    
    with open(f"{project_name}_generate.log", "w"):
        pass
    
    
    def log(tolog):
        with log_lock:
            with open(f"{project_name}_generate.log", "a") as log_file:
                print(tolog, file=log_file)

    log(datetime.now())
    
    if ref_sequence.endswith(".txt"):
        with open(ref_sequence, 'r') as file:
            ref_sequence = file.readline().strip()
        
    
    if isinstance(mutated_sequences, str):
        mutated_data = pd.read_csv(mutated_sequences)
        mutated_data.columns = ['ID', 'sequence', 'clinical_sig']

    else:
        ids, sequences, clinical_sigs = mutated_sequences
        mutated_data = pd.DataFrame({
            'ID': ids,
            'sequence': sequences,
            'clinical_sig': clinical_sigs
        })
      
    # a bit ugly but ok  
    true_vals = {'true', '1', 'yes'}
    mutated_data['clinical_sig'] = mutated_data['clinical_sig'].map(
        lambda x: str(x).strip().lower() in true_vals
    )

    # adding reference seq to dataset.
    mutated_data.loc[len(mutated_data)] = ["ref", ref_sequence, False]
    
    unique_benign, unique_disease = split_and_unique(mutated_data)

    
    log(f"{len(unique_benign)} Benign Sequences Found")
    log(f"{len(unique_disease)} Benign Sequences Found")
    
    patient_data = generate_synth(unique_benign, unique_disease, num_datapoints=num_datapoints)
    log(f"{num_datapoints} rows generated")
    
    patient_data.loc[patient_data["Sequence"] == ref_sequence, "ID"] = "ref"
    
    
    # remember to normalize
    scores = []
    if sequence_score_method.lower() == "proteinbert":
        proteinbert_pkl = f"pickles/{project_name}_proteinbert_scores.pkl"
        if os.path.exists(proteinbert_pkl):
            log(f"Reading ProteinBERT scores from file {proteinbert_pkl}")
            with open(proteinbert_pkl, "rb") as f:
                scores = pickle.load(f)
            log(f"{len(scores)} ProteinBERT scores read in")
        else:
            log("Computing ProteinBERT Scores")
            scores = protein_bert_scores(ref_sequence, patient_data["Sequence"], protein_length)
            with open(proteinbert_pkl, "wb") as f:
                pickle.dump(scores, f)
            log(f"{len(scores)} Scores Computed")
    elif sequence_score_method.lower() == "blosum":
        blosum_pkl = f"pickles/{project_name}_blosum_scores.pkl"
        if os.path.exists(blosum_pkl):
            log(f"Reading BLOSUM scores from file {blosum_pkl}")
            with open(blosum_pkl, "rb") as f:
                scores = pickle.load(f)
            log(f"{len(scores)} BLOSUM scores read in")
        else:
            log("Computing BLOSUM Scores")
            scores = blosum_scores(ref_sequence, patient_data["Sequence"], log)
            with open(blosum_pkl, "wb") as f:
                pickle.dump(scores, f)
            log(f"{len(scores)} Scores Computed")
    else:
        scores = [0] * num_datapoints
    
    np_scores = np.array(scores)
    normalized_scores = 0 if np.std(np_scores) == 0 else normalize(np_scores)
    patient_data["Sequence_Score"] = normalized_scores
    
    fold_all(patient_data["ID"], patient_data["Sequence"], "b", project_name, log, num_workers=7)
    
    align_pkl = f"pickles/{project_name}_alignments.pkl"
    if os.path.exists(align_pkl):
        log(f"Reading alignments from file {align_pkl}")
        with open(align_pkl, "rb") as f:
            alignments = pickle.load(f)
        log(f"{len(alignments)} alignments read in")
        
    else:
        log("Aligning all proteins")
        alignments = align_all("ref", patient_data["ID"], project_name, log)
        with open(align_pkl, "wb") as f:
            pickle.dump(alignments, f)
        log(f"Aligned {len(alignments)} proteins")
    
    np_alignments = np.array(alignments)
    
    #optional normalize align scores
    # np_alignments = normalize(np_alignments)
    
    patient_data["Align_Score"] = np_alignments
    
    log("Calculating Disease Column")
        
    logit = (
        coeffs["Intercept"]
        + coeffs["Age"] * patient_data["Age"]
        + coeffs["Race"] * patient_data["Race"]
        + coeffs["Lifestyle"] * patient_data["Lifestyle"]
        + coeffs["Sequence_Score"] * patient_data["Sequence_Score"]
        + coeffs["Align_Score"] * patient_data["Align_Score"]
        
    )
    
    # print(logit)
    prob_Y = 1 / (1 + np.exp(-logit))  # Sigmoid function
    patient_data["Disease"] = np.random.binomial(1, prob_Y)

    log(f"Mean of Probability of Disease: {np.mean(patient_data['Disease'])}")
    
    patient_data.to_csv(f'{project_name}_cov_data.csv', index=False)
    
    # Now we generate the dose response curve
    


def main():
    parser = argparse.ArgumentParser(description="Generate prion data with given coefficients.")
    parser.add_argument("--project", default="prion411", help="Name of the project")
    parser.add_argument("--ref", default="prion_ref.txt", help="Reference file")
    parser.add_argument("--input", default="prion_input.csv", help="Input file")
    parser.add_argument("--metric", default="blosum", help="Score type (e.g., blosum)")

    args = parser.parse_args()

    coeffs = {"Intercept": -3, "Age": 0.5, "Race": 0.1, "Lifestyle": 0.3, "Sequence_Score": 0.7, "Align_Score": 0.8}

    generate_data(args.ref, args.input, args.metric, args.project, coeffs)


if __name__ == "__main__":
    main()
