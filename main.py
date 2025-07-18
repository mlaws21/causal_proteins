NUM_ROWS = 1000
NUM_SUBSET = None
OUTCOME = "Grade"

# controller -- just calls into all of my bs

# Todo this could be a chatbot lol
# Another idea could be LLM aided generation of mutation parsing/extracting code

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # 0 = all logs, 1 = filter INFO, 2 = filter INFO+WARNING, 3 = filter all but ERROR

# generate data
import pandas as pd
import numpy as np
import threading


import pickle
import argparse
from datetime import datetime
from typing import List, Dict, Optional, Tuple

from helper.generate_data import generate_synth, split_and_unique
from helper.fold import fold_all
from helper.proteinbert import protein_bert_scores
from helper.blosum import blosum_scores
from helper.dose_response import generate_dr_curve
from alignment_functions.alignment_functions import align_all, tm_align_all, tm_align_wrapper, rmsd_align, findPandQ_refonly, findPandQ_averaged
from helper.intervention import process_mutations
from helper.calc_effect import process_files
from helper.results import ordering, generate_boxplot, generate_pr_curve, generate_summary

np.random.seed(42)
GPUS = 8
log_lock = threading.Lock()


def normalize(nums):
    return np.round((nums - np.mean(nums)) / np.std(nums), decimals=4)
# before this need to create this list
# note need to ensure ids are filename safe
# name = name.replace(":", "_").replace(">", "_").lower()
# FIXME change back to 1000
def generate_data(ref_sequence: str, mutated_sequences: str | Tuple[list[str], list[str], list[bool]], sequence_score_method, project_name: str, protein_name: str, coeffs, protein_length=253, num_datapoints=NUM_ROWS):
    '''
    TODO

    Parameters:
    - ref_sequence (str): A reference sequence or a path to a file containing it.
    - mutated_sequences (str | Tuple[list[str], list[str], list[str]]): Either a csv file containing 
      mutated sequences, clinical significance, or a tuple of the list of mutation ids, mutated 
      sequences, and significance (disease=True) directly.
    - sequence_score_method (str): ["proteinbert", "blosum"] or default to none
    - project_name (str)
    - protein_name (str)
    - coeffs
    - protein_length (int)
    - num_datapoints (int, default: 1000): number of datapoints to generate

    Returns:
    TODO
    '''
    
    with open(f"outputs/{project_name}/generate.log", "w"):
        pass
    
    
    def log(tolog):
        with log_lock:
            with open(f"outputs/{project_name}/generate.log", "a") as log_file:
                print(tolog, file=log_file)

    log(datetime.now())
    
    if ref_sequence.endswith(".txt"):
        with open(ref_sequence, 'r') as file:
            ref_sequence = file.readline().strip()
    
    data_filepath = f'outputs/{project_name}/data.csv'
    
    numerical_cols = list(coeffs.keys())
    numerical_cols.remove("Intercept")
    
    if os.path.exists(data_filepath):
        patient_data = pd.read_csv(data_filepath)
        log(f"Reading Data from {data_filepath}")
    else:
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
            proteinbert_pkl = f"outputs/{project_name}/pickles/proteinbert_scores.pkl"
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
            blosum_pkl = f"outputs/{project_name}/pickles/blosum_scores.pkl"
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
        
        fold_all(patient_data["ID"], patient_data["Sequence"], protein_name, log, num_workers=GPUS)
        
        align_pkl = f"outputs/{project_name}/pickles/alignments.pkl"
        if os.path.exists(align_pkl):
            log(f"Reading alignments from file {align_pkl}")
            with open(align_pkl, "rb") as f:
                alignments = pickle.load(f)
            log(f"{len(alignments)} alignments read in")
            
        else:
            log("Aligning all proteins")
            # alignments = align_all("ref", patient_data["ID"], protein_name, rmsd_align, findPandQ_refonly, confidence_weighted=False, reweighted=True, log_fn=log)
            alignments = tm_align_all("ref", patient_data["ID"], protein_name, log)
            with open(align_pkl, "wb") as f:
                pickle.dump(alignments, f)
            log(f"Aligned {len(alignments)} proteins")
        
        np_alignments = np.array(alignments)
        
        #optional normalize align scores
        # np_alignments = normalize(np_alignments)
        
        patient_data["Align_Score"] = np_alignments
        
        log(f"Calculating {OUTCOME} Column")
        
        logit = coeffs["Intercept"] + sum([coeffs[x] * patient_data[x] for x in numerical_cols])
        
        # print(logit)
        prob_Y = 1 / (1 + np.exp(-logit))  # Sigmoid function
        patient_data[OUTCOME] = np.random.binomial(1, prob_Y)

        log(f"Mean of Probability of OUTCOME: {np.mean(patient_data[OUTCOME])}")
        
        patient_data.to_csv(f'outputs/{project_name}/data.csv', index=False)
    
    # Now we generate the dose response curve
    
    # Take all columns from patient_data and drop any non-numerical ones
    numerical_cols.append(OUTCOME)
    numerical_data = patient_data[numerical_cols]
    
    
    treatment = "Align_Score"
    # outcome = "Disease"
    generate_dr_curve(numerical_data, project_name, coeffs, log, treatment, OUTCOME)
    
    log("Data Generation Complete")
    
def use_existing_data(ref_sequence, seqs_and_cov_data, sequence_score_method, project_name: str, protein_name: str, protein_length=414):

    '''
    TODO

    Parameters:
    - ref_sequence (str): A reference sequence or a path to a file containing it.
    - mutated_sequences (str | Tuple[list[str], list[str], list[str]]): Either a csv file containing 
      mutated sequences, clinical significance, or a tuple of the list of mutation ids, mutated 
      sequences, and significance (disease=True) directly.
    - sequence_score_method (str): ["proteinbert", "blosum"] or default to none
    - project_name (str)
    - protein_name (str)
    - coeffs
    - protein_length (int)
    - num_datapoints (int, default: 1000): number of datapoints to generate

    Returns:
    TODO
    '''
    
    with open(f"outputs/{project_name}/generate.log", "w"):
        pass
    
    
    def log(tolog):
        with log_lock:
            with open(f"outputs/{project_name}/generate.log", "a") as log_file:
                print(tolog, file=log_file)

    log(datetime.now())
    
    if ref_sequence.endswith(".txt"):
        with open(ref_sequence, 'r') as file:
            ref_sequence = file.readline().strip()
    
    data_filepath = f'outputs/{project_name}/data.csv'
    
    if os.path.exists(data_filepath):
        patient_data = pd.read_csv(data_filepath)
        log(f"Reading Data from {data_filepath}")
    else:
        
        patient_data = pd.read_csv(seqs_and_cov_data)
        
        if "ID" in patient_data.columns:
            patient_data["ID"] = patient_data["ID"].astype(str).str.lower()
        else:
            raise ValueError("ID column not found in data file")
        
        num_datapoints = len(patient_data)
        
        
        # remember to normalize
        scores = []
        if sequence_score_method.lower() == "proteinbert":
            proteinbert_pkl = f"outputs/{project_name}/pickles/proteinbert_scores.pkl"
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
            blosum_pkl = f"outputs/{project_name}/pickles/blosum_scores.pkl"
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
        
        fold_all(patient_data["ID"], patient_data["Sequence"], protein_name, log, num_workers=GPUS)
        
        
        align_pkl = f"outputs/{project_name}/pickles/alignments.pkl"
        if os.path.exists(align_pkl):
            log(f"Reading alignments from file {align_pkl}")
            with open(align_pkl, "rb") as f:
                alignments = pickle.load(f)
            log(f"{len(alignments)} alignments read in")
            
        else:
            log("Aligning all proteins")
            # alignments = align_all("ref", patient_data["ID"], protein_name, rmsd_align, findPandQ_refonly, confidence_weighted=False, reweighted=True, log_fn=log)
            alignments = tm_align_all("ref", patient_data["ID"], protein_name, log)
            with open(align_pkl, "wb") as f:
                pickle.dump(alignments, f)
            log(f"Aligned {len(alignments)} proteins")
        
        np_alignments = np.array(alignments)
        
        #optional normalize align scores
        # np_alignments = normalize(np_alignments)
        
        patient_data["Align_Score"] = np_alignments


        log(f"Mean of Probability of OUTCOME: {np.mean(patient_data[OUTCOME])}")
        
    
    
    # Now we generate the dose response curve
    
    numerical_data = patient_data.select_dtypes(include=[np.number])
    log(f"Numerical Columns: {numerical_data.columns}")
    # outcome = "Disease"
    # numerical_data = numerical_data.sample(frac=1, random_state=42).reset_index(drop=True)
    # numerical_data = numerical_data.head(300)
    if "Ground" not in patient_data.columns:
        patient_data["Ground"] = 0
    patient_data.to_csv(f'outputs/{project_name}/data.csv', index=False)
    
    
    treatment = "Align_Score"
    
    generate_dr_curve(numerical_data, project_name, None, log, treatment, OUTCOME)
    
    log("Data Generation Complete")

# FIXME CHANGE BACK TO 100
def analyze_data(data_csv, spline_pkl, ref_sequence, project_name, protein_name, chosen_mutations=None, subset=NUM_SUBSET):
    
    with open(f"outputs/{project_name}/analysis.log", "w"):
        pass
    
    
    def log(tolog):
        with log_lock:
            with open(f"outputs/{project_name}/analysis.log", "a") as log_file:
                print(tolog, file=log_file)

    log(datetime.now())
    

    
    if ref_sequence.endswith(".txt"):
        with open(ref_sequence, 'r') as file:
            ref_sequence = file.readline().strip()
    # pandq finder, confidence weighted, reweight
    # FIXME
    alignment_args = ()
    process_mutations(data_csv, spline_pkl, ref_sequence, subset, project_name, protein_name, tm_align_wrapper, alignment_args, log, chosen_mutations, num_workers=GPUS)
    log("Analysis Complete")
    
def calculate_effect(project_name):
    
    with open(f"outputs/{project_name}/effect.log", "w"):
        pass
    
    
    def log(tolog):
        with log_lock:
            with open(f"outputs/{project_name}/effect.log", "a") as log_file:
                print(tolog, file=log_file)
    
    log(datetime.now())
    treatment = "Align_Score"
    log("Effect Calculation Complete")
    
    data = process_files(project_name, treatment, log)
    
    ordering(project_name, log)
    log("Generating Boxplot")
    generate_boxplot(project_name)
    log("Generating Precision-Recall Curve")
    generate_pr_curve(project_name)
    log("Generating Summary")
    generate_summary(project_name, threshold=0.0)
    log("Result Generation Complete")
    
    
def main():
    # NOTE can add defaults
    parser = argparse.ArgumentParser(description="Generate prion data with given coefficients.")
    parser.add_argument("--project", help="Name of the project")
    parser.add_argument("--protein", help="Name of Protein")
    parser.add_argument("--ref", help="Reference file")
    parser.add_argument("--input", help="Input file")
    parser.add_argument("--metric", default="none", help="Score type (e.g., blosum)")
    parser.add_argument('--no-interactive', action='store_false', dest='interactive', help='Disable UI')


    args = parser.parse_args()
    
    coeffs = {"Intercept": -1.75, "Age": 0.5, "Race": 0.1, "Lifestyle": 0.3, "Sequence_Score": 0.7, "Align_Score": 0.8}
    
    response = ""
    
    # initialize all the necesary directories
    if not os.path.exists(f"outputs/{args.project}"):
        os.makedirs(f"outputs/{args.project}")
        
    if not os.path.exists(f"outputs/{args.project}/pickles"):
        os.makedirs(f"outputs/{args.project}/pickles")
        
    if not os.path.exists(f"outputs/{args.project}/results"):
        os.makedirs(f"outputs/{args.project}/results")
        
    if not os.path.exists(f"outputs/{args.project}/intervention_data"):
        os.makedirs(f"outputs/{args.project}/intervention_data")
        
    if not os.path.exists(f"slurm/"):
        os.makedirs(f"slurm/")
    
    if args.interactive:
        while True:
            print(f"Welcome to Project [{args.project}]! Would you like to:")
            print("- Generate Data [1a]")
            print("- Use Existing Data [1b]")
            print("- Analyze Data [2]")
            print("- Calculate Effect [3]")
            print("- Help [4]")
            print("- quit [5]")
            response = input("> ")
            
            if response == '1a':
                print(f"Generating Data Effects and Logging Output to [outputs/{args.project}/generate.log]")
                generate_data(args.ref, args.input, args.metric, args.project, args.protein, coeffs)
                
            elif response == '1b':
                print(f"Using Existing Data and Logging Output to [outputs/{args.project}/generate.log]")
                use_existing_data(args.ref, args.input, args.metric, args.project, args.protein)
                    
            elif response == '2':
                chosen_mutations = input("Enter chosen mutations (comma-separated), or press Enter to use all: ").strip()
                if not chosen_mutations:
                    chosen_mutations = None
                else:
                    chosen_mutations = [m.lower().strip() for m in chosen_mutations.split(",")]
                print(f"Analyzing Data and Logging Output to [outputs/{args.project}/analysis.log]")
                analyze_data(f"outputs/{args.project}/data.csv", f"outputs/{args.project}/pickles/spline.pkl", args.ref, args.project, args.protein, chosen_mutations)
                
            elif response == '3':
                print(f"Calculating Effects and Logging Output to [outputs/{args.project}/effect.log]")
                calculate_effect(args.project)
            elif response == '4':
                print("HELP MESSAGE")
            elif response == '5':
                print(f"Exiting Project [{args.project}]")
                exit(0)
            else:
                print("Invalid Request: please enter a number 1-5")
    else: # not interactive
        generate_data(args.ref, args.input, args.metric, args.project, args.protein, coeffs)
        analyze_data(f"outputs/{args.project}/data.csv", f"outputs/{args.project}/pickles/spline.pkl", args.ref, args.project, args.protein)
        calculate_effect(args.project)
        
    

if __name__ == "__main__":
    main()
