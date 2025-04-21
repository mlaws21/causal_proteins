# this is AL(ref, germline-mut) vs AL(ref, germline+mut)

import pandas as pd
import os
import time
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import itertools



# from mutgen import generate_brca_mutations, generate_prion_mutation
# FIXME this is direct pathing call from main
from alphafold_driver.alt_align import tm_align_rmsd
# from predict import calc_point_and_conf


# - mutation should be in XNY form 
# - data needs to include sequences, mutations
# - dr_curve_func should be the dose response curve -- takes in values [0, inf)
#   that represent align scores, then outputs number [0, 1] as p(cancer)


# have to recompute the sequences

def check_job(job_id):
    while True:
        try:
            # Check if the job is in the queue
            result = subprocess.run(
                ["squeue", "-j", job_id],
                capture_output=True,
                text=True
            )
            # If the job ID is in the output, the job is not done
            if job_id not in result.stdout:
                return True
        except Exception as e:
            print(f"Squeue failed: {e}")
            exit(1)
    
        time.sleep(5)
        
def fold(name, seq, partition, protein_name, log_fn):

    if name == f"{protein_name}_":
        name = f"{protein_name}_ref"
    
    muts = name.split("_")[1:]
    perms = list(itertools.permutations(muts))
    
    
    for p in perms:
        new_name = f"{protein_name}_{'_'.join(p)}"
        if os.path.exists(f"/shared/25mdl4/af_output/{new_name}"):
            log_fn(f"SKIPPING: [{new_name}] has already been folded")
            return new_name
    else:
        log_fn(f"Folding {name}")

    start_time = time.time()
    # datet = datetime.now().strftime("%m.%d_%H.%M.%S")
    
    inp_data = {
        "name": name,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq
                }
            }
        ],
        "modelSeeds": [2],
        "dialect": "alphafold3",
        "version": 1
    }

    my_filename = f"{name}.json"
    
    with open("/shared/25mdl4/af_input/" + my_filename, "w") as json_file:
        json.dump(inp_data, json_file, indent=2)
        
    
    run_cmd = f"/shared/25mdl4/thesis/alphafold_driver/run-alpha-{partition}.sh"
    
    result = subprocess.run(["sbatch", run_cmd, my_filename], capture_output=True, text=True)

    
    # print("Script Output:", result.stdout)
    if result.stderr:
        log_fn(f"Slurm Errors: {result.stderr}")
    # print("Return Code:", result.returncode)
    
    try:
        job_num = result.stdout.strip().split()[-1]
    except Exception as e:
        print(job_num)
        print(name)
        exit(1)
    # print(job_num)
    
    
    if check_job(job_num):
        
        end_time = time.time()
        with log_lock:
            log(f"{name} folded successfully. in {end_time-start_time}s.")
        return name
    return "BAD"
        # job has finished
        # read the output and then do the comparison
        
def compute_align(ref_id, seq_id, protein_name, log_fn):
    
    out_root = "/shared/25mdl4/af_output/"
    
    # ID,Old,White,Unhealthy,Align Score,Cancer,Sequence,Ground
    # ids = data["ID"]
    

    
    # name = name.replace(":", "_").replace(">", "_").lower()
    

    folder1 = os.path.join(out_root, f"{seq_id}")
    folder2 = os.path.join(out_root, f"{ref_id}")
    
    if not os.path.isdir(folder1):
        log_fn(f"WARNING: {seq_id} invalid path -- skipping" )
        return None
        
    cif_file1 = f'{folder1}/seed-2_sample-0/model.cif'
    cif_file2 = f'{folder2}/seed-2_sample-0/model.cif'
    json_file1 = f'{folder1}/seed-2_sample-0/confidences.json'  # Replace with your first JSON file
    json_file2 = f'{folder2}/seed-2_sample-0/confidences.json'
    align_score = tm_align_rmsd(cif_file1, cif_file2, json_file1, json_file2)
    
    if align_score is None:
        log_fn(f"WARNING: nonsense mutation -- skipping")
        return None
    
    return align_score
    
    
def fold_with_wo_and_score(with_data, wo_data, protein_name, partition, log_fn):
    # presumably either prev or post should already be folded so we can use 1 thread for this
    
    
    with_fold = fold(f"{protein_name}_{with_data[0]}", with_data[1], partition, protein_name, log_fn)
    wo_fold = fold(f"{protein_name}_{wo_data[0]}", wo_data[1], partition, protein_name, log_fn)
    
    assert with_fold is not None
    assert wo_fold is not None
    
    ref_name = f"{protein_name}_ref"
    with_align = compute_align(ref_name, with_fold, protein_name, log_fn)
    wo_align = compute_align(ref_name, wo_fold, protein_name, log_fn)
    
    # align_score = with_align - wo_align
    
    # if align_score < 0:
    #     print(f"WARNING: negative score {align_score}")
    #     align_score = 0
    
    return with_align, wo_align
    
    
    
    
def mutate_protein(ref, mutations, log_fn):
    # mutations should be a list of XNX format mutations
    # you can also give mutation groups which are XNX mutations 
    # separate by "_": XNX_XNX...
    # stop codon is denoted "-"
    
    ref = ref.lower()
    
    for mut in mutations:
        
        outchar = mut[0]
        inchar = mut[-1]

        position = int(mut[1:-1])
        # Char doesnt match
        # assert ref[position-1] == outchar
        if ref[position-1] != outchar:
            log_fn("ERROR")
            log_fn(mut)
            log_fn(mutations)
            log_fn(position)
            log_fn(f"{ref[position-1]},{outchar}")
            exit(1)
        
        # this is a truncation
        if inchar == "-":
            ref = ref[:position-1] 
        else:
            ref = ref[:position-1] + inchar + ref[position:]
    

    return ref.upper()

def match_mut_location(new_muts, germline):
    
    matches = set()
    for nmut in new_muts:
        nmut_loc = nmut[1:-1]
    
        for i in germline:
            germ_mut_loc = i[1:-1]
            if nmut_loc == germ_mut_loc:
                matches.add(i)
                break
    return matches

def query(mut_data, data, ref_seq, num_individuals, partition, project_name, protein_name, log_fn):
    # mut_data is mut_data.ID, mut_data.Ground
    mutation = mut_data.ID
    

    # we take a random sample to make shit faster
    if num_individuals > len(data):
        print(f"Not Enough data for sample size {num_individuals}")
        exit(1)
    
    # this is the full data
    population = data.sample(n=num_individuals, replace=False, random_state=2)
    
    post_treatment = []
    
    if mutation == "ref":
        new_mutations = set()
    else:
        new_mutations = set(mutation.split("_"))
    for index, row in population.iterrows():

        # TODO if the new mutation is in the same location --> overwrite the old one
        
        # we want to check if the induced mutation is equal to the old mutation

        # individual_muts = set(row["ID"].split("_")
        
        if row["ID"] == "ref":
            germline = set()
        else:
            germline = set(row["ID"].split("_"))
        
        
        if new_mutations.issubset(germline):
            with_mut = germline
            wo_mut = germline.difference(new_mutations)
        else:
            position_matches = match_mut_location(new_mutations, germline)
            temp_with = germline.difference(position_matches)
            with_mut = temp_with.union(new_mutations)
            wo_mut = germline
            
            
        # print(with_mut, wo_mut)
        
        # for nmut in new_mutations:
        #     if nmut in individual_muts:
        #         with_mut = set(germline)
        #         individual_muts.remove(nmut) 
        #         wo_mut = set(individual_muts)
        #         # uninduce mutation
        #     elif len((matches := match_mut_location(nmut, germline))) == 0:
        #         wo_mut = set(germline)
        #         individual_muts.append(nmut)
        #         for i in matches:
        #             individual_muts.remove(i) 
        #         with_mut = set(individual_muts)    
                
        #     else:
        #         wo_mut = set(germline)
        #         individual_muts.append(nmut)
        #         # FIXME: this is _ for prion (xnx) or , for BRCA (hgvs)
        #         with_mut = set(individual_muts)
        #         # induce mutation
        
        # print(wo_mut, with_mut)
        with_data = ("_".join(with_mut), mutate_protein(ref_seq, with_mut, log_fn))
        wo_data = ("_".join(wo_mut), mutate_protein(ref_seq, wo_mut, log_fn))
        
        # with_data = mutation_generator([with_mut], [row["Ground"]])[0]
        # wo_data = mutation_generator([wo_mut], [row["Ground"]])[0]
        
        
        
        # print("prev", prev_data)
        # print("post", post_data)
        # print("------------------")
        
        # TODO: right now we are gonna hang here -- need to parallelize
        with_align, wo_align = fold_with_wo_and_score(with_data, wo_data, protein_name, partition, log_fn)
        
        # in loop to see progress
        
        new_row = {
            "Age": row["Age"],
            "Race": row["Race"],
            "Lifestyle": row["Lifestyle"],
            "Sequence_Score": row["Sequence_Score"],
            "Align_Score_do_mut": with_align,
            "Align_Score_do_no_mut": wo_align,
            "Ground": row["Ground"],
        
        }
        
        post_treatment.append(new_row)
        
        post_treatment_df = pd.DataFrame(post_treatment)
        pt_filename = f'intervention_data/{project_name}/post_treatment_{mutation}_{"p" if mut_data.Ground else "b"}.csv'
        
        log_fn(f"Post Treatment Data Created and saved to {pt_filename}")
        # post_treatment_df.to_csv(f'post_treatment_{mutation.replace(":", "_").replace(">", "_")}_{"p" if mut_data.Ground else "b"}.csv', index=False)
        post_treatment_df.to_csv(pt_filename, index=False)
        
        # return post_treatment_df
        

def thread(mut, data, ref_seq, partition, subset, log_fn):
    
    query(mut, data, ref_seq, subset, partition=partition)
    
    filename = f'post_treatment_{mut.ID}_{"p" if mut_data.Ground else "b"}.csv'
    # log_fn(f'Saved data to post_treatment_{mut.ID}_{"p" if mut_data.Ground else "b"}.csv')
    # with open("prion_output.txt", "a") as file:
    #     print(f'{mut} [{mut.Ground}]: {calc_point_and_conf(filename, spline_filename)}', file=file)

         
        
def process_mutations(data_filename, spline_filename, ref_seq, subset, project_name, protein_name, log_fn, partition, num_workers=None):
    data = pd.read_csv(data_filename)

    all_mutations = data[['ID', 'Ground']].drop_duplicates(subset=['ID'])
    
    # print(sum(all_mutations['Ground']))
    if num_workers is None:
        num_workers = 3 if partition == "a" else 8
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # for mut in all_mutations.itertuples(index=False):
        #     executor.submit(thread, mut, partition)
        # # print(mut)
        
        futures = [executor.submit(query, mut, data, ref_seq, subset, partition, project_name, protein_name, log_fn) for mut in all_mutations.itertuples(index=False)]
        
        # Wait for all futures to complete
        for future in futures:
            future.result()  # This ensures exceptions are raised if any thread fails


        
        

        
        

    

def main():
    
    process_mutations("prion_data.csv", "prion_spline.pkl")
    # # mutation = "chr17:43063882:C>T"
    # if len(sys.argv) < 2:
    #     print("Usage: python query [HGVS mutation]")
    
    # mutation = sys.argv[1]
    # print(mutation)
    
    # query(mutation, data, None, partition="a")
    # # print(generate_brca_mutations([""], ["path"]))
    
    
    # # print(generate_brca_mutations(["chr17:43094113:T>A,chr17:43091792:C>T,chr17:43094495:G>A"], ["Pathogenic"]))
    


if __name__ == "__main__":
    main()