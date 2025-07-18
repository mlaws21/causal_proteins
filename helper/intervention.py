# this is AL(ref, germline-mut) vs AL(ref, germline+mut)

import pandas as pd
import os
import time
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import itertools
import traceback

# - mutation should be in XNY form 
# - data needs to include sequences, mutations
# - dr_curve_func should be the dose response curve -- takes in values [0, inf)
#   that represent align scores, then outputs number [0, 1] as p(cancer)


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

def node_has_free_gpu(node):
    result = subprocess.check_output([
        "squeue",
        "--noheader",
        "--format=%N"
    ], text=True)

    used = 0
    for line in result.strip().split("\n"):
        if line == node:
            used += 1
    return used
      
def get_run_cmd():
    if node_has_free_gpu("gpmoo-a1") < 4:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'a'}.sh"
    elif node_has_free_gpu("gpmoo-b1") < 7:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'b1'}.sh"
    else:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'b2'}.sh"
    
    return run_cmd

def fold(name, seq, protein_name, log_fn):

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
        
    run_cmd = get_run_cmd()

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

        log_fn(f"{name} folded successfully. in {end_time-start_time}s.")
        return name
    return "BAD"
        # job has finished
        # read the output and then do the comparison
        
def compute_align(ref_id, seq_id, protein_name, alignment_function, alignment_args, log_fn):
    
    out_root = "/shared/25mdl4/af_output/"
    

    folder1 = os.path.join(out_root, f"{seq_id}")
    folder2 = os.path.join(out_root, f"{ref_id}")
    
    if not os.path.isdir(folder1):
        log_fn(f"WARNING: {seq_id} invalid path -- skipping" )
        return None
        
    cif_file1 = f'{folder1}/seed-2_sample-0/model.cif'
    cif_file2 = f'{folder2}/seed-2_sample-0/model.cif'
    # FIXME
    # json_file1 = f'{folder1}/seed-2_sample-0/confidences.json'  # Replace with your first JSON file
    # json_file2 = f'{folder2}/seed-2_sample-0/confidences.json'
    align_score = alignment_function(cif_file1, cif_file2, *alignment_args)
    
    if align_score is None:
        log_fn(f"WARNING: nonsense mutation -- skipping")
        return None
    
    return align_score
    
    
def fold_with_wo_and_score(with_data, wo_data, protein_name, alignment_function, alignment_args, log_fn):
    # presumably either prev or post should already be folded so we can use 1 thread for this
    
    
    with_fold = fold(f"{protein_name}_{with_data[0]}", with_data[1], protein_name, log_fn)
    wo_fold = fold(f"{protein_name}_{wo_data[0]}", wo_data[1], protein_name, log_fn)
    
    assert with_fold is not None
    assert wo_fold is not None
    
    ref_name = f"{protein_name}_ref"
    with_align = compute_align(ref_name, with_fold, protein_name, alignment_function, alignment_args, log_fn)
    wo_align = compute_align(ref_name, wo_fold, protein_name, alignment_function, alignment_args, log_fn)
    
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

def query(mut_data, data, ref_seq, num_individuals, project_name, protein_name, alignment_function, alignment_args, log_fn):
    # mut_data is mut_data.ID, mut_data.Ground
    mutation = mut_data.ID
    
    if num_individuals is None:
        num_individuals = len(data)
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
            
            
        with_data = ("_".join(with_mut), mutate_protein(ref_seq, with_mut, log_fn))
        wo_data = ("_".join(wo_mut), mutate_protein(ref_seq, wo_mut, log_fn))
        
        # with_data = mutation_generator([with_mut], [row["Ground"]])[0]
        # wo_data = mutation_generator([wo_mut], [row["Ground"]])[0]
        
        
        
        # print("prev", prev_data)
        # print("post", post_data)
        # print("------------------")
        
        # TODO: right now we are gonna hang here -- need to parallelize
        with_align, wo_align = fold_with_wo_and_score(with_data, wo_data, protein_name, alignment_function, alignment_args, log_fn)
        
        # in loop to see progress
        
        # new_row = {
        #     "Age": row["Age"],
        #     "Race": row["Race"],
        #     "Lifestyle": row["Lifestyle"],
        #     "Sequence_Score": row["Sequence_Score"],
        #     "Align_Score_do_mut": with_align,
        #     "Align_Score_do_no_mut": wo_align,
        #     "Ground": row["Ground"],
        
        # }
        # Just pull the existing covariates from the row, and add the new align scores
        new_row = {
            "Align_Score_do_mut": with_align,
            "Align_Score_do_no_mut": wo_align,
            "Ground": row["Ground"],
        }
        
        post_treatment.append(new_row)
        
        post_treatment_df = pd.DataFrame(post_treatment)
        pt_filename = f'outputs/{project_name}/intervention_data/{mutation}_{"p" if mut_data.Ground else "b"}.csv'
        
        log_fn(f"Post Treatment Data Created and saved to {pt_filename}")
        # post_treatment_df.to_csv(f'post_treatment_{mutation.replace(":", "_").replace(">", "_")}_{"p" if mut_data.Ground else "b"}.csv', index=False)
        post_treatment_df.to_csv(pt_filename, index=False)

         
        
def process_mutations(data_filename, spline_filename, ref_seq, subset, project_name, protein_name, alignment_function, alignment_args, log_fn, chosen_mutations=None, num_workers=None):
    data = pd.read_csv(data_filename)

    if chosen_mutations is None:
        all_mutations = data[['ID', 'Ground']].drop_duplicates(subset=['ID'])
    else:
        # Only process the specific IDs provided in chosen_mutations
        all_mutations = data[['ID', 'Ground']].drop_duplicates(subset=['ID'])
        all_mutations = all_mutations[all_mutations['ID'].isin(chosen_mutations)]
    
    log_fn(f"Processing: {all_mutations}")
    # print(sum(all_mutations['Ground']))
    if num_workers is None:
        num_workers = 1
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for mut in all_mutations.itertuples(index=False):
            # time.sleep(0.5)
            futures.append(executor.submit(query, mut, data, ref_seq, subset, project_name, protein_name, alignment_function, alignment_args, log_fn))
        
        for future in as_completed(futures):
            try:
                # Wait for each future to finish and check for errors
                result = future.result()  # This will raise an exception if the thread encountered one
            except Exception as e:
                # Log the exception or print it as needed
                print(f"Error occurred in thread: {e}")
                print(traceback.format_exc())
                


def process_mutations_serial(data_filename, spline_filename, ref_seq, subset, project_name, protein_name, alignment_function, log_fn, num_workers=None):
    data = pd.read_csv(data_filename)

    all_mutations = data[['ID', 'Ground']].drop_duplicates(subset=['ID'])
    

    for mut in all_mutations.itertuples(index=False):
        # time.sleep(0.5)
        query(mut, data, ref_seq, subset, project_name, protein_name, alignment_function, log_fn)
    


def main():
    pass


if __name__ == "__main__":
    main()