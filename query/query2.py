# this is AL(ref, germline-mut) vs AL(ref, germline+mut)

import pandas as pd
import os
import time
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor



from mutgen import generate_brca_mutations, generate_prion_mutation
from alt_align import compute_align
from predict import calc_point_and_conf


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
        
def fold(name, seq, partition):
    
    if name == "":
        name = "prion_ref"
        
        
    name = name.replace(":", "_").replace(">", "_").replace(",", "-").lower()
    if os.path.exists(f"/shared/25mdl4/af_output/{name}"):
        print(f"SKIPPING: [{name}] has already been folded")
        return name
    else:
        print(f"Folding {name}")
    
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
        
    
    run_cmd = f"./run-alpha-{partition}.sh"
    
    result = subprocess.run(["sbatch", run_cmd, my_filename], capture_output=True, text=True)

    
    # print("Script Output:", result.stdout)
    # print("Script Errors:", result.stderr)
    # print("Return Code:", result.returncode)
    
    
    job_num = result.stdout.strip().split()[-1]
    # print(job_num)
    
    
    if check_job(job_num):
        
        end_time = time.time()
        print(f"{name} folded successfully. in {end_time-start_time}s.")
        return name
    return None
        # job has finished
        # read the output and then do the comparison
    
def fold_with_wo_and_score(with_data, wo_data, partition):
    # presumably either prev or post should already be folded so we can use 1 thread for this
    
    with_fold = fold(with_data[0], with_data[1], partition)
    wo_fold = fold(wo_data[0], wo_data[1], partition)
    
    assert with_fold is not None
    assert wo_fold is not None
    
    # FIXME make dynamic
    with_align = compute_align("prion_ref", with_fold)
    wo_align = compute_align("prion_ref", wo_fold)
    
    # align_score = with_align - wo_align
    
    # if align_score < 0:
    #     print(f"WARNING: negative score {align_score}")
    #     align_score = 0
    
    return with_align, wo_align
    
    
    
    
    

        
def query(mut_data, data, mutation_generator, num_individuals=100, partition="b"):
    
    mutation = mut_data.ID
    
    if num_individuals > len(data):
        print(f"Not Enough data for sample size {num_individuals}")
        exit(1)
    
    population = data.sample(n=num_individuals, replace=False, random_state=2)
    
    # population.to_csv(f'pre_treatment_{mutation.replace(":", "_").replace(">", "_")}.csv', index=False)
    
    post_treatment = []
    
    for index, row in population.iterrows():
        
        individual_muts = row["ID"].split(",")
        germline = row["ID"]
        
        if mutation in individual_muts:
            with_mut = germline
            individual_muts.remove(mutation) 
            wo_mut = "_".join(individual_muts)
            # uninduce mutation
        else:
            wo_mut = germline
            individual_muts.append(mutation)
            # FIXME: this is _ for prion (xnx) or , for BRCA (hgvs)
            with_mut = "_".join(individual_muts)
            # induce mutation
        
        # print(wo_mut, with_mut)
        with_data = mutation_generator([with_mut], [row["is_pathogenic"]])[0]
        wo_data = mutation_generator([wo_mut], [row["is_pathogenic"]])[0]
        
        
        
        # print("prev", prev_data)
        # print("post", post_data)
        # print("------------------")
        
        # TODO: right now we are gonna hang here -- need to parallelize
        with_align, wo_align = fold_with_wo_and_score(with_data, wo_data, partition)
        
        # in loop to see progress
        
        new_row = {
            "Old": row["Old"],
            "White": row["White"],
            "Unhealthy": row["Unhealthy"],
            "Align_Score_do_mut": with_align,
            "Align_Score_do_no_mut": wo_align,
            "Ground": row["is_pathogenic"],
        
        }
        
        post_treatment.append(new_row)
        
        post_treatment_df = pd.DataFrame(post_treatment)
        # post_treatment_df.to_csv(f'post_treatment_{mutation.replace(":", "_").replace(">", "_")}_{"p" if mut_data.is_pathogenic else "b"}.csv', index=False)
        post_treatment_df.to_csv(f'post_treatment_{mutation}_{"p" if mut_data.is_pathogenic else "b"}.csv', index=False)
        

def thread(mut, data, partition):
    
    query(mut, data, generate_prion_mutation, partition=partition)
    
    filename = f'post_treatment_{mut.ID}_{"p" if mut_data.is_pathogenic else "b"}.csv'
    
    # with open("prion_output.txt", "a") as file:
    #     print(f'{mut} [{mut.is_pathogenic}]: {calc_point_and_conf(filename, spline_filename)}', file=file)

         
        
def process_mutations(data_filename, spline_filename, partition="b"):
    data = pd.read_csv(data_filename)

    all_mutations = data[['ID', 'is_pathogenic']].drop_duplicates(subset=['ID'])
    
    # print(sum(all_mutations['is_pathogenic']))
    
    num_workers = 3 if partition == "a" else 8
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # for mut in all_mutations.itertuples(index=False):
        #     executor.submit(thread, mut, partition)
        # # print(mut)
        
        futures = [executor.submit(thread, mut, data, partition) for mut in all_mutations.itertuples(index=False)]
        
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