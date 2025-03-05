import pandas as pd
import os
import time
import json
import subprocess
import sys


from mutgen import generate_brca_mutations
from align import compute_align


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
        name = "reference"
        
        
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
    
def fold_prev_post_and_align(prev_data, post_data, partition):
    # presumably either prev or post should already be folded so we can use 1 thread for this
    
    prev_fold = fold(prev_data[0], prev_data[1], partition)
    post_fold = fold(post_data[0], post_data[1], partition)
    
    assert prev_fold is not None
    assert post_fold is not None
    
    align_score = compute_align(prev_fold, post_fold)
    
    return align_score
    
    
    
    
    

        
def query(mutation, data, dr_curve_func, num_individuals=100, partition="b"):
    
    if num_individuals > len(data):
        print(f"Not Enough data for sample size {num_individuals}")
        exit(1)
    
    population = data.sample(n=num_individuals, replace=False, random_state=2)
    
    population.to_csv(f'pre_treatment_{mutation.replace(":", "_").replace(">", "_")}.csv', index=False)
    
    post_treatment = []
    
    for index, row in population.iterrows():
        
        individual_muts = row["ID"].split(",")
        prev = row["ID"]
        
        if mutation in individual_muts:
            # uninduce mutation
            individual_muts.remove(mutation) 
        else:
            # induce mutation
            individual_muts.append(mutation)

        post = ",".join(individual_muts)
        

        prev_data = generate_brca_mutations([prev], [row["is_pathogenic"]])[0]
        post_data = generate_brca_mutations([post], [row["is_pathogenic"]])[0]
        
        
        
        # print("prev", prev_data)
        # print("post", post_data)
        # print("------------------")
        
        # TODO: right now we are gonna hang here -- need to parallelize
        align_score = fold_prev_post_and_align(prev_data, post_data, partition)
        
        # in loop to see progress
        
        new_row = {
            "Old": row["Old"],
            "White": row["White"],
            "Unhealthy": row["Unhealthy"],
            "Align_Score": align_score,
            "Ground": row["is_pathogenic"],
        
        }
        
        post_treatment.append(new_row)
        
        post_treatment_df = pd.DataFrame(post_treatment)
        post_treatment_df.to_csv(f'post_treatment_{mutation.replace(":", "_").replace(">", "_")}.csv', index=False)
        
        

    

def main():
    data = pd.read_csv("final_full.csv")
    # mutation = "chr17:43063882:C>T"
    if len(sys.argv) < 2:
        print("Usage: python query [HGVS mutation]")
    
    mutation = sys.argv[1]
    print(mutation)
    
    query(mutation, data, None, partition="a")
    # print(generate_brca_mutations([""], ["path"]))
    
    
    # print(generate_brca_mutations(["chr17:43094113:T>A,chr17:43091792:C>T,chr17:43094495:G>A"], ["Pathogenic"]))
    


if __name__ == "__main__":
    main()