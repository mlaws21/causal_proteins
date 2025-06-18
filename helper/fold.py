
# we want to take in a sequence, format it properly, write the json to be read by alphafold, then read the output when complete then return the align score
# this will be fun

from datetime import datetime
import json
import subprocess
import time
import signal
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import sys

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
        # print("running")

# def cleanup():




def fold(name, seq, partition, log_fn):

    if os.path.exists(f"/shared/25mdl4/af_output/{name}"):
        log_fn(f"SKIPPING: [{name}] has already been folded")
        return True
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
        
    # TODO FIXME
    if node_has_free_gpu("gpmoo-a1") < 4:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'a'}.sh"
    elif node_has_free_gpu("gpmoo-b2") < 8:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'b2'}.sh"
    else:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'b1'}.sh"
        
    
    result = subprocess.run(["sbatch", run_cmd, my_filename], capture_output=True, text=True)

    # print("Script Output:", result.stdout)
    # print("Script Errors:", result.stderr)
    # print("Return Code:", result.returncode)
    
    job_num = result.stdout.strip().split()[-1]
    
    
    if check_job(job_num):
        
        end_time = time.time()
        log_fn(f"{name} folded successfully. in {end_time-start_time}s.")
        return True
    return False
        # job has finished
        # read the output and then do the comparison
        
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
def fold_all(ids, seqs, partition, protein_name, log_fn, num_workers=None):
    # Example function to be executed by threads
    
    log_fn("Beginning Fold Routine")
    
    if num_workers == None:
        num_workers = 3 if partition == "a" else 7
    log_fn(f"Using Partition {partition} with {num_workers} workers")
# Create a ThreadPoolExecutor with 4 threads
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []

        for i, (row_id, row_seq) in enumerate(zip(ids, seqs)):

            futures.append(executor.submit(fold, f"{protein_name}_{row_id}", row_seq, partition, log_fn))
    
        for future in as_completed(futures):
            try:
                # Wait for each future to finish and check for errors
                result = future.result()  # This will raise an exception if the thread encountered one
            except Exception as e:
                # Log the exception or print it as needed
                print(f"Error occurred in thread: {e}")
def main():

    if len (sys.argv) < 2:
        print("Usage python driver.py [a/b]")
        exit(1)
        
    assert sys.argv[1] in ["a", "b", "b1"]
    
    
    fold_all_from_datafile("data.json", sys.argv[1])
    
    
if __name__ == "__main__":
    main()
