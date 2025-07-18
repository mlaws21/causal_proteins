
# we want to take in a sequence, format it properly, write the json to be read by alphafold, then read the output when complete then return the align score
# this will be fun

from datetime import datetime
import json
import subprocess
import time
import signal
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
import traceback

def check_job(job_id):
    # print(f"Checking job {job_id}")
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

def choose_run_cmd():
    # TODO FIXME
    
    if node_has_free_gpu("gpmoo-a1") < 4:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'a'}.sh"
    elif node_has_free_gpu("gpmoo-b1") < 7:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-{'b1'}.sh"
    else:
        run_cmd = f"/shared/25mdl4/causal_proteins/helper/shell/run-alpha-b2.sh"

    return run_cmd

def fold(name, seq, log_fn):

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
    with open("/shared/25mdl4/af_input/" + my_filename, "w+") as json_file:
        json.dump(inp_data, json_file, indent=2)
    
    run_cmd = choose_run_cmd()
        
    
    result = subprocess.run(["sbatch", run_cmd, my_filename], capture_output=True, text=True)

    # print("Script Output:", result.stdout)
    # print("Script Errors:", result.stderr)
    # print("Return Code:", result.returncode)
    
    job_num = result.stdout.strip().split()[-1]
    
    # print(f"Job Number: {job_num}")
    
    
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
def fold_all(ids, seqs, protein_name, log_fn, num_workers=None):
    # Example function to be executed by threads
    
    log_fn("Beginning Fold Routine")
    
    if num_workers == None:
        num_workers = 1
    log_fn(f"Using {num_workers} workers")   
# Create a ThreadPoolExecutor with 4 threads
    started = set()
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []

        for i, (row_id, row_seq) in enumerate(zip(ids, seqs)):
            if f"{protein_name}_{row_id}" not in started:
                log_fn(f"SUBMITTING: [{protein_name}_{row_id}]")
                futures.append(executor.submit(fold, f"{protein_name}_{row_id}", row_seq, log_fn))
                started.add(f"{protein_name}_{row_id}")
            else:
                log_fn(f"SKIPPING: [{protein_name}_{row_id}] has already been submitted")
    
        for future in as_completed(futures):
            try:
                # Wait for each future to finish and check for errors
                result = future.result()  # This will raise an exception if the thread encountered one
            except Exception as e:
                # Log the exception or print it as needed
                print(f"Error occurred in thread: {e}")
                print(traceback.format_exc())
def main():

    fold("idh1_ref", "MSKKISGGSVVEMQGDEMTRIIWELIKEKLIFPYVELDLHSYDLGIENRDATNDQVTKDAAEAIKKHNVGVKCATITPDEKRVEEFKLKQMWKSPNGTIRNILGGTVFREAIICKNIPRLVSGWVKPIIIGRHAYGDQYRATDFVVPGPGKVEITYTPSDGTQKVTYLVHNFEEGGGVAMGMYNQDKSIEDFAHSSFQMALSKGWPLYLSTKNTILKKYDGRFKDIFQEIYDKQYKSQFEAQKIWYEHRLIDDMVAQAMKSEGGFIWACKNYDGDVQSDSVAQGYGSLGMMTSVLVCPDGKTVEAEAAHGTVTRHYRMYQKGQETSTNPIASIFAWTRGLAHRAKLDNNKELAFFANALEEVSIETIEAGFMTKDLAACIKGLPNVQRSDYLNTFEFMDKLGENLKIKLAQAKL", print)
    
if __name__ == "__main__":
    main()
