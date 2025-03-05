
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




def fold(name, seq, partition):
    
    name = name.replace(":", "_").replace(">", "_").lower()
    if os.path.exists(f"/shared/25mdl4/af_output/{name}"):
        print(f"SKIPPING: [{name}] has already been folded")
        return True
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
        return True
    return False
        # job has finished
        # read the output and then do the comparison
        


def fold_all(data_file, partition):
    df = pd.read_csv("prion_uniprot.csv")[::-1]
    
    
    # Example function to be executed by threads

    
    num_workers = 3 if partition == "a" else 7
# Create a ThreadPoolExecutor with 4 threads
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        
        for index, row in df.iterrows():
            executor.submit(fold, row["ID"], row["Sequence"], partition)
    # tasks = [executor.submit(worker_function, i) for i in range(10)]

    
        
        



# def setup_sigint_handler(runtime_param):
#     def handle_sigint(signal_number, frame):
#         print("\nInterrupt Detected")
#         print("Canceling Job...")
        
#         result = subprocess.run(
#                     ["scancel", JOB_NUM],
#                     capture_output=True,
#                     text=True
#                 )
#         print("Script Output:", result.stdout)
#         print("Script Errors:", result.stderr)

#         exit(1)  # Exit gracefully
#     return handle_sigint
   
def main():
    # signal.signal(signal.SIGINT, setup_sigint_handler(JOB_NUM))
    
    # signal.signal(signal.SIGINT, setup_sigint_handler(JOB_NUM))
#     file_path = "reference.txt"

# # Open the file and read its content into a string
#     with open(file_path, "r") as file:
#         file_content = file.read()
        
#     file_content.strip()
#     fold("reference", file_content)
    if len (sys.argv) < 2:
        print("Usage python driver.py [a/b]")
        exit(1)
        
    assert sys.argv[1] in ["a", "b", "b1"]
    
    
    fold_all("data.json", sys.argv[1])
    
    
if __name__ == "__main__":
    main()
