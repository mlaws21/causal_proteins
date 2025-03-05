
# - mutation should be in XNY form 
# - data needs to include sequences, mutations
# - dr_curve_func should be the dose response curve -- takes in values [0, inf)
#   that represent align scores, then outputs number [0, 1] as p(cancer)


# have to recompute the sequences


def fold(name, seq, partition):
    
    name = name.replace(":", "_").replace(">", "_").replace(",", "-").lower()
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
    


        
def query(mutation, data, dr_curve_func, num_individuals=100):
    
    if num_individuals > len(data):
        print(f"Not Enough data for sample size {num_individuals}")
        exit(1)
    
    population = data.sample(n=num_individuals, replace=False, random_state=2)
    
    
    for index, row in df.iterrows():
        if row["ID"] == mutation:
            # uninduce mutation
            name = name.replace(":", "_").replace(">", "_").replace(",", "-").lower()
            
        else:
            # induce mutation
    
    
    

def main():
    pass


if __name__ == "__main__":
    main()