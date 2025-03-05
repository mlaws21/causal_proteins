#!/bin/sh
#SBATCH -c 8               # Request 1 CPU core
#SBATCH -t 0-02:00          # Runtime in D-HH:MM, minimum of 10 mins
                                # (this requests 2 hours)
#SBATCH --partition=gpmoo-b          # Partition to submit to
#SBATCH --mem=80G           # Request 10G of memory
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written
                                # (%j inserts jobid)
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written
                                # (%j inserts jobid)
#SBATCH --gres=gpu:1        # Request two GPUs           

# Command you want to run on the cluster
# Notice, you must set-up testEval correctly as a conda virtual environment
# Calling this full path makes sure you are running the correct package versions

singularity exec \
    --nv \
    --bind $HOME/af_input:/25mdl4/af_input \
    --bind $HOME/af_output:/25mdl4/af_output \
    --bind $HOME/model_params/:/25mdl4/models \
    --bind $HOME/public_databases:/25mdl4/public_databases \
    $HOME/alphafold3/alphafold3.sif \
    python $HOME/alphafold3/run_alphafold.py \
    --json_path=/25mdl4/af_input/$1 \
    --model_dir=/25mdl4/models \
    --db_dir=/25mdl4/public_databases \
    --output_dir=/25mdl4/af_output
