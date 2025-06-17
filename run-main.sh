#!/bin/sh
#SBATCH -c 16             # Request 1 CPU core
#SBATCH -t 0-02:00          # Runtime in D-HH:MM, minimum of 10 mins
                                # (this requests 2 hours)
#SBATCH --partition=gpmoo-a          # Partition to submit to
#SBATCH --mem=80G           # Request 10G of memory
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written
                                # (%j inserts jobid)
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written
                                # (%j inserts jobid)
#SBATCH --gres=gpu:1        # Request two GPUs           

# Command you want to run on the cluster
# Notice, you must set-up testEval correctly as a conda virtual environment
# Calling this full path makes sure you are running the correct package versions

python main.py --project prionblosum2 --ref prion_ref.txt --input prion_input_uni.csv --metric blosum  --protein prion