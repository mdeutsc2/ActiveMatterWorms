#!/bin/bash

#SBATCH --job-name=amatter_sim_chpl
#SBATCH --account=pgs0213               # account for mdeutsc2
#SBATCH --time=6:00:00                  # timeout for the batch job
#SBATCH --nodes=1                       # requesting number of nodes
#SBATCH --ntasks-per-node=16            # requesting number of cpus/node
#SBATCH --mail-type=END,FAIL            # send email to submitting user on job completion or failure

set -e
#loading modules
module load cmake/3.25.2
module load cuda/11.8.0

#executable section
source ~/chapel-1.31.0/util/setupchplenv.bash
make clean
make
make run
