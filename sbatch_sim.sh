#!/bin/bash

#SBATCH --job-name=amatter_sim_chpl
#SBATCH --account=pgs0213               # account for mdeutsc2
#SBATCH --time=72:00:00                  # timeout for the batch job
#SBATCH --nodes=1                       # requesting number of nodes
#SBATCH --ntasks-per-node=16            # requesting number of cpus/node
#SBATCH --mail-type=END,FAIL            # send email to submitting user on job completion or failure

set -e
#loading modules
module load cmake/3.25.2
#module load cuda/11.8.0
module load gnu/12.3.0
#executable section
lscpu
source /users/PGS0213/mdeutsc2/chapel-1.33.0_cpu/util/setchplenv.bash
make clean
make
./runsim.sh
