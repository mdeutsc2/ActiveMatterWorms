#!/bin/bash

module load cuda
SRC_DIR=../simulation
HERE=$(pwd)
cd $SRC_DIR && make && cd $HERE

for activity in 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
do
    sbatch --export=NW_ARGS="-o=nw_a${activity}.xyz -activity=${activity}" $SRC_DIR/slurm_job.sh
done