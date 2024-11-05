#!/bin/bash
##############################
# sub multipule jobs in loop #
##############################

for YDIM in 10 20 30 40 50 60 70 80 90 100
do
#    for KBT in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0   
    for KBT in 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1.0 1.125 1.25
    do
	qsub -v kbt=${KBT},ka=5,np=10,xdim=15,ydim=${YDIM},zdim=3,xbox=130,ybox=130 kbt_density.pbs
    done
done
