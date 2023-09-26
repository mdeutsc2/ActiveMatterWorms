#!/bin/bash
##############################
# sub multipule jobs in loop #
##############################

for YDIM in 10 20 30 40 50 60 70 80 90 100 110 120
do
    for KBT in 1.5
    do
	qsub -v kbt=${KBT},ka=5,np=10,xdim=15,ydim=${YDIM},zdim=3,xbox=130,ybox=130 kbt_density.pbs
    done
done
