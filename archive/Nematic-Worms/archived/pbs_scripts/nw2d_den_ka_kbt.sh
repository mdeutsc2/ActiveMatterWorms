#!/bin/bash
#################################################
### experiment to study cross-linking effects ###
#################################################

#	Uses the -L flag on the 'nw3d' simulation
#    invokation to set the cross-linker density.
#    This is an attemt to model C. Marchetti's
#    work on cross-linked microtubules on a 2D
#    interface.  Minimal noise used.

###################################################
### make director for experiment data if needed ###
###################################################

datadir="/nfs/02/ksu0236/mvarga/worms/data/2d_den_ka_kbt"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
##############################
xdim=
for kbt in 0.05 0.4 0.8
do
    for ka in 0.5 2 4 6.5 9
    do   
    	for ydim in 40 80 120 160 200
    	do
	    name="kbt${kbt}_ka${ka}_ydim${ydim}"
	    qsub -v NAME=$name,DATADIR=$datadir,DT=0.001,NP=10,NSTEPS=5000000,FRAMERATE=50000,XBOX=260,YBOX=260,XDIM=30,YDIM=${ydim},KBT=${kbt},KA=${ka} ~/mvarga/worms/exp/nw2d_run.pbs
	done
    done
done
