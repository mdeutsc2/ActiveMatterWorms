#!/bin/bash
#####################################################
### experiment studying isolated active filaments ###
#####################################################

#	Uses -noint flag in .pbs script to remove
#   inter-filament interactions.  Many particles can
#   be simulated but they will not interact sterically
#   with each other.

###################################################
### make director for experiment data if needed ###
###################################################

datadir="/nfs/02/ksu0236/mvarga/worms/data/singles"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
#############################
for ka in 1 2 3 5 6 7 8 9 10 15 20 30
do   
    for drive in 0 2.5 5 7.5 10 15 20
    do
	# multiple sims run in one job for all temperatures needed
	#for kbt in 0 1 2.5 5 7.5 10 20
	name="ka${ka}_drive${drive}" #_kbt${kbt} added in job
	qsub -v NAME=$name,DATADIR=$datadir,DT=0.005,NP=15,NSTEPS=1500000,FRAMERATE=200,XBOX=200,YBOX=200,XDIM=5,YDIM=5,KA=${ka},DRIVE=${drive} ~/mvarga/worms/exp/nw2d_single_run.pbs
    	
    done
done

