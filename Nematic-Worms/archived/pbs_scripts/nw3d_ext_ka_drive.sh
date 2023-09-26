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

datadir="/nfs/02/ksu0236/mvarga/worms/data/q3d_extensile"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
#############################
for ka in 0.5 2 4 6.5 9
do   
    for drive in 0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
    do
    	name="ka${ka}_drive${drive}"
	qsub -v NAME=$name,DATADIR=$datadir,DT=0.0025,NP=15,NSTEPS=2000000,FRAMERATE=20000,XBOX=260,YBOX=260,ZBOX=3,XDIM=20,YDIM=220,ZDIM=2,KBT=0.25,KA=${ka},DRIVE=${drive} ~/mvarga/worms/exp/nw3d_run.pbs
    done
done

