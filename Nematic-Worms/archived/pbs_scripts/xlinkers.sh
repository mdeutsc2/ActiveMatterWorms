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

datadir="/nfs/02/ksu0236/mvarga/worms/data/xlinkers"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
##############################

for xlink in 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5
do   
    for drive in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    do
	name="xlink${xlink}_drive${drive}"
	qsub -v NAME=$name,DATADIR=$datadir,DT=0.002,NSTEPS=1000000,FRAMERATE=20000,XBOX=130,YBOX=130,XDIM=10,YDIM=60,ZDIM=3,XLINK=$xlink,DRIVE=$drive,KBT=0.1,KX=10.0,LX=1.122 ~/mvarga/worms/exp/xlinkers.pbs
    done
done
