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

for kbt in 0.0 0.25 0.5 0.75 1.0
do   
    for drive in  0.2 0.4 0.6 0.8 1.0
    do
	name="ramp_xlink_kbt{kbt}_drive${drive}"
	qsub -v NAME=$name,DATADIR=$datadir,DT=0.002,NSTEPS=2000000,FRAMERATE=20000,XBOX=130,YBOX=130,XDIM=10,YDIM=60,ZDIM=3,XLINK=0.9,DRIVE=$drive,KBT=$kbt,KX=10.0,LX=1.122 ~/mvarga/worms/exp/xlinkers_ramp.pbs
    done
done
