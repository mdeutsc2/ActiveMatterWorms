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
xdim=
for np in 5 10 20 30 40
do   
    for drive in  0.2 0.4 0.6 0.8 1.0
    do
	if [ "$np" -eq "5" ]; then xdim=30; fi
	if [ "$np" -eq "10" ]; then xdim=15; fi
	if [ "$np" -eq "20" ]; then xdim=8; fi
	if [ "$np" -eq "30" ]; then xdim=5; fi
	if [ "$np" -eq "40" ]; then xdim=4; fi
	
	name="ramp_xlink_np${np}_drive${drive}"
	qsub -v NAME=$name,DATADIR=$datadir,DT=0.002,NP=$np,NSTEPS=2000000,FRAMERATE=20000,XBOX=140,YBOX=140,XDIM=$xdim,YDIM=60,ZDIM=3,XLINK=0.9,DRIVE=$drive,KBT=0.25,XSTART=200000,KX=10.0,LX=1.5 ~/mvarga/worms/exp/xlinkers_ramp.pbs
    done
done
