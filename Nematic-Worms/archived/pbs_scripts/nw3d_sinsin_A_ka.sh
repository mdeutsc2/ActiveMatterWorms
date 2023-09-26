#!/bin/bash
#####################################################
### experiment studying flexible active filaments ###
### in a flexible encapsulation                   ###
#####################################################

#	Uses the (-A|-landscale) flag to set a the
#   amplitude of a Sin(x)*Sin(y) curvature on simulation
#   domain.  Surface parameterized by:
#
#	r(x,y) = { x, y, z = f(x,y) }
#
#   where
#
#	f(x,y) = A sin(qx*x) sin(qy*y).
#
#   Tangents of this curve are used to project forces
#   in 3D onto surface where particle positions are
#   updated as normal once this constraint is imposed.

####################################################
### make directory for experiment data if needed ###
####################################################

datadir="/nfs/02/ksu0236/mvarga/worms/data/sinsin"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
#############################
for A in 0 10 20 30 40 50
do
    for ka in 5 25
    do    
	for drive in 0 1 2
	do
	    name="sinsin_A${A}_ka${ka}_drive${drive}"
	    qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.0025,DRIVE=${drive},NP=10,NSTEPS=2000000,FRAMERATE=2000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin_run.pbs
    	done
    done
done
