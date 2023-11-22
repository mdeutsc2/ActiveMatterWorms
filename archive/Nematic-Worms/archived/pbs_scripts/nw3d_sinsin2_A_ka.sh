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

datadir="/nfs/02/ksu0236/mvarga/worms/data/sinsin2"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
#############################
#for A in 0 5 10 15 20 25
#do
#    for ka in 2.5 25
#    do    
#	for drive in 0 1 2
#	do
#	    name="sinsin2_A${A}_ka${ka}_drive${drive}"
#	    qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs
#    	done
#    done
#done

# comment to execute rerun section
#exit(0)

#################
# RERUN SECTION #
#################

#1
A=0
ka=2.5
drive=1
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#2
drive=2
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#3
ka=25
drive=0
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#4
drive=2
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#5
A=5
drive=0
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#6
drive=2
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#7
A=10
ka=2.5
drive=1
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#8
ka=25
drive=0
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#9
A=15
ka=2.5
drive=1
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#10
ka=25
drive=0
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#11
A=20
ka=2.5
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#12
ka=25
drive=1
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#13
A=25
ka=2.5
drive=0
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#14
drive=2
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs

#15
ka=25
drive=1
name="sinsin2_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=5000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=${ka} ~/mvarga/worms/exp/nw3d_sinsin2_run.pbs
