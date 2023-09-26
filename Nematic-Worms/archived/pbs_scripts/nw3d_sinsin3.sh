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

datadir="/nfs/02/ksu0236/mvarga/worms/data/sinsin3"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
#############################
A=25
#for ka in 1 10 30
#do    
#    for drive in 0 2.5 5 7.5 10
#    do
#        name="sinsin3_A${A}_ka${ka}_drive${drive}"
#        qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs
#    done
#done

# comment to execute rerun section
#exit 0

#################
# RERUN SECTION #
#################

#1
ka=10
drive=10
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs

#2
drive=5
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs


#3
drive=7.5
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs


#4
ka=1
drive=10
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs


#5
drive=5
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs


#6
drive=7.5
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs


#7
ka=30
drive=10
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs


#8
drive=7.5
name="sinsin3_A${A}_ka${ka}_drive${drive}"
qsub -v NAME=$name,DATADIR=$datadir,WGAM=10,AMP=$A,DT=0.002,DRIVE=${drive},NP=10,NSTEPS=7000000,FRAMERATE=10000,XBOX=200,YBOX=200,XDIM=23,YDIM=180,KA=25 ~/mvarga/worms/exp/nw3d_sinsin3_run.pbs
