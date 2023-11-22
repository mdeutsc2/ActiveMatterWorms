#!/bin/bash
#########################################################
### experiment studyindg mixing dead and living cells ###
#########################################################


### make directory for experiment data if needed 

datadir="/nfs/02/ksu0236/mvarga/worms/data/Hoeger4"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
##############################
drive=10
for rev in 1 5 10
do
for alive in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 1.0
do    
	name="deadrev_alive${alive}_drive${drive}_taorev${rev}_kbt0.005_gamma20"
	qsub -v NAME=$name,DATADIR=$datadir,ALIVE=${alive},TAOREV=${rev},DT=0.001,NP=6,NSTEPS=40000000,FRAMERATE=20000,XBOX=250,YBOX=250,DRIVE=${drive},XDIM=42,YDIM=220,KA=15,KBT=0.005,WEPS=0.75,WGAM=20 ~/mvarga/worms/exp/nw2d_deadcells_run.pbs
done
done
