#!/bin/bash
#####################################################
### experiment studying flexible active filaments ###
### in a flexible encapsulation                   ###
#####################################################

#	Uses the -encapsilate2d flag to initiate
#   simulation in a finite flexible ring of particles.
#   Filaments are initiallized in concentric cirlces
#   facing radial (+/- r_hat).  Make sure box size is
#   larger enough that a circle of diameter = box_width / 2
#   will fit all filaments.  The encapsulation is held
#   fix until nsteps/5, when it begins to shrink to a
#   natural neighbor bond length of 0.35.

####################################################
### make directory for experiment data if needed ###
####################################################

datadir="/nfs/02/ksu0236/mvarga/worms/data/encap_kak2_2"
[ -d $datadir ] || mkdir $datadir
cd $datadir

##############################
# sub multipule jobs in loop #
#############################
for ka in 1 2 3 4 5 7 9 12 16 20
do
for k2_encap in 100 500 1000
do    
	name="encap_ka${ka}_k2_${k2_encap}"
	qsub -v NAME=$name,DATADIR=$datadir,DT=0.001,NP=15,NSTEPS=20000000,FRAMERATE=20000,XBOX=1500,YBOX=1500,DRIVE=2,XDIM=11,YDIM=74,KA=${ka},K2=${k2_encap} ~/mvarga/worms/exp/nw2d_encap_run.pbs
done
done
