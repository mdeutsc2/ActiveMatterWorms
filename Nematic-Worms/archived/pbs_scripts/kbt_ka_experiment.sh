####################################
### linkage to folders and files ###
####################################

binfolder="~/mvarga/worms/bin"
datafolder="~/mvarga/worms/data"
srcfile="~/mvarga/worms/repo/fluid_substrate/main.cu"

#########################################
# check that program is build and ready #
#########################################

if [ ! -e $binfolder/nw3d ]; then
#	module load cuda
#	nvcc $srcfile -o $binfolder/nw3d -lcurand
	echo "Something"
fi

############################
# setup nw3d options array #
############################

ka_vals=$(awk 'BEGIN{for(i=0.25;i<=10;i+=0.25)print i}')
kbt_vals=$(awk 'BEGIN{for(i=0.0;i<1;i+=0.025)print i}')

id=0

for ka in $ka_vals
do
    for kbt in $kbt_vals
    do
	options[$id]="-o ${datafolder}/sim_${ka}_${kbt}.xyz -kbt ${kbt} -ka ${ka}"
	id=$id+1
    done
done

#####################################################
# using qsub to run nw3d instances with all options #
#####################################################

for opt in ${options[@]}
do
    echo ${opt}
done
