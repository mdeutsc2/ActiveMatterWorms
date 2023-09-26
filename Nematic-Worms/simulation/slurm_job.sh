#!/bin/bash
#SBATCH --job-name=active-nematic-simulation
#SBATCH --account=PGS0213
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mvarga6@kent.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus=1
#SBATCH --time=5:00:00

echo "Starting active filaments job as user '$USER'"
pwd; hostname; date

# The location of the simulation binary on head node
LOCAL_BIN=/users/PGS0213/ksu0236/mvarga/Nematic-Worms/simulation

# Get name of the temporary directory working directory, physically on the compute node
WORK_DIR="${TMPDIR}"

# Get submit directory
# (every file/folder below this directory is copied to the compute node)
SUBMIT_DIR="${SLURM_SUBMIT_DIR}"

echo "LOCAL_BIN:  $LOCAL_BIN"
echo "WORK_DIR:   $WORK_DIR"
echo "SUBMIT_DIR: $SUBMIT_DIR"
echo "NW_ARGS:    $NW_ARGS"

echo
echo "======================"
echo "  SLURM ENVIRONMENT"
echo "====================="
printenv | grep SLURM

# Copy simluation binary to compute-node
cp "$LOCAL_BIN/nw" "$WORK_DIR"
cd "$WORK_DIR"

# define clean-up function
function clean_up {
  # - delete temporary files from the compute-node, before copying
  rm -r "$WORK_DIR/nw"
  # - compress XYZ files
  while IFS= read -r file_name; do
    base_name=$(basename "${file_name}" .xyz)
    tar -cJvf "${base_name}.tar.xz" "${file_name}"
  done < <( ls {*.xyz,*.xyzv} 2> /dev/null )
  # - change directory to the location of the sbatch command (on the head node)
  cd "${SUBMIT_DIR}"
  # - copy compressed files from the temporary directory on the compute-node
  cp -prf "${WORK_DIR}"/*.tar.xz .
  # - erase the temporary directory from the compute-node
  rm -rf "${workdir}"
  # - exit the script
  exit
}

trap 'clean_up' EXIT

module load cuda
time ./nw $NW_ARGS