#!/bin/bash

# check if the commands.txt file is provided as a command line argument
if [ $# -ne 1 ]; then
	echo "Usage $0 <commands.txt>"
	exit 1
fi

# read the commands file from the command-line argument
COMMANDS_FILE="$1"

# check if the specified file exists
if [ ! -f "$COMMANDS_FILE" ]; then
	echo "Error: the specified commands file does not exist."
	exit 1
fi

# Define the program you want to run
PROGRAM="./amatter.x"

# create a directory to store the results
RESULTS_DIR="data_04_10_2023"
mkdir -p "$RESULTS_DIR"

# determine the number of commands in the file and set --jobs accordingly
NUM_COMMANDS=$(wc -l < "$COMMANDS_FILE")
echo "Number of jobs: $NUM_COMMANDS"
cat "$COMMANDS_FILE"
echo "----------STARTING------------"

# Use cat and pipe to read commands from the file and run them in parallel
# cat "$COMMANDS_FILE" | parallel --jobs "$NUM_COMMANDS" --line-buffer --halt now,fail --results "$RESULTS_DIR/results_{#}" ::: \
# taskset -c "$(( {#} - 1 ))" {}
(cat "$COMMANDS_FILE" | while read -r command; do
	(i=0
	export j=$((i+1))
	output_dir="$RESULTS_DIR/results_$j"
	mkdir -p "$output_dir" # create a directory for the job's output
	(cd "$output_dir" && taskset -c "$i" "$command") & # run the job in its output directory
	pid=$! # get the pid of the background job
	# cpu_core=$i # calculate the cpu core based on the job number
	# taskset -c "$cpu_core" -p "$pid" #pin the job to the cpu core
	echo "$command"
	echo "$pid is on $i"
	)
done; wait)

echo "done!"
