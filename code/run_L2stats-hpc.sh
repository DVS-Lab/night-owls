#!/bin/bash

# Ensure paths are correct irrespective from where user runs the script

projectdir=/gpfs/scratch/tug87422/smithlab-shared/night-owls
scriptdir=$projectdir/code
basedir="$(dirname "$scriptdir")"
nruns=2 

mapfile -t myArray < "${scriptdir}/sublist.txt" 

# grab the first n elements
ntasks=1
counter=0

#tasks=("SR" "MID")
tasks=("SR")

#for task in 'sharedreward'; do
#for task in  'trust' ; do
		
while [ $counter -lt ${#myArray[@]} ]; do
	subjects=${myArray[@]:$counter:$ntasks}
	let counter=$counter+$ntasks

		# Loop over each task script and submit with the same subject chunk
	for task in "${tasks[@]}"; do
		script="L2stats-hpc-${task}.sh"
		qsub -v subjects="${subjects[@]}" "$script"
        echo $subjects $script

        script-subj="L2stats-hpc-${task}-subj.sh"
		qsub -v subjects="${subjects[@]}" "$script"
        echo $subjects $script-subj
	done
done
