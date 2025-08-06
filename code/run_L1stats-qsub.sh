#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

mapfile -t myArray < "${scriptdir}/sublist.txt" 

# grab the first n elements
ntasks=1
counter=0
		
while [ $counter -lt ${#myArray[@]} ]; do
	subjects=${myArray[@]:$counter:$ntasks}
	let counter=$counter+$ntasks

	# Loop over each task script and submit with the same subject chunk
	qsub -v subjects="${subjects[@]}" L1stats.qsub
	qsub -v subjects="${subjects[@]}" L1stats-t1w.qsub
done
