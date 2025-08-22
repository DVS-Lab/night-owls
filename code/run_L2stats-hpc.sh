#!/bin/bash

# Ensure paths are correct irrespective from where user runs the script

projectdir=/gpfs/scratch/tug87422/smithlab-shared/night-owls
scriptdir=$projectdir/code
basedir="$(dirname "$scriptdir")"

mapfile -t myArray < "${scriptdir}/sublist.txt" 

# grab the first n elements
ntasks=1
counter=0
		
while [ $counter -lt ${#myArray[@]} ]; do
	subjects=${myArray[@]:$counter:$ntasks}
	let counter=$counter+$ntasks

    script="L2stats-hpc.sh"
    qsub -v subjects="${subjects[@]}" "$script"
    echo $subjects $script
    
    #script_subj="L2stats-hpc-subj.sh"
    #qsub -v subjects="${subjects[@]}" "$script_subj"
    #echo $subjects $script_subj
done
