#!/bin/bash

# ensure paths are correct
maindir=/gpfs/scratch/tug87422/smithlab-shared/night-owls #this should be the only line that has to change if the rest of the script is set up correctly
scriptdir=$maindir/code


mapfile -t myArray < ${scriptdir}/sublist.txt


# grab the first 10 elements
ntasks=1
counter=0
while [ $counter -lt ${#myArray[@]} ]; do
	subjects=${myArray[@]:$counter:$ntasks}
	echo $subjects
	let counter=$counter+$ntasks
	qsub -v subjects="${subjects[@]}" warpkit-hpc.sh
done
