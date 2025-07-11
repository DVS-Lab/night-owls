#!/bin/bash

# ensure paths are correct
maindir=/gpfs/scratch/tug47822/night-owls #this should be the only line that has to change if the rest of the script is set up correctly
scriptdir=$maindir/code


mapfile -t myArray < ${scriptdir}/sublist.txt


ntasks=9
counter=0
while [ $counter -lt ${#myArray[@]} ]; do
	subjects=${myArray[@]:$counter:$ntasks}
	echo $subjects
	let counter=$counter+$ntasks
	qsub -v subjects="${subjects[@]}" tedana-hpc-v2.sh
done
