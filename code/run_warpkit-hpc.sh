#!/bin/bash

# ensure paths are correct
maindir=/gpfs/scratch/tug87422/smithlab-shared/night-owls #this should be the only line that has to change if the rest of the script is set up correctly
scriptdir=$maindir/code

pairs=()
while read -r sub ses; do
  pairs+=( "${sub}:${ses}" )
done < "$scriptdir/sublist.txt"

# join into one space-delimited string
PAIR_LIST="${pairs[*]}"

echo "Submitting warpkit for: ${PAIR_LIST}"
qsub -v PAIRS="$PAIR_LIST" warpkit-hpc.sh