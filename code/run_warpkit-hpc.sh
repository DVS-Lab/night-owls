#!/bin/bash

# ensure paths are correct
maindir=/gpfs/scratch/tug87422/smithlab-shared/night-owls #this should be the only line that has to change if the rest of the script is set up correctly
scriptdir=$maindir/code

mapfile -t lines < "$scriptdir/sublist.txt"
pairs=()
for line in "${lines[@]}"; do
  # split into sub and ses
  sub=${line%%[[:space:]]*}
  ses=${line##*[[:space:]]}
  [[ -z $sub || -z $ses ]] && continue
  pairs+=( "${sub}:${ses}" )
done


# join into one space-delimited string
PAIR_LIST="${pairs[*]}"

echo "Submitting warpkit for: ${PAIR_LIST}"
qsub -v PAIRS="$PAIR_LIST" warpkit-hpc.sh