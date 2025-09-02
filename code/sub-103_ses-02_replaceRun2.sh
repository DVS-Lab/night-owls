#!/bin/bash

# Define the directory
DIR="/gpfs/scratch/tug87422/smithlab-shared/night-owls/bids/sub-103/ses-02/func"

# Loop over all run-3 files in that directory
for f in "$DIR"/*run-3*; do
    # Construct the run-2 filename
    target="${f/run-3/run-2}"
    
    # Only proceed if the run-2 file exists
    if [[ -e "$target" ]]; then
        echo "Replacing $target with $f"
        cp -f "$f" "$target"
    else
        echo "No run-2 match for $f, skipping."
    fi
done
