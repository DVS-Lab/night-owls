#!/bin/bash

DIR="/ZPOOL/data/projects/night-owls/bids/sub-103/ses-02/func"

cd "$DIR" || { echo "Directory not found"; exit 1; }

shopt -s nullglob

for f in *run-3*; do
    target="${f/run-3/run-2}"

    if [[ -f "$target" ]]; then
        if mv -f "$f" "$target"; then
            echo "Moved $f → $target"
        else
            echo "ERROR: Failed to move $f → $target" >&2
        fi
    else
        echo "Skipping $f: no matching run-2" >&2
    fi
done
