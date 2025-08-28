
#!/usr/bin/env bash

umask 0000

#for sub in 101 103 104 105; do
for sub in 104; do
  for ses in $(seq -w 01 12); do
    qsub fmriprep-anat/fmriprep_sub-${sub}_ses-${ses}.sh
  done
done
