#!/usr/bin/env bash

# this script will convert your BIDS *events.tsv files into the 3-col format for FSL
# it relies on Tom Nichols' converter, which has been copied to our scriptdir to preserve modularity
# https://github.com/bids-standard/bidsutils



scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
baseout=${maindir}/derivatives/fsl/EVFiles
if [ ! -d ${baseout} ]; then
  mkdir -p $baseout
fi

sub=$1
ses=$2

for task in mid sharedreward; do
  for run in 1 2; do
    input=${maindir}/bids/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_events.tsv
    output=${baseout}/sub-${sub}/ses-${ses}/${task}/run-${run}

    if [ -e $input ]; then
    	  mkdir -p $output
      bash $scriptdir/BIDSto3col.sh $input ${output}/
    else
      echo "PATH ERROR: cannot locate ${input}."
      continue
    fi
  done
done