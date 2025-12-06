#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"


declare -A SKIP=(
  ["101:04"]=1
  ["101:05"]=1
  ["101:12"]=1
  ["103:12"]=1
)

#for sub in 101 103 104 105; do
for sub in 101; do
  for acq in multiecho single; do
    for confounds in tedana base; do
      for task in mid sharedreward; do
        for ses in {01..12}; do
          if [[ -n "${SKIP[$sub:$ses]}" ]]; then
            echo "Skipping sub $sub ses $ses"
            continue
          fi
          qsub -v task="$task",sub="$sub",confounds="$confounds",acq="$acq",ses="$ses" \
               "${scriptdir}/L1stats_LSS.qsub"
        done
      done
    done
  done
done
