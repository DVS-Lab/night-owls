#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"



for sub in 101 103 104 105; do
	for acq in multiecho single; do 	
		for confounds in tedana base; do  
			for task in mid sharedreward; do
				for ses in {01..12}; do
					qsub -v task=${task},sub=${sub},confounds=${confounds},acq=${acq},ses=${ses} ${scriptdir}/L1stats_LSS.qsub
				done
			done	
		done
	done
done