#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# grab the first n elements
ntasks=1
counter=0

#ntasks_sub=1
#counter_sub=0

#still decide to manually add the sub-ses pair because of some sub has different number of sessions
#add total session number of each p in 
for sub in 101 103 104 105; do
	for acq in 'multiecho' 'single';do 	
		for confounds in 'tedana' 'base'; do  #need to figure out is the based confound generated at all
			for task in mid sharedreward; do
			scriptname=${scriptdir}/L1statsSingleTrial-${task}.sh

			if [ $task == "sharedreward" ]; then
			t=54
			else 
			t=56
			fi

			myArray=($(seq 1 $t))				
			#for trial in 2; do #test
				
			while [ $counter -lt ${#myArray[@]} ]; do
			trial=${myArray[@]:$counter:$ntasks}
			let counter=$counter+$ntasks
			# Loop over each task script and submit with the same subject chunk
			qsub -v scriptname=${scriptname},task=${task},sub=${sub},confounds=${confounds},acq=${acq},trial=${trial@}  L1stats_LSS.qsub
			done
	
			done	
		done
	done
done