#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# grab the first n elements
ntasks=28
counter=0

#still decide to manually add the sub-ses pair because of some sub has different number of sessions
#add total session number of each p in 
for subinfo in "101 01" "101 02" "101 03" "101 04" "101 05" "101 06" "101 07" "101 08" "101 09" "101 10" "101 11"; do  #\ after other subject's fmriprep done, uncomment the other
	# "103 01" "103 02" "103 03" "103 04" "103 05" "103 06" "103 07" "103 08" "103 09" "103 10" "103 11" "103 12"\
 #"104 01" "104 02" "104 03" "104 04" "104 05" "104 06" "104 07" "104 08" "104 09" "104 10" "104 11" "104 12"; do 
 #"105 01" "105 02" "105 03" "105 04" "105 05" "105 06" "105 07" "105 08" "105 09" "105 10" "105 11" "105 12"

#for subinfo in "101 03";do #test
	# split subinfo variable
	set -- $subinfo
	sub=$1
	ses=$2

	for task in mid sharedreward; do
			scriptname=${scriptdir}/L1statsSingleTrial-${task}.sh

		for run in 1 2; do
				if [ $task == "sharedreward" ]; then
				t=54
				else 
				t=56
				fi

				myArray=($(seq 1 $trial))				
				#for trial in 2; do #test
			while [ $counter -lt ${#myArray[@]} ]; do
				trial=${myArray[@]:$counter:$ntasks}
				let counter=$counter+$ntasks
				# Loop over each task script and submit with the same subject chunk
				qsub -v scriptname=${scriptname},task=${task},sub=${sub},ses=${ses},run=${run},trial=${trial@}  L1stats_LSS.qsub
			#unsure whether these three should be here or .qsub: acq=${acq}, space=${space}, confounds=${confounds}
			done	
		done
	done
done