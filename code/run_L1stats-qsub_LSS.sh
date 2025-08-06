#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

#mapfile -t myArray < "${scriptdir}/sublist.txt" 

# grab the first n elements
ntasks=1
counter=0

#still decide to manually add the sub-ses pair because of some sub has different number of sessions
#add total session number of each p in 
#for subinfo in "101 01" "101 02" "101 03" "101 04" "101 05" "101 06" "101 07" "101 08" "101 09" "101 10" "101 11"\
# "103 01" "103 02" "103 03" "103 04" "103 05" "103 06" "103 07" "103 08" "103 09" "103 10" "103 11" "103 12"\
 #"104 01" "104 02" "104 03" "104 04" "104 05" "104 06" "104 07" "104 08" "104 09" "104 10" "104 11" "104 12"; do 
 #"105 01" "105 02" "105 03" "105 04" "105 05" "105 06" "105 07" "105 08" "105 09" "105 10" "105 11" "105 12"
for subinfo in "101 01";do #test
# set is not working properly in seperating subinfo
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

				#myArray=($(seq 1 $trial))

				#for trial in $(seq 1 $t); do
				for trial in 1; do #test
				#while [ $counter -lt ${#myArray[@]} ]; do
				#trial=${myArray[@]:$counter:$ntasks}
				

					#dealing with bad trials, bad runs: dealing with miss_outcome and when loop till that one trial, continue
					##NOTEongoing task for: record the number of the missed trial that is spit out by the matlab EV gen script
					# and plug in the task-sub-ses-run-trial as a new elif-branch here,

						if [[ $task == "sharedreward"  &&  $sub -eq 101  &&  $ses -eq 3   &&  $run -eq 1   &&  $trial -eq 2 ]]; then
						continue
						elif [[ $task == "sharedreward"  &&  $sub -eq 101  &&  $ses -eq 4  &&  $run -eq 4   && ( $trial -eq 15 || $trial -eq 32 ) ]]; then
						continue
						elif [[ $task == "sharedreward"  &&  $sub -eq 103  &&  $ses -eq 2   &&  $run -eq 1   &&  $trial -eq 36 ]]; then
						continue				
						fi

					#let counter=$counter+$ntasks

							# Loop over each task script and submit with the same subject chunk
							qsub -v scriptname=${scriptname},sub=${sub},ses=${ses},run=${run},trial=${trial}  L1stats_LSS.qsub
				done

			done
		done
done

