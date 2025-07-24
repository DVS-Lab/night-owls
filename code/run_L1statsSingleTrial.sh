#!/bin/bash
# Ensure paths are correct irrespective from where user runs the script

#this should be the only line that has to change if the rest of the script is set up correctly

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

#1 task here need to be conditional on task, sub-ses-run-trial; 5 trial as one job - how? 

# grab the first 5 elements
ntasks=12
counter=0

for subinfo in "101 03" "101 06" "103 02"; do 

# set is not working properly in seperating subinfo
	# split subinfo variable
	set -- $subinfo
	sub=$1
	ses=$2

	for task in "sharedreward"; do # mid
		for run in 1 2; do
			#directory for LSS single-trial events and regular event			
			EVDIR_LSS=${maindir}/derivatives/fsl/EVfiles/sub-${sub}/singletrial/ses-${ses}/${task}/run-${run}
			#EVDIR= ${maindir}/derivatives/fsl/EVfiles/sub-${sub}/singletrial/ses-${ses}/${task}/run-${run}			
 			
 			#get the max number of trial in that sub's that ses that task that run
			#t=$(ls ${EVDIR_LSS}/ses-0${ses}run-${run}_SingleTrial*.txt| sed -n 's/.*SingleTrial\([0-9]\+\)\.txt/\1/p' | sort ->

			if [ $task == "sharedreward" ]; then
			t=54
			else 
			t=56
			fi
			 
			for trial in $(seq 1 $t); do 

			#dealing with bad trials, bad runs
			#record the number of the trial that is spit out by the matlab script, plug in the if-statement here,
			# that is the miss_outcome and when loop till that one trial, continue
  											
			#need to manually add in the specificity of the missed trial
				if [[ $task == "sharedreward"  &&  $sub -eq 101  &&  $ses -eq 3   &&  $run -eq 1   &&  $trial -eq 2 ]]; then
				continue
				elif [[ $task == "sharedreward"  &&  $sub -eq 101  &&  $ses -eq 4  &&  $run -eq 4   && ( $trial -eq 15 || $trial -eq 32 ) ]]; then
				continue
				elif [[ $task == "sharedreward"  &&  $sub -eq 103  &&  $ses -eq 2   &&  $run -eq 1   &&  $trial -eq 36 ]]; then
				continue				
				fi

			script=${scriptdir}/L1statsSingleTrial-${task}.sh
			NCORES=20 
				while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
				sleep 5s
				done
			bash $script $sub $ses $run $trial  #$td?
			sleep 5s

			done
		done
	done
done
