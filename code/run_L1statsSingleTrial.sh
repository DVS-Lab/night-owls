#!/bin/bash
# Ensure paths are correct irrespective from where user runs the script

#this should be the only line that has to change if the rest of the script is set up correctly
maindir=/gpfs/scratch/tug87422/smithlab-shared/night-owls
scriptdir=$maindir/code
basedir="$(dirname "$scriptdir")"
logs=$basedir/logs
logfile=${logs}/rerunL1-LSS_date-`date +"%FT%H%M"`.log
nruns=2 

mapfile -t myArray < "${scriptdir}/sublist-ses.txt" # ???? what should be in this one? this guys is in code directory and the script is run from this exact path as well

#1 task here need to be conditional on task, sub-ses-run-trial; 5 trial as one job - how? 

# grab the first 5 elements
ntasks=5
counter=0

#for task in 'sharedreward'; do
#for task in  'mid' ; do
      # for run in 1 2; do

                while [ $counter -lt ${#myArray[@]} ]; do
                        subjects=${myArray[@]:$counter:$ntasks}
                        echo $subjects
                        let counter=$counter+$ntasks
                        qsub -v subjects="${subjects[@]}" L1stats-hpc.sh
                done
				#	how to pass $task $sub $ses $run $trial to L1stats-hpc.sh with qsub??
                        sleep 1s
                #done     
        # done
#done
############## work in progress, with the specificity of trial number and missed trial in certain sub's session and certain run,
# don't know what to pass to qsub and what's the best way################

for subinfo in "101 01" "101 02" "101 03" "101 06"; do 

	# split subinfo variable
	set -- $subinfo
	sub=$1
	ses=$2

	for task in sharedreward mid; do
		for run in 1 2; do
			#directory for LSS single-trial events and regular event			
			EVDIR_LSS=${maindir}/derivatives/fsl/EVfiles/sub-${sub}/singletrial/ses-${ses}/${task}/run-${run}
			#EVDIR= ${maindir}/derivatives/fsl/EVfiles/sub-${sub}/singletrial/ses-${ses}/${task}/run-${run}			
 			
 			#get the max number of trial in that sub's that ses that task that run
			#t=$(ls ${EVDIR_LSS}/ses-0${ses}run-${run}_SingleTrial*.txt| sed -n 's/.*SingleTrial\([0-9]\+\)\.txt/\1/p' | sort ->

			if [ $task==sharedreward ]; then
			t=54
			else t=56
			fi
			 
			for trial in = $(seq 1 $t); do 

			#dealing with bad trials, bad runs
			#record the number of the trial that is spit out by the matlab script, plug in the if-statement here,
			# that is the miss_outcome and when loop till that one trial, continue
  											
			#need to manually add in the specificity of the missed trial
				if [$task == 'sharedreward'] [ $sub -eq 101] && [ $ses -eq 3]  && [ $run -eq 1]  && [ $trial -eq 2] ; then
				continue
				elif [$task == 'sharedreward'] [ $sub -eq 101] && [ $ses -eq 4]  && [ $run -eq 4]  && {[ $trial -eq 15] ||[ $trial -eq 32] }; then
				continue
				elif [$task == 'sharedreward'] [ $sub -eq 103] && [ $ses -eq 2]  && [ $run -eq 1]  && [ $trial -eq 36]; then
				fi

			script=${scriptdir}/L1statsSingleTrial-${task}.sh
			NCORES=20 
				while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
				sleep 5s
				done
			bash $script $sub $ses $task $run $trial  #$td?
			sleep 5s

			done
		done
	done
done
