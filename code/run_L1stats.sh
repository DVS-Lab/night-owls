#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

task=sharedreward # edit if necessary

for denoise in "base" "tedana"; do # "base" = aCompCor confounds; "tedana" = aCompCor + tedana
	for ppi in 0 "VS_thr5"; do #"VS_thr5"; do #"VS_thr5"; do # putting 0 first will indicate "activation" "VS_thr5"
		for model in 1 2 3; do

			for sub in 10085 10094; do
			#for sub in `cat ${scriptdir}/sublist-sp.txt`; do # `ls -d ${basedir}/derivatives/fmriprep/sub-*/`
				sub=${sub#*sub-}
				sub=${sub%/}
			
			if [[ $sub == *sp ]]; then
				acqs=("mb2me4" "mb3me1fa50" "mb3me3" "mb3me3ip0" "mb3me4" "mb3me4fa50")
				#acqs="mb3me3ip0"
			else
				acqs=("mb1me1" "mb1me4" "mb3me1" "mb3me4" "mb6me1" "mb6me4")
			fi
			
			#for mbme in mb1me1 mb1me4 mb3me1 mb3me4 mb6me1 mb6me4 mb2me4 mb3me1fa50 mb3me3 mb3me3ip0 mb3me4fa50; do
			for mbme in "${acqs[@]}"; do

			  	# Manages the number of jobs and cores
			  	SCRIPTNAME=${basedir}/code/L1stats.sh
			  	NCORES=15
			  	while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
			    		sleep 5s
			  	done
			  	bash $SCRIPTNAME $model $sub $mbme $ppi $denoise &
				echo $SCRIPTNAME $model $sub $mbme $ppi $denoise &
					sleep 1s

			    	# done
			  done
			done
		done
	done
done
