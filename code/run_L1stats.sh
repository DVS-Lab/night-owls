#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

td=1 # temporal derivatives -- 1 is on, 0 is off, everything else is invalid

for subinfo in "101 01" "101 02" "101 03" "101 06"; do

	# split subinfo variable
	set -- $subinfo
	sub=$1
	ses=$2

	for run in 1 2; do

		script=${scriptdir}/L1stats-sharedreward.sh
		NCORES=10 
		while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
			sleep 5s
		done
		bash $script $sub $ses $run $td &
		sleep 5s


		script=${scriptdir}/L1stats-mid.sh
		NCORES=10 
		while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
			sleep 5s
		done
		bash $script $sub $ses $run $td &
		sleep 5s

	done
done

