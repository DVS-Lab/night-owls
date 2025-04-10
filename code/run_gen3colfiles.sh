#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

for subinfo in "101 01" "101 02" "101 03" "101 04" "101 05" "101 06" "101 07" "101 08" "101 09" "103 01" "103 02"; do

	# split subinfo variable
	set -- $subinfo
	sub=$1
	ses=$2

	for run in 1 2; do

		script=${scriptdir}/gen3colfiles.sh
		NCORES=10
		while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
			sleep 5s
		done
		bash $script $sub $ses &
		sleep 5s
	done
	
done
