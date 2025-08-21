#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

for subinfo in "101 01" "101 02" "101 03" "101 04" "101 05" "101 06" "101 07" "101 08" "101 09" "101 10" "101 11" "103 01" "103 02" "103 03" "103 04" "103 05" "103 06" "103 07" "103 08" "103 09" "103 10" "103 11" "103 12"; do

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
