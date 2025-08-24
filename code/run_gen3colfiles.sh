#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

for subinfo in "105 01" "105 02" "105 03" "105 04" "105 05" "105 06" "105 07" "105 08" "105 09" "105 10" "105 11" "105 12"; do

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
