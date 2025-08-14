#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

for subinfo in "104 01" "104 02" "104 03" "104 04" "104 05" "104 06" "104 07" "104 08" "104 09" "104 10" "104 11" "104 12"; do

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
