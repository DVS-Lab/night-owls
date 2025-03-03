#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for subinfo in "101 01" "101 02"; do

	# split subinfo variable
	set -- $subinfo
	sub=$1
	ses=$2

	script=${scriptdir}/fmriprep.sh
	NCORES=3 # max should be 3 on Smith Lab Linux Box. Each session 4 runs of data with 4 echoes (will need 16 processors per session)
	while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
		sleep 5s
	done
	bash $script $sub $ses &
	sleep 5s
	
done

