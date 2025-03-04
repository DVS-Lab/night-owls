#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for subinfo in 101; do

	# split subinfo variable
	set -- $subinfo
	sub=$1

	script=${scriptdir}/fmriprep.sh
	NCORES=1 # Tricky. Each session 4 runs of data with 4 echoes (will need 16 processors per session). Will do all sessions at once. :-()
	while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
		sleep 5s
	done
	bash $script $sub
	sleep 5s
	
done

