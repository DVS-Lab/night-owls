#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


# these define the subject number (sub) and session id (ses)
for subinfo in "101 06"; do

	# split subinfo variable
	set -- $subinfo
	sub=$1
	ses=$2


	script=${scriptdir}/prepdata.sh
	NCORES=20
	while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
		sleep 5s
	done
   	bash $script $sub $ses &
	sleep 5s
	
done
