#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#for sub in `cat ${scriptdir}/newsubs.txt` ; do
for sub in 10317 10584 10590 10636 10638 10641 10642 10644 10647 10649 10652 10657 10659 10673 10674 10685 10690 10691 10700 10701 10716 10720 10723 10741 10774 10777 10781 10783 10800 10804; do

	script=${scriptdir}/fmriprep.sh
	NCORES=1 # need to do on OwlsNest with datalad since each sub has 8 runs of data with 4 echoes (needs 32 processors per sub)
	while [ $(ps -ef | grep -v grep | grep $script | wc -l) -ge $NCORES ]; do
		sleep 5s
	done
	bash $script $sub &
	sleep 5s
	
done

