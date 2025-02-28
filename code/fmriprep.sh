#!/bin/bash

# example code for FMRIPREP
# runs FMRIPREP on input subject and session
# usage: bash fmriprep.sh subject session
# example: bash fmriprep.sh 102 01

sub=$1
ses=$2 

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# make derivatives folder if it doesn't exist.
# let's keep this out of bids for now
if [ ! -d $maindir/derivatives ]; then
	mkdir -p $maindir/derivatives
fi

scratchdir=/ZPOOL/data/scratch/`whoami`
if [ ! -d $scratchdir ]; then
	mkdir -p $scratchdir
fi


TEMPLATEFLOW_DIR=/ZPOOL/data/tools/templateflow
export APPTAINERENV_TEMPLATEFLOW_HOME=/opt/templateflow

apptainer run --cleanenv \
-B ${TEMPLATEFLOW_DIR}:/opt/templateflow \
-B $maindir:/base \
-B /ZPOOL/data/tools/licenses:/opts \
-B $scratchdir:/scratch \
/ZPOOL/data/tools/fmriprep-24.1.1.simg \
/base/bids /base/derivatives/fmriprep \
participant --participant_label $sub \
--session-label $ses \
--stop-on-first-crash \
--me-output-echos \
--output-spaces MNI152NLin6Asym \
--bids-filter-file /base/code/fmriprep_config.json \
--fs-no-reconall --fs-license-file /opts/fs_license.txt -w /scratch


## Assistance from ChatGPT:
# Ensure your data follows BIDS format with ses-01, ses-02 directories.
# Run fMRIPrep without specifying a session to process all at once.
# Use --session-label ses-XX to process a single session at a time.
# Use --longitudinal to prevent anatomical image averaging across sessions.
# Expect separate preprocessing outputs for each session.

# --session-label $ses \
# --longitudinal \


