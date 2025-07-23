#!/usr/bin/env bash


# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# study-specific inputs
TASK=sharedreward
sub=$1
ses=$2
run=$3
me=$4  # 1 on, 0 off

# zeropad the session number
ses=`zeropad ${ses} 2`

# use the optimally-combined data or the the second echo
if [ $me -eq 1 ]; then 
	INDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
	OUTDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_${sm}mm.nii.gz
	MASK=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-brain_mask.nii.gz
elif [ $me -eq 0 ]; then 
	# sub-101_ses-01_task-sharedreward_run-2_echo-2_part-mag_desc-preproc_bold.nii.gz
	# will need to move these to space-MNI152NLin6Asym (or whatever we want to use)
	INDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_echo-2_part-mag_desc-preproc_bold.nii.gz
	OUTDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_echo-2_part-mag_desc-preproc_bold_${sm}mm.nii.gz
	MASK=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_desc-brain_mask.nii.gz
else
	echo "exiting: invalid option for multiecho (me) variable. use 1 for on, 0 for off."
	exit
fi


if [ ! -e $INDATA ]; then
	echo "NO DATA: ${INDATA}"
	exit
fi

#only run if we're missing output
if [ -e $OUTDATA ]; then
	exit
else
	3dBlurToFWHM -FWHM $sm -input $INDATA -prefix $OUTDATA -mask $MASK
fi

# not yet sure how to suppress or control this output, but it conflicts with other processes (no overwrite)
rm -rf ${scriptdir}/3dFWHMx.1D ${scriptdir}/3dFWHMx.1D.png
