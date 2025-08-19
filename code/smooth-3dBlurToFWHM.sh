#!/usr/bin/env bash


# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# study-specific inputs
sub=$1
ses=$2
TASK=$3
run=$4
me=$5  # 1 on, 0 off
space=$6
sm=5  # 5mm smoothing

# zeropad the session number
ses=`zeropad ${ses} 2`

# use the optimally-combined data or the the second echo
if [ $me -eq 1 ]; then 
	INDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-${space}_desc-preproc_bold.nii.gz
	OUTDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-${space}_desc-preproc_bold_${sm}mm.nii.gz
elif [ $me -eq 0 ]; then 
	INDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold.nii.gz
	OUTDATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold_${sm}mm.nii.gz
else
	echo "exiting: invalid option for multiecho (me) variable. use 1 for on, 0 for off."
	exit
fi

# only two spaces, so two masks per run
MASK=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-${space}_desc-brain_mask.nii.gz


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
