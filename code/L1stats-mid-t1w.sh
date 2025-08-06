#!/usr/bin/env bash

# This script will perform Level 1 statistics in FSL.

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# study-specific inputs
TASK=mid
sm=0 # still needs to be done with afni (DVS needs to do this)
sub=$1
ses=$2
run=$3
td=0 # 1 for on, 0 for off (temporal derivatives)
model=1 # everyone should just have one model

# set inputs and general outputs (should not need to change across studies in Smith Lab)
MAINOUTPUT=${maindir}/derivatives/fsl/space-t1w/sub-${sub}/ses-${ses}
mkdir -p $MAINOUTPUT
DATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-T1w_desc-preproc_bold_5mm.nii.gz
if [ ! -e $DATA ]; then
	echo " Exiting -- missing data: ${DATA}"
	exit
fi

# check in template
NVOLUMES=`fslnvols $DATA`
#CONFOUNDEVS=${maindir}/derivatives/fsl/confounds/sub-${sub}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-fslConfounds.tsv
CONFOUNDEVS=${maindir}/derivatives/fsl/confounds_tedana/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-TedanaPlusConfounds.tsv  # switch to this later

if [ ! -e $CONFOUNDEVS ]; then
	echo "missing confounds: sub-${sub}_ses-${ses}_run-${run}"
	echo "missing confounds: sub-${sub}_ses-${ses}_run-${run}" >> ${maindir}/re-runL1.log
	exit # exiting to ensure nothing gets run without confounds
fi

EVDIR=${maindir}/derivatives/fsl/EVFiles/sub-${sub}/ses-${ses}/${TASK}/run-${run}/
if [ ! -e $EVDIR ]; then
	echo ${sub} ${acq} "EVDIR missing"
	echo "missing events files: $EVDIR " >> ${maindir}/re-runL1.log
	exit # exiting to ensure nothing gets run without confounds
fi



# set output based in whether it is activation and/or temporal derivatives
TYPE=act
if [ $td -eq 1 ]; then
	OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_run-${run}_td
	ITEMPLATE=${maindir}/templates/L1_task-${TASK}_model-${model}_type-${TYPE}_td.fsf
	OTEMPLATE=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_run-${run}_td.fsf
elif [ $td -eq 0 ]; then
	OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_run-${run}
	ITEMPLATE=${maindir}/templates/L1_task-${TASK}_model-${model}_type-${TYPE}.fsf
	OTEMPLATE=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_run-${run}.fsf
else
	echo "invalid parameter for temporal derivatives; it can only be 0 or 1."
	exit
fi

# check for output and skip existing
if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
	exit
else
	echo "missing feat output: $OUTPUT " >> ${maindir}/re-runL1.log
	rm -rf ${OUTPUT}.feat
fi
   
# create template and run analyses
sed -e 's@OUTPUT@'$OUTPUT'@g' \
-e 's@DATA@'$DATA'@g' \
-e 's@EVDIR@'$EVDIR'@g' \
-e 's@SMOOTH@'$sm'@g' \
-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
-e 's@NVOLUMES@'$NVOLUMES'@g' \
<$ITEMPLATE> $OTEMPLATE
feat $OTEMPLATE

# fix registration as per NeuroStars post:
# https://neurostars.org/t/performing-full-glm-analysis-with-fsl-on-the-bold-images-preprocessed-by-fmriprep-without-re-registering-the-data-to-the-mni-space/784/3
mkdir -p ${OUTPUT}.feat/reg
cp $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/example_func2standard.mat
cp $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/standard2example_func.mat
cp ${OUTPUT}.feat/mean_func.nii.gz ${OUTPUT}.feat/reg/standard.nii.gz

# delete unused files
rm -rf ${OUTPUT}.feat/stats/res4d.nii.gz
rm -rf ${OUTPUT}.feat/stats/corrections.nii.gz
rm -rf ${OUTPUT}.feat/stats/threshac1.nii.gz
rm -rf ${OUTPUT}.feat/filtered_func_data.nii.gz
