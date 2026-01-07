#!/usr/bin/env bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"



# study-specific inputs
sm=0 # check templates to ensure no additional smoothing is being applied
sub=$1
ses=`zeropad $2 2`
TASK=sharedreward
run=$3
trial=`zeropad $4 2`
acq=$5
space=$6
confounds=$7
MODEL=LSS # everyone should just have one model
TYPE=act


# set inputs and general outputs (should not need to chage across studies in Smith Lab)
MAINOUTPUT=${maindir}/derivatives/fsl/sub-${sub}/LSS-FLOBS/ses-$ses/$TASK
mkdir -p $MAINOUTPUT
if [ "${acq}" == "multiecho" ]; then
	DATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-${space}_desc-preproc_bold.nii.gz
elif [ "${acq}" == "single" ]; then
	DATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold.nii.gz
else
	exit
fi

NVOLUMES=`fslnvols ${DATA}`

if [ "$confounds" == "tedana" ]; then
	CONFOUNDEVS=${maindir}/derivatives/fsl/confounds_tedana/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-TedanaPlusConfounds.tsv
elif [ "$confounds" == "base" ]; then
	CONFOUNDEVS=${maindir}/derivatives/fsl/confounds_tedana/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-fslConfounds.tsv
else
	exit
fi

if [ ! -e $CONFOUNDEVS ]; then
	echo "missing confounds: sub-${sub}_ses-${ses}_run-${run}"
	echo "missing: $CONFOUNDEVS" >> ${maindir}/re-runL1-LSS-${TASK}.log
	exit # exiting to ensure nothing gets run without confounds
fi

# EV files
EVDIR=${maindir}/derivatives/fsl/EVFiles/sub-${sub}/ses-${ses}/${TASK}/run-${run}/
SSLEVDIR=${maindir}/derivatives/fsl/EVFiles/sub-${sub}/singletrial/ses-${ses}/${TASK}/run-${run}/
SINGLETRIAL=${SSLEVDIR}run-${run}_SingleTrial${trial}.txt
OTHERTRIAL=${SSLEVDIR}run-${run}_OtherTrials${trial}.txt


# empty EVs (specific to sharedreward) don't work with basis functions in FEAT, so need separate template
EV_MISSED_DEC=${EVDIR}/_miss_decision.txt
if [ -e $EV_MISSED_DEC ]; then
	ITEMPLATE=${maindir}/templates/L1_task-${TASK}_model-${MODEL}_type-${TYPE}_FLOBS.fsf
else
	ITEMPLATE=${maindir}/templates/L1_task-${TASK}_model-${MODEL}_type-${TYPE}_FLOBS_NoMisses.fsf
fi


# set output based in whether it is activation or ppi
OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${MODEL}-type-${TYPE}_sub-${sub}_ses-${ses}_run-${run}_sm-${sm}_trial-${trial}_acq-${acq}_space-${space}_confounds-${confounds}_FLOBS


# check for output and skip existing
if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
    exit
else
    echo "running: $OUTPUT " >> "${maindir}/re-runL1-LSS-${TASK}.log"
    rm -rf "${OUTPUT}.feat"
fi

# create template and run analyses	
OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_task-${TASK}_model-${MODEL}_type-${TYPE}_ses-${ses}_run-${run}_sm-${sm}_trial-${trial}_acq-${acq}_space-${space}_confounds-${confounds}.fsf
sed -e 's@OUTPUT@'$OUTPUT'@g' \
-e 's@EVDIR@'$EVDIR'@g' \
-e 's@DATA@'$DATA'@g' \
-e 's@SINGLETRIAL@'$SINGLETRIAL'@g' \
-e 's@OTHERTRIAL@'$OTHERTRIAL'@g' \
-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
-e 's@NVOLUMES@'$NVOLUMES'@g' \
-e 's@FSLDIR@'$FSLDIR'@g' \
<$ITEMPLATE> $OTEMPLATE
feat $OTEMPLATE

# delete unused files
rm -rf ${OUTPUT}.feat/stats/res4d.nii.gz
rm -rf ${OUTPUT}.feat/stats/corrections.nii.gz
rm -rf ${OUTPUT}.feat/stats/threshac1.nii.gz
rm -rf ${OUTPUT}.feat/filtered_func_data.nii.gz
rm -rf ${OTEMPLATE}

exit
