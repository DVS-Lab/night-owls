#!/usr/bin/env bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"



# study-specific inputs
sm=0 # this is already hard coded into all fsf files
sub=$1
ses=`zeropad $2 2`
TASK=mid
run=$3
trial=`zeropad $4 2`
acq=$5
space=$6
confounds=$7
MODEL=LSS # everyone should just have one model
TYPE=act

			
# set inputs and general outputs (should not need to chage across studies in Smith Lab)
MAINOUTPUT=${maindir}/derivatives/fsl/sub-${sub}
mkdir -p $MAINOUTPUT

if [ "${acq}" == single ]; then
DATA=${projectdir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-${space}_desc-preproc_bold.nii.gz
else
DATA=${projectdir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold.nii.gz
fi

NVOLUMES=`fslnvols ${DATA}`

#####confoundev need 1. ask if based confound generated from Matt 2. if statement based on value of ${confounds}
CONFOUNDEVS=${maindir}/derivatives/fsl/confounds_tedana/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-TedanaPlusConfounds.tsv


if [ ! -e $CONFOUNDEVS ]; then
	echo "missing confounds: sub-${sub}_ses-${ses}_run-${run}"
	echo "missing: $CONFOUNDEVS " >> ${maindir}/re-runL1_midLSS.log
	exit # exiting to ensure nothing gets run without confounds
fi

# EV files
EVDIR=${maindir}/derivatives/fsl/EVFiles/sub-${sub}/ses-${ses}/${TASK}/run-${run}/
SSLEVDIR=${maindir}/derivatives/fsl/EVFiles/sub-${sub}/singletrial/ses-${ses}/${TASK}/run-${run}/
SINGLETRIAL=${SSLEVDIR}run-${run}_SingleTrial${trial}.txt
OTHERTRIAL=${SSLEVDIR}run-${run}_OtherTrials${trial}.txt

# create common directory for zstat outputs
zoutdir=${MAINOUTPUT}/LSS_task-${TASK}_sub-${sub}_ses-${ses}_run-${run}_sm-${sm}
if [ ! -d $zoutdir ]; then
	mkdir -p $zoutdir
fi

# set output based in whether it is activation or ppi
OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${MODEL}-type-${TYPE}_sub-${sub}_ses-${ses}_run-${run}_sm-${sm}_trial-${trial}_acq-${acq}_space-${space}_confounds-${confounds}


# check for output and skip existing
if [ -e ${zoutdir}/zstat_trial-${trial}.nii.gz ]; then
	exit
else
	echo "running: $OUTPUT " >> ${maindir}/re-runL1-LSS-${TASK}.log
	rm -rf ${OUTPUT}.feat
fi

# create template and run analyses
ITEMPLATE=${maindir}/templates/L1_task-${TASK}_model-${MODEL}_type-${TYPE}.fsf
OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_task-${TASK}_model-${MODEL}_type-${TYPE}_ses-${ses}_run-${run}_sm-${sm}_trial-${trial}.fsf
sed -e 's@OUTPUT@'$OUTPUT'@g' \
-e 's@EVDIR@'$EVDIR'@g' \
-e 's@DATA@'$DATA'@g' \
-e 's@SINGLETRIAL@'$SINGLETRIAL'@g' \
-e 's@OTHERTRIAL@'$OTHERTRIAL'@g' \
-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
-e 's@NVOLUMES@'$NVOLUMES'@g' \
<$ITEMPLATE> $OTEMPLATE
feat $OTEMPLATE 


# fix registration as per NeuroStars post:
# https://neurostars.org/t/performing-full-glm-analysis-with-fsl-on-the-bold-images-prepro>
mkdir -p ${OUTPUT}.feat/reg
cp $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/example_func2standard.mat
cp $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/standard2example_func.mat
cp ${OUTPUT}.feat/mean_func.nii.gz ${OUTPUT}.feat/reg/standard.nii.gz

# delete unused files
rm -rf ${OUTPUT}.feat/stats/res4d.nii.gz
rm -rf ${OUTPUT}.feat/stats/corrections.nii.gz
rm -rf ${OUTPUT}.feat/stats/threshac1.nii.gz
rm -rf ${OUTPUT}.feat/filtered_func_data.nii.gz



# copy zstat image to common output folder and delete feat output
cp ${OUTPUT}.feat/stats/zstat1.nii.gz ${zoutdir}/zstat_trial-${trial}.nii.gz
chmod -R 777 ${zoutdir}
rm -rf ${OUTPUT}.feat
chmod -R 777 ${MAINOUTPUT}
