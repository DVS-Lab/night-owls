#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -N L1stats-SR-LSS
#PBS -q normal
#PBS -m ae
#PBS -M shenghan.wang@temple.edu
#PBS -l nodes=1:ppn=28

# load modules and go to workdir
# module load fsl/6.0.2
# source $FSLDIR/etc/fslconf/fsl.sh
#cd $PBS_O_WORKDIR

# ensure paths are correct
#shareddir=/gpfs/scratch/tug87422/smithlab-shared
#maindir=$shareddir/night-owls #this should be the only line that has to change if the rest of the script is set up correctly
#scriptdir=$maindir/code
#bidsdir=$maindir/bids
#logdir=$maindir/logsf
#mkdir -p $logdir

#rm -f $logdir/cmd_feat_${PBS_JOBID}.txt
#touch $logdir/cmd_feat_${PBS_JOBID}.txt

#rm -f L1stats-SR-LSS.o*
#rm -f L1stats-SR-LSS.e*

###########up here all hpc ##############

#1.need to think of  way to loop through sessions (which is different for each subject)
#2. loop the run

#!/usr/bin/env bash

# This script will perform Level 1 statistics in FSL.
# Rather than having multiple scripts, we are merging three analyses
# into this one script:
#		1) activation
#		2) seed-based ppi
#		3) network-based ppi
# Note that activation analysis must be performed first.
# Seed-based PPI and Network PPI should follow activation analyses.

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"



# study-specific inputs
sm=5 # this is already hard coded into all fsf files
subject=$1
ses=$2
TASK=sharedreward
run=$4
trial=`zeropad $5 2`
TYPE=act
#td=$5 # 1 for on, 0 for off (temporal derivatives)
MODEL=LSS # everyone should just have one model

#for sub in ${subjects[@]}; do 
#for ses in ${ses[@]}; do
#for trial in ${trial[@]}; do
#for run in 1 2; do
# set inputs and general outputs (should not need to chage across studies in Smith Lab)

MAINOUTPUT=${maindir}/derivatives/fsl/sub-${sub}
mkdir -p $MAINOUTPUT
DATA=${maindir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz

if [ ! -e $DATA ]; then
        echo " Exiting -- missing data: ${DATA}"
        exit
fi


NVOLUMES=`fslnvols ${DATA}`
CONFOUNDEVS=${maindir}/derivatives/fsl/confounds/sub-${sub}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-fslConfounds.tsv
#CONFOUNDEVS=${maindir}/derivatives/fsl/confounds_tedana/sub-${sub}/sub-${sub}_task-${TASK}_acq-${acq}_desc-TedanaPlusConfounds.tsv # switch to this later


if [ ! -e $CONFOUNDEVS ]; then
        echo "missing confounds: sub-${sub}_ses-${ses}_run-${run}"
        echo "missing confounds: sub-${sub}_ses-${ses}_run-${run}" >> ${maindir}/re-runL1_srLLS.log
        exit # exiting to ensure nothing gets run without confounds
fi


# EV files
EVDIR=${maindir}/derivatives/fsl/EVfiles/sub-${sub}/singletrial/ses-${ses}/${TASK}/run-${run}
if [ ! -e $EVDIR ]; then
        echo ${sub} ${acq} "EVDIR missing"
        echo "missing events files: $EVDIR " >> ${maindir}/re-runL1.log
        exit # exiting to ensure nothing gets run without confounds
fi


# empty EVs (specific to this study)
EV_MISSED_DEC=${EVDIR}/_miss_decision.txt
if [ -e $EV_MISSED_DEC ]; then
	SHAPE_MISSED_DEC=3
else
	SHAPE_MISSED_DEC=10
fi
EV_MISSED_OUTCOME=${EVDIR}/_miss_outcome.txt
if [ -e $EV_MISSED_OUTCOME ]; then
	SHAPE_MISSED_OUTCOME=3
else
	SHAPE_MISSED_OUTCOME=10
fi


SINGLETRIAL=${EVDIR}/ses-0${ses}run-${run}_SingleTrial${trial}.txt
OTHERTRIAL=${EVDIR}/run-0${run}_OtherTrials${trial}.txt


	OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_ses-${ses}_run-${run}_sm-${sm}_trial-${trial}


	# create common directory for zstat outputs
	zoutdir=${MAINOUTPUT}/LSS_task-${TASK}_model-_run-${run}_sm-${sm}
	if [ ! -d $zoutdir ]; then
	mkdir -p $zoutdir
	fi
	
	# check for output and skip existing
	if [ -e ${zoutdir}/zstat_trial-${trial}.nii.gz ]; then
		exit
	else
		echo "running: $OUTPUT " >> $logfile
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
		-e 's@SHAPE_MISSED_DEC@'$SHAPE_MISSED_DEC'@g' \
		-e 's@SHAPE_MISSED_OUTCOME@'$SHAPE_MISSED_OUTCOME'@g'\
		<$ITEMPLATE> $OTEMPLATE
    # add feat cmd to submission script
      echo feat $OTEMPLATE >>$logdir/cmd_feat_${PBS_JOBID}.txt

 #       done
#done

#torque-launch -p $logdir/chk_feat_${PBS_JOBID}.txt $logdir/cmd_feat_${PBS_JOBID}.txt

#end



############down here need work#####

# fix registration as per NeuroStars post:
# https://neurostars.org/t/performing-full-glm-analysis-with-fsl-on-the-bold-images-prepro>
#mkdir -p ${OUTPUT}.feat/reg
cp  $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/example_func2standard.mat
cp  $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/standard2example_func.mat
cp  ${OUTPUT}.feat/mean_func.nii.gz ${OUTPUT}.feat/reg/standard.nii.gz

# delete unused files
#rm -rf ${OUTPUT}.feat/stats/res4d.nii.gz
rm -rf ${OUTPUT}.feat/stats/corrections.nii.gz
rm -rf ${OUTPUT}.feat/stats/threshac1.nii.gz
rm -rf ${OUTPUT}.feat/filtered_func_data.nii.gz




# copy zstat image to common output folder and delete feat output
cp ${OUTPUT}.feat/stats/zstat1.nii.gz ${zoutdir}/zstat_trial-${trial}.nii.gz
rm -rf ${OUTPUT}.feat


