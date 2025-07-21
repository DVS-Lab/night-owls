
#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N L1stats-SR
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=28

# load modules and go to workdir
module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
cd $PBS_O_WORKDIR

# ensure paths are correct
shareddir=/gpfs/scratch/tug87422/smithlab-shared
projectdir=$shareddir/night-owls
scriptdir=$projectdir/code
bidsdir=$projectdir/bids
logdir=$projectdir/logs
mkdir -p $logdir

rm -f $logdir/cmd_feat_${PBS_JOBID}.txt
touch $logdir/cmd_feat_${PBS_JOBID}.txt

rm -f L1stats-SR.o*
rm -f L1stats-SR.e*


TASK=sharedreward
sm=5 #mid & sr
TYPE=act
ppi=0
sm=5 #trust & sr

rm $logdir/re-runL1.log

for sub in ${subjects[@]}; do
	for ses in {01..12}; do
		for run in 1 2; do

			# set inputs and general outputs 
			MAINOUTPUT=${projectdir}/derivatives/fsl/level-run/space-MNI/sub-${sub}/ses-${ses}/
			mkdir -p $MAINOUTPUT
			DATA=${projectdir}/derivatives/fmriprep/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
		
			#Check for input data, skip if missing
			if [ ! -e $DATA ]; then
				echo "MISSING FMRIPREP INPUT ${sub}, run ${run}: $DATA" >> $logdir/re-runL1.log
				continue
			fi
			

			NVOLUMES=$(fslnvols $DATA)
			CONFOUNDEVS=${projectdir}/derivatives/fsl/confounds_tedana/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-TedanaPlusConfounds.tsv
			if [ ! -e $CONFOUNDEVS ]; then
				echo "missing: $CONFOUNDEVS " >> $logdir/re-runL1.log
				continue # exiting/continuing to ensure nothing gets run without confounds
			fi
			
			EVDIR=${projectdir}/derivatives/fsl/EVFiles/sub-${sub}/ses-${ses}/${TASK}/run-${run}/ # don't zeropad here since only 2 runs at most
			if [ ! -d "${projectdir}/derivatives/fsl/EVFiles/sub-${sub}/ses-${ses}/${TASK}/run-${run}/" ]; then
				echo "missing EVFiles: $EVDIR " >> $logdir/re-runL1.log
				continue # skip these since some won't exist yet
			fi

			
			# check for BOTH  empty missed EVs (specific to this study)
	
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
			
			OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-1_type-${TYPE}_run-${run}_sm-${sm}

			# check for output and skip existing
			if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
				continue
			else
			echo "missing: $OUTPUT " >> $logdir/re-runL1.log
				rm -rf ${OUTPUT}.feat
			fi

			# create template and run analyses
			ITEMPLATE=${projectdir}/templates/L1_task-${TASK}_model-1_type-${TYPE}.fsf
			OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_ses-${ses}_task-${TASK}_model-1_type-${TYPE}_run-${run}.fsf
				sed -e 's@OUTPUT@'$OUTPUT'@g' \
					-e 's@DATA@'$DATA'@g' \
					-e 's@EVDIR@'$EVDIR'@g' \
					-e 's@SHAPE_MISSED_DEC@'$SHAPE_MISSED_DEC'@g' \
					-e 's@SHAPE_MISSED_OUTCOME@'$SHAPE_MISSED_OUTCOME'@g' \
					-e 's@EV_SHAPE@'$EV_SHAPE'@g' \
					-e 's@SMOOTH@'$sm'@g' \
					-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
					-e 's@NVOLUMES@'$NVOLUMES'@g' \
					<$ITEMPLATE >$OTEMPLATE					


			# add feat cmd to submission script
			echo feat $OTEMPLATE >>$logdir/cmd_feat_${PBS_JOBID}.txt

		done
	done
done

torque-launch -p "$logdir/chk_feat_${PBS_JOBID}.txt" "$logdir/cmd_feat_${PBS_JOBID}.txt"
