
#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N L1stats-SR
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=28

# load modules and go to workdir
# module load fsl/6.0.2
# source $FSLDIR/etc/fslconf/fsl.sh
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
ppi=0 # dmn #dmn #VS  #dmn #0
sm=5 #trust & sr


for sub in ${subjects[@]}; do
	for ses in {01..12}; do
		for run in 1 2; do

			# set inputs and general outputs 
			MAINOUTPUT=${projectdir}/derivatives/fsl/level-run/space-MNI/sub-${sub}/ses-${ses}
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
			
			EVDIR=${projectdir}/derivatives/fsl/EVfiles/sub-${sub}/ses-${ses}/${TASK}/run-${run} # don't zeropad here since only 2 runs at most
			if [ ! -d "${projectdir}/derivatives/fsl/EVfiles/sub-${sub}/ses-${ses}/${TASK}/run-${run}" ]; then
				echo "missing EVfiles: $EVDIR " >> $logdir/re-runL1.log
				continue # skip these since some won't exist yet
			fi

			# Add variable naming for missed EV files based on task
			if [ $TASK == "trust" ]; then
				EVTITLE="missed_trial"
			else
				EVTITLE="missed_trial"
			fi

			# check for empty EVs (extendable to other studies)
			MISSED_TRIAL=${EVDIR}_${EVTITLE}.txt
			if [ -e $MISSED_TRIAL ]; then
				EV_SHAPE=3
			else
				EV_SHAPE=10
			fi

			# if network (ecn or dmn), do nppi; otherwise, do activation or seed-based ppi
			if [ "$ppi" == "ecn" -o "$ppi" == "dmn" ]; then

				# check for output and skip existing
				OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-1_type-nppi-${ppi}_run-${run}_sm-${sm}
				if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
					continue
				else
					echo "missing: $OUTPUT " >> $logdir/re-runL1.log
					rm -rf ${OUTPUT}.feat
				fi

				# network extraction. need to ensure you have run Level 1 activation
				MASK=${MAINOUTPUT}/L1_task-${TASK}_model-1_type-act_run-${run}_sm-${sm}.feat/mask
				if [ ! -e ${MASK}.nii.gz ]; then
					echo "cannot run nPPI because you're missing $MASK"
					continue
				fi
				for net in $(seq 0 9); do
					NET=${projectdir}/masks/networkmasks/seed-net${net}.nii.gz # this need to be changed to my project dir's naming of network masks
					TSFILE=${MAINOUTPUT}/ts_task-${TASK}_melodic-114_net${net}_nppi-${ppi}_run-${run}.txt
					fsl_glm -i $DATA -d $NET -o $TSFILE --demean -m $MASK
					eval INPUT${net}=$TSFILE
				done

				# set names for network ppi (we generally only care about ECN and DMN)
				DMN=$INPUT3
				ECN=$INPUT7
				if [ "$ppi" == "dmn" ]; then
					MAINNET=$DMN
					OTHERNET=$ECN
				else
					MAINNET=$ECN
					OTHERNET=$DMN
				fi

				# create template and run analyses
				ITEMPLATE=${projectdir}/templates/L1_task-${TASK}_model-1_type-nppi.fsf
				OTEMPLATE=${MAINOUTPUT}/L1_task-${TASK}_model-1_seed-${ppi}_run-${run}.fsf
				sed -e 's@OUTPUT@'$OUTPUT'@g' \
					-e 's@DATA@'$DATA'@g' \
					-e 's@EVDIR@'$EVDIR'@g' \
					-e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
					-e 's@EV_SHAPE@'$EV_SHAPE'@g' \
					-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
					-e 's@MAINNET@'$MAINNET'@g' \
					-e 's@OTHERNET@'$OTHERNET'@g' \
					-e 's@INPUT0@'$INPUT0'@g' \
					-e 's@INPUT1@'$INPUT1'@g' \
					-e 's@INPUT2@'$INPUT2'@g' \
					-e 's@INPUT4@'$INPUT4'@g' \
					-e 's@INPUT5@'$INPUT5'@g' \
					-e 's@INPUT6@'$INPUT6'@g' \
					-e 's@INPUT8@'$INPUT8'@g' \
					-e 's@INPUT9@'$INPUT9'@g' \
					-e 's@NVOLUMES@'$NVOLUMES'@g' \
					<$ITEMPLATE >$OTEMPLATE

			else # otherwise, do activation and seed-based ppi

				# set output based in whether it is activation or ppi
				if [ "$ppi" == "0" ]; then
					TYPE=act
					OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-1_type-${TYPE}_run-${run}_sm-${sm}
				else
					TYPE=ppi
					OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-1_type-${TYPE}_seed-${ppi}_run-${run}_sm-${sm}
				fi

				# check for output and skip existing
				if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
					continue
				else
					echo "missing: $OUTPUT " >> $logdir/re-runL1.log
					rm -rf ${OUTPUT}.feat
				fi

				# create template and run analyses
				ITEMPLATE=${projectdir}/templates/L1_task-${TASK}_model-1_type-${TYPE}.fsf
				if [ "$TASK" == "sharedreward" ]; then
					OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_ses-${ses}_task-${TASK}_model-1_type-${TYPE}_run-${run}.fsf
					if [ "$ppi" == "0" ]; then
						sed -e 's@OUTPUT@'$OUTPUT'@g' \
							-e 's@DATA@'$DATA'@g' \
							-e 's@EVDIR@'$EVDIR'@g' \
							-e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
							-e 's@EV_SHAPE@'$EV_SHAPE'@g' \
							-e 's@SMOOTH@'$sm'@g' \
							-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
							<$ITEMPLATE >$OTEMPLATE
					else
						PHYS=${MAINOUTPUT}/ts_task-${TASK}_mask-${ppi}_run-${run}.txt
						MASK=${projectdir}/masks/seed-${ppi}.nii.gz
						fslmeants -i $DATA -o $PHYS -m $MASK --eig
						sed -e 's@OUTPUT@'$OUTPUT'@g' \
							-e 's@DATA@'$DATA'@g' \
							-e 's@EVDIR@'$EVDIR'@g' \
							-e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
							-e 's@EV_SHAPE@'$EV_SHAPE'@g' \
							-e 's@PHYS@'$PHYS'@g' \
							-e 's@SMOOTH@'$sm'@g' \
							-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
							<$ITEMPLATE >$OTEMPLATE
					fi
				else
					OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_ses-${ses}_task-${TASK}_model-1_seed-${ppi}_run-${run}.fsf
					if [ "$ppi" == "0" ]; then
						sed -e 's@OUTPUT@'$OUTPUT'@g' \
							-e 's@DATA@'$DATA'@g' \
							-e 's@EVDIR@'$EVDIR'@g' \
							-e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
							-e 's@EV_SHAPE@'$EV_SHAPE'@g' \
							-e 's@SMOOTH@'$sm'@g' \
							-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
							<$ITEMPLATE >$OTEMPLATE
							#NVOLUMES hard coded in trust template
					else
						PHYS=${MAINOUTPUT}/ts_task-${TASK}_mask-${ppi}_run-${run}.txt
						MASK=${projectdir}/masks/seed-${ppi}.nii.gz
						fslmeants -i $DATA -o $PHYS -m $MASK
						sed -e 's@OUTPUT@'$OUTPUT'@g' \
							-e 's@DATA@'$DATA'@g' \
							-e 's@EVDIR@'$EVDIR'@g' \
							-e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
							-e 's@EV_SHAPE@'$EV_SHAPE'@g' \
							-e 's@PHYS@'$PHYS'@g' \
							-e 's@SMOOTH@'$sm'@g' \
							-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
							<$ITEMPLATE >$OTEMPLATE

					fi
				fi

			fi

			# add feat cmd to submission script
			echo feat $OTEMPLATE >>$logdir/cmd_feat_${PBS_JOBID}.txt

		done
	done
done

torque-launch -p "$logdir/chk_feat_${PBS_JOBID}.txt" "$logdir/cmd_feat_${PBS_JOBID}.txt"
