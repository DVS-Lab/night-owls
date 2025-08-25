#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N L2stats-subj
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=14

# load modules and go to workdir
# module load fsl/6.0.2
# source $FSLDIR/etc/fslconf/fsl.sh
cd $PBS_O_WORKDIR
umask 0000

# ensure paths are correct
shareddir=/gpfs/scratch/tug87422/smithlab-shared
projectdir=$shareddir/night-owls
scriptdir=$projectdir/code
bidsdir=$projectdir/bids
logdir=$projectdir/logs
mkdir -p $logdir

rm -f $logdir/cmd_L2_${PBS_JOBID}.txt
touch $logdir/cmd_L2_${PBS_JOBID}.txt

rm -f L2stats-sub.o*
rm -f L2stats-sub.e*



type="act"               # "act" or "ppi" (or "nppi-dmn")
#sm=5                    # smoothing kernel label
model=1                 # first-level model number
tasks=("sharedreward" "mid")
spaces=(mni t1w)
echos=(single-echo multi-echo)
confounds=(cnfds-fmriprep cnfds-tedana)

for task in "${tasks[@]}"; do

    for sub in ${subjects[@]}; do
        for space in "${spaces[@]}"; do
            for echo in "${echos[@]}"; do
                for confound in "${confounds[@]}"; do

        rm -f $logdir/re-runL2subj_sub-${sub}.log
        MAINOUTPUT=${projectdir}/derivatives/fsl/sub-${sub}
        
        # Initialize arrays to store all available inputs for this subject
        all_inputs=()
        input_labels=()
        
        # Collect all available L1 outputs across all sessions and runs
        for ses in {01..12}; do
            SESDIR=${MAINOUTPUT}/ses-${ses}
            
            # skip if the session folder itself doesn't exist
            if [ ! -d "${SESDIR}" ]; then
                echo "SKIP sub-${sub} ses-${ses} ${task}: session directory not found" >> $logdir/re-runL2subj_sub-${sub}.log
                continue
            fi
            
            # Check each run for this session
            for run in 1 2; do
                INPUT=${SESDIR}/L1_task-${task}_model-${model}_type-${type}_run-${run}_space-${space}_${echo}_${confound}.feat
                
                if [ -d "${INPUT}" ]; then
                    all_inputs+=(${INPUT})
                    input_labels+=("ses-${ses}_run-${run}")
                else
                    echo "MISSING sub-${sub} ses-${ses} run-${run} ${task} ${space} ${echo} ${confound}: ${INPUT}" >> $logdir/re-runL2subj_sub-${sub}.log
                fi
            done
        done
        
        # Skip subject if fewer than 2 valid inputs found (need at least 2 for L2 analysis)
        if [ ${#all_inputs[@]} -lt 2 ]; then
            echo "SKIP sub-${sub} ${task} ${space} ${echo} ${confound}: insufficient inputs (${#all_inputs[@]} found), need at least 2" >> $logdir/re-runL2subj_sub-${sub}.log
            continue
        fi
        
        # Set NSES to number of inputs
        NSES=${#all_inputs[@]}

        # Set output path for subject-level analysis
        mkdir ${MAINOUTPUT}/subject-level/
        OUTPUT=${MAINOUTPUT}/subject-level/L2_subj-${sub}_task-${task}_model-${model}_type-${type}_space-${space}_${echo}_${confound}
        
        # skip if output already exists
        if [ -e ${OUTPUT}.gfeat/cope${NCOPES}.feat/cluster_mask_zstat1.nii.gz ]; then
            echo "SKIP sub-${sub} ${task}: L2 already done" >> $logdir/re-runL2subj.log
            continue
        fi
        
        # set template based on type
        #if [ "${type}" == "act" ]; then
        #    ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_type-act_subject-level.fsf
        #    NCOPES=34
        #else
        #    ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_type-nppi-dmn_subject-level.fsf
        #    NCOPES=$((34 + 1))
        #fi
        
        ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_type-act_subject-level.fsf
        # Create dynamic FSF template
        OTEMPLATE=${MAINOUTPUT}/subject-level/L2_task-${task}_model-${model}_type-${type}_subj-${sub}_space-${space}_${echo}_${confound}.fsf
        
        # Start with base template and modify for dynamic inputs
        cp ${ITEMPLATE} ${OTEMPLATE}
        
        # Update basic parameters
        sed -i "s@OUTPUT@${OUTPUT}@g" ${OTEMPLATE}
        sed -i "s@NSES@${NSES}@g" ${OTEMPLATE}
            
        # Add inputs
        for i in "${!all_inputs[@]}"; do
            feat_index=$((i + 1))
            echo "set feat_files(${feat_index}) \"${all_inputs[i]}\"" >> ${OTEMPLATE}
            echo "set fmri(evg${feat_index}.1) 1.0" >> ${OTEMPLATE}
            echo "set fmri(groupmem.${feat_index}) 1" >> ${OTEMPLATE}
        done
        
        # Log what we're processing
        echo "PROCESSING sub-${sub}: ${#all_inputs[@]} inputs (${input_labels[*]})" >> $logdir/re-runL2subj_sub-${sub}.log
        
        # Add to command file
        echo feat $OTEMPLATE >>$logdir/cmd_L2_${PBS_JOBID}.txt

    done
done



torque-launch -p "$logdir/chk_L2_${PBS_JOBID}.txt" "$logdir/cmd_L2_${PBS_JOBID}.txt"


# delete unused files
for task in "${tasks[@]}"; do
    for sub in ${subjects[@]}; do
            for space in "${spaces[@]}"; do
            for echo in "${echos[@]}"; do
            for confound in "${confounds[@]}"; do
        OUTPUT=${projectdir}/derivatives/fsl/sub-${sub}/subject-level/L2_subj-${sub}_task-${task}_model-${model}_type-${type}_space-${space}_${echo}_${confound}
        
        # Loop through cope numbers based on task
        if [[ "$task" == "sharedreward" ]]; then
            for cope_num in {1..15}; do
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/stats/res4d.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/stats/corrections.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/stats/threshac1.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/filtered_func_data.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/var_filtered_func_data.nii.gz
            done
        elif [[ "$task" == "mid" ]]; then
            for cope_num in {1..10}; do
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/stats/res4d.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/stats/corrections.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/stats/threshac1.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/filtered_func_data.nii.gz
                rm -f ${OUTPUT}.gfeat/cope${cope_num}.feat/var_filtered_func_data.nii.gz
            done
        fi
    done
done
