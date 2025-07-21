#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N L2stats-SR
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=28

# load modules and go to workdir
module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
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
rm -f $logdir/cmd_L2subj_${PBS_JOBID}.txt
touch $logdir/cmd_L2subj_${PBS_JOBID}.txt

rm -f L2stats-SR.o*
rm -f L2stats-SR.e*

rm $logdir/re-runL2subj.log


type="act"               # "act" or "ppi" (or "nppi-dmn")
task=sharedreward       # edit if necessary
sm=5                    # smoothing kernel label
model=1                 # first-level model number
NCOPES=34               # base number of copes for act


for sub in ${subjects[@]}; do

    MAINOUTPUT=${projectdir}/derivatives/fsl/space-MNI/sub-${sub}
    
    # Initialize arrays to store all available inputs for this subject
    all_inputs=()
    input_labels=()
    
    # Collect all available L1 outputs across all sessions and runs
    for ses in {01..12}; do
        SESDIR=${MAINOUTPUT}/ses-${ses}
        
        # skip if the session folder itself doesn't exist
        if [ ! -d "${SESDIR}" ]; then
            echo "SKIP sub-${sub} ses-${ses}: session directory not found" >> $logdir/re-runL2subj.log
            continue
        fi
        
        # Check each run for this session
        for run in 1 2; do
            INPUT=${SESDIR}/L1_task-${task}_model-${model}_type-${type}_run-${run}_sm-${sm}.feat
            
            if [ -d "${INPUT}" ]; then
                all_inputs+=(${INPUT})
                input_labels+=("ses-${ses}_run-${run}")
            else
                echo "MISSING sub-${sub} ses-${ses} run-${run}: ${INPUT}" >> $logdir/re-runL2subj.log
            fi
        done
    done
    
    # Skip subject if fewer than 2 valid inputs found (need at least 2 for L2 analysis)
    if [ ${#all_inputs[@]} -lt 2 ]; then
        echo "SKIP sub-${sub}: insufficient inputs (${#all_inputs[@]} found), need at least 2" >> $logdir/re-runL2subj.log
        continue
    fi
    
    # Set output path for subject-level analysis
    OUTPUT=${MAINOUTPUT}/L2_task-${task}_model-${model}_type-${type}_subj-${sub}_sm-${sm}
    
    # skip if output already exists
    if [ -e ${OUTPUT}.gfeat/cope${NCOPES}.feat/cluster_mask_zstat1.nii.gz ]; then
        echo "SKIP sub-${sub}: L2 already done" >> $logdir/re-runL2subj.log
        continue
    fi
    
    # set template based on type
    if [ "${type}" == "act" ]; then
        ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_type-act_subject-level.fsf
        NCOPES=34
    else
        ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_type-nppi-dmn_subject-level.fsf
        NCOPES=$((34 + 1))
    fi
    
    # Create dynamic FSF template
    OTEMPLATE=${MAINOUTPUT}/L2_task-${task}_model-${model}_type-${type}_subj-${sub}.fsf
    
    # Start with base template and modify for dynamic inputs
    cp ${ITEMPLATE} ${OTEMPLATE}
    
    # Update basic parameters
    sed -i "s@OUTPUT@${OUTPUT}@g" ${OTEMPLATE}
    sed -i "s@set fmri(npts) .*@set fmri(npts) ${#all_inputs[@]}@g" ${OTEMPLATE}
        
    # Add inputs
    for i in "${!all_inputs[@]}"; do
        feat_index=$((i + 1))
        echo "set feat_files(${feat_index}) \"${all_inputs[i]}\"" >> ${OTEMPLATE}
        echo "set fmri(evg${feat_index}.1) 1.0" >> ${OTEMPLATE}
    done
    
    # Log what we're processing
    echo "PROCESSING sub-${sub}: ${#all_inputs[@]} inputs (${input_labels[*]})" >> $logdir/re-runL2subj.log
    
    # Add to command file
    echo feat $OTEMPLATE >>$logdir/cmd_L2_${PBS_JOBID}.txt

done

torque-launch -p "$logdir/chk_L2_${PBS_JOBID}.txt" "$logdir/cmd_L2_${PBS_JOBID}.txt"

# delete unused files
for sub in ${subjects[@]}; do
    OUTPUT=${projectdir}/derivatives/fsl/space-MNI/sub-${sub}/L2_task-${task}_model-${model}_type-${type}_subj-${sub}_sm-${sm}
    rm -f ${OUTPUT}.gfeat/cope*.feat/stats/res4d.nii.gz
    rm -f ${OUTPUT}.gfeat/cope*.feat/stats/corrections.nii.gz
    rm -f ${OUTPUT}.gfeat/cope*.feat/stats/threshac1.nii.gz
    rm -f ${OUTPUT}.gfeat/cope*.feat/filtered_func_data.nii.gz
    rm -f ${OUTPUT}.gfeat/cope*.feat/var_filtered_func_data.nii.gz
done