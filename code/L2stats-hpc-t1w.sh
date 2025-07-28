
#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N L2stats-t1w
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=14

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

rm -f $logdir/cmd_L2-t1w_${PBS_JOBID}.txt
touch $logdir/cmd_L2-t1w_${PBS_JOBID}.txt

type="act"               # "act" or "ppi" (or "nppi-dmn")
sm=5                    # smoothing kernel label
model=1                 # first-level model number
tasks=("sharedreward" "mid")

for sub in ${subjects[@]}; do

    MAINOUTPUT=${projectdir}/derivatives/fsl/space-t1w/sub-${sub}
    rm -f L2stats-t1w_sub-${sub}.o*
    rm -f L2stats-t1w_sub-${sub}.e*

    rm -f $logdir/re-runL2_sub-${sub}-t1w.log

    for ses in {01..12}; do

        SESDIR=${MAINOUTPUT}/ses-${ses}

        for task in "${tasks[@]}"; do

            # skip if the session folder itself doesn’t exist
            if [ ! -d "${SESDIR}" ]; then
                echo "SKIP sub-${sub} ses-${ses} ${task}: session directory not found" >> $logdir/re-runL2_sub-${sub}-t1w.log
                continue
            fi

            ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_2-runs.fsf
            INPUT1=${SESDIR}/L1_task-${task}_model-${model}_type-${type}_run-1_sm-${sm}.feat
            INPUT2=${SESDIR}/L1_task-${task}_model-${model}_type-${type}_run-2_sm-${sm}.feat
            OUTPUT=${SESDIR}/L2_task-${task}_model-${model}_type-${type}_ses-${ses}_sm-${sm}
            #NCOPES=30

            # skip if either run folder is missing
            missing=()
                [ ! -d "${INPUT1}" ] && missing+=(run-1)
                [ ! -d "${INPUT2}" ] && missing+=(run-2)
            if [ ${#missing[@]} -gt 0 ]; then
                echo "SKIP sub-${sub} ses-${ses} ${task}: missing ${missing[*]}" >> $logdir/re-runL2_sub-${sub}-t1w.log
                continue
            fi

            # skip if output already exists
            #if [ -e ${OUTPUT}.gfeat/cope${NCOPES}.feat/cluster_mask_zstat1.nii.gz ]; then
            #    echo "SKIP sub-${sub} ses-${ses} ${task}: L2 already done" >> $logdir/re-runL2_sub-${sub}-t1w.log
            #    continue
            #fi

            # build and run session FSF
            OTEMPLATE=${SESDIR}/L2_task-${task}_model-${model}_type-${type}.fsf
            sed -e 's@OUTPUT@'$OUTPUT'@g' \
            -e 's@INPUT1@'$INPUT1'@g' \
            -e 's@INPUT2@'$INPUT2'@g' \
            <$ITEMPLATE >$OTEMPLATE
            echo feat $OTEMPLATE >>$logdir/cmd_L2_${PBS_JOBID}.txt

        done
    done
done

torque-launch -p "$logdir/chk_L2-t1w_${PBS_JOBID}.txt" "$logdir/cmd_L2-t1w_${PBS_JOBID}.txt"

# delete unused files
#for sub in ${subjects[@]}; do
#    for ses in {01..12}; do
#        OUTPUT=${projectdir}/derivatives/fsl/space-MNI/sub-${sub}/L2_task-${task}_model-${model}_type-${type}_ses-${ses}_sm-${sm}
#        rm -f ${OUTPUT}.gfeat/cope*.feat/stats/res4d.nii.gz
#        rm -f ${OUTPUT}.gfeat/cope*.feat/stats/corrections.nii.gz
#        rm -f ${OUTPUT}.gfeat/cope*.feat/stats/threshac1.nii.gz
#        rm -f ${OUTPUT}.gfeat/cope*.feat/filtered_func_data.nii.gz
#        rm -f ${OUTPUT}.gfeat/cope*.feat/var_filtered_func_data.nii.gz
#    done
#done