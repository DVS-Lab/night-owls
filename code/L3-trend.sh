
#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N L3stats-trend
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=28

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

rm -f $logdir/cmd_L3_sub-${sub}_trend_${PBS_JOBID}.txt
touch $logdir/cmd_L3_sub-${sub}_trend_${PBS_JOBID}.txt

type="act"               # "act" or "ppi" (or "nppi-dmn")
#sm=5                    # smoothing kernel label
model=1                 # first-level model number
tasks=("sharedreward" "mid")
#echos=(single-echo multi-echo)
#confounds=(cnfds-fmriprep cnfds-tedana)
space=t1w
echos=(multi-echo)
confounds=(cnfds-tedana)

for sub in ${subjects[@]}; do

    #Choose itemplate n
    if [ "$sub" == "101" ]; then
        nses=9
    elif [ "$sub" == "103" ]; then
        nses=11
    else 
        nses=12
    fi

    MAINOUTPUT=${projectdir}/derivatives/fsl/sub-${sub}
    rm -f L3stats_sub-${sub}_trend.o*
    rm -f L3stats_sub-${sub}_trend.e*

    rm -f $logdir/re-runL3_sub-${sub}.log

    for task in "${tasks[@]}"; do
        for echo in "${echos[@]}"; do
            for confound in "${confounds[@]}"; do

    copen=$([ "$task" = "mid" ] && echo 7 || echo 11)
    ITEMPLATE=${projectdir}/templates/L3_sessionTrend_n${nses}.fsf

    if [ "$sub" == "101" ]; then
        INPUT01=${MAINOUTPUT}/ses-01/L2_task-${task}_model-${model}_type-${type}_ses-01_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT02=${MAINOUTPUT}/ses-02/L2_task-${task}_model-${model}_type-${type}_ses-02_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT03=${MAINOUTPUT}/ses-03/L2_task-${task}_model-${model}_type-${type}_ses-03_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT04=${MAINOUTPUT}/ses-06/L2_task-${task}_model-${model}_type-${type}_ses-06_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT05=${MAINOUTPUT}/ses-07/L2_task-${task}_model-${model}_type-${type}_ses-07_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT06=${MAINOUTPUT}/ses-08/L2_task-${task}_model-${model}_type-${type}_ses-08_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT07=${MAINOUTPUT}/ses-09/L2_task-${task}_model-${model}_type-${type}_ses-09_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT08=${MAINOUTPUT}/ses-10/L2_task-${task}_model-${model}_type-${type}_ses-10_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT09=${MAINOUTPUT}/ses-11/L2_task-${task}_model-${model}_type-${type}_ses-11_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
    elif [ "$sub" == "103" ]; then
        INPUT01=${MAINOUTPUT}/ses-01/L2_task-${task}_model-${model}_type-${type}_ses-01_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT02=${MAINOUTPUT}/ses-02/L2_task-${task}_model-${model}_type-${type}_ses-02_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT03=${MAINOUTPUT}/ses-03/L2_task-${task}_model-${model}_type-${type}_ses-03_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT04=${MAINOUTPUT}/ses-04/L2_task-${task}_model-${model}_type-${type}_ses-04_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT05=${MAINOUTPUT}/ses-05/L2_task-${task}_model-${model}_type-${type}_ses-05_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT06=${MAINOUTPUT}/ses-06/L2_task-${task}_model-${model}_type-${type}_ses-06_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT07=${MAINOUTPUT}/ses-07/L2_task-${task}_model-${model}_type-${type}_ses-07_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT08=${MAINOUTPUT}/ses-08/L2_task-${task}_model-${model}_type-${type}_ses-08_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT09=${MAINOUTPUT}/ses-09/L2_task-${task}_model-${model}_type-${type}_ses-09_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT10=${MAINOUTPUT}/ses-10/L2_task-${task}_model-${model}_type-${type}_ses-10_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT11=${MAINOUTPUT}/ses-11/L2_task-${task}_model-${model}_type-${type}_ses-11_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
    else 
        INPUT01=${MAINOUTPUT}/ses-01/L2_task-${task}_model-${model}_type-${type}_ses-01_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT02=${MAINOUTPUT}/ses-02/L2_task-${task}_model-${model}_type-${type}_ses-02_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT03=${MAINOUTPUT}/ses-03/L2_task-${task}_model-${model}_type-${type}_ses-03_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT04=${MAINOUTPUT}/ses-04/L2_task-${task}_model-${model}_type-${type}_ses-04_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT05=${MAINOUTPUT}/ses-05/L2_task-${task}_model-${model}_type-${type}_ses-05_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT06=${MAINOUTPUT}/ses-06/L2_task-${task}_model-${model}_type-${type}_ses-06_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT07=${MAINOUTPUT}/ses-07/L2_task-${task}_model-${model}_type-${type}_ses-07_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT08=${MAINOUTPUT}/ses-08/L2_task-${task}_model-${model}_type-${type}_ses-08_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT09=${MAINOUTPUT}/ses-09/L2_task-${task}_model-${model}_type-${type}_ses-09_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT10=${MAINOUTPUT}/ses-10/L2_task-${task}_model-${model}_type-${type}_ses-10_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT11=${MAINOUTPUT}/ses-11/L2_task-${task}_model-${model}_type-${type}_ses-11_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
        INPUT12=${MAINOUTPUT}/ses-12/L2_task-${task}_model-${model}_type-${type}_ses-12_space-${space}_${echo}_${confound}.feat/cope${copen}.feat/stats/cope${copen}.nii.gz
    fi

    OUTPUT=${MAINOUTPUT}/subject-level/L3_task-${task}_trend_space-${space}_${echo}_${confound}
    #NCOPES=30

    NCOPES=$([ "$task" = "mid" ] && echo 10 || echo 15)

    # skip if output already exists
    if [ -e ${OUTPUT}.gfeat/cope${NCOPES}.feat/cope${copen}.feat/cluster_mask_zstat1.nii.gz ]; then
        echo "SKIP sub-${sub} ${task}: L3 already done" >> $logdir/re-runL3_sub-${sub}.log
        continue
    fi

    # build and run session FSF
    OTEMPLATE=${MAINOUTPUT}/L3_task-${task}_trend_space-${space}_${echo}_${confound}.fsf

    if [ "$sub" == "101" ]; then
        sed -e 's@OUTPUT@'$OUTPUT'@g' \
        -e 's@INPUT01@'$INPUT01'@g' \
        -e 's@INPUT02@'$INPUT02'@g' \
        -e 's@INPUT03@'$INPUT03'@g' \
        -e 's@INPUT04@'$INPUT04'@g' \
        -e 's@INPUT05@'$INPUT05'@g' \
        -e 's@INPUT06@'$INPUT06'@g' \
        -e 's@INPUT07@'$INPUT07'@g' \
        -e 's@INPUT08@'$INPUT08'@g' \
        -e 's@INPUT09@'$INPUT09'@g' \
        <$ITEMPLATE >$OTEMPLATE
        echo feat $OTEMPLATE >>$logdir/cmd_L3_sub-${sub}_trend_${PBS_JOBID}.txt
    elif [ "$sub" == "103" ]; then
        sed -e 's@OUTPUT@'$OUTPUT'@g' \
        -e 's@INPUT01@'$INPUT01'@g' \
        -e 's@INPUT02@'$INPUT02'@g' \
        -e 's@INPUT03@'$INPUT03'@g' \
        -e 's@INPUT04@'$INPUT04'@g' \
        -e 's@INPUT05@'$INPUT05'@g' \
        -e 's@INPUT06@'$INPUT06'@g' \
        -e 's@INPUT07@'$INPUT07'@g' \
        -e 's@INPUT08@'$INPUT08'@g' \
        -e 's@INPUT09@'$INPUT09'@g' \
        -e 's@INPUT10@'$INPUT10'@g' \
        -e 's@INPUT11@'$INPUT11'@g' \
        <$ITEMPLATE >$OTEMPLATE
        echo feat $OTEMPLATE >>$logdir/cmd_L3_sub-${sub}_trend_${PBS_JOBID}.txt
    else 
        sed -e 's@OUTPUT@'$OUTPUT'@g' \
        -e 's@INPUT01@'$INPUT01'@g' \
        -e 's@INPUT02@'$INPUT02'@g' \
        -e 's@INPUT03@'$INPUT03'@g' \
        -e 's@INPUT04@'$INPUT04'@g' \
        -e 's@INPUT05@'$INPUT05'@g' \
        -e 's@INPUT06@'$INPUT06'@g' \
        -e 's@INPUT07@'$INPUT07'@g' \
        -e 's@INPUT08@'$INPUT08'@g' \
        -e 's@INPUT09@'$INPUT09'@g' \
        -e 's@INPUT10@'$INPUT10'@g' \
        -e 's@INPUT11@'$INPUT11'@g' \
        -e 's@INPUT12@'$INPUT12'@g' \
        <$ITEMPLATE >$OTEMPLATE
        echo feat $OTEMPLATE >>$logdir/cmd_L3_sub-${sub}_trend_${PBS_JOBID}.txt
    fi

    done
    done
done
done

torque-launch -p "$logdir/chk_L3_sub-${sub}_trend_${PBS_JOBID}.txt" "$logdir/cmd_L3_sub-${sub}_trend_${PBS_JOBID}.txt"
