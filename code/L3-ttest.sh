
#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N L3stats
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

rm -f $logdir/cmd_L3_sub-${sub}_${PBS_JOBID}.txt
touch $logdir/cmd_L3_sub-${sub}_${PBS_JOBID}.txt

type="act"               # "act" or "ppi" (or "nppi-dmn")
#sm=5                    # smoothing kernel label
model=1                 # first-level model number
tasks=("sharedreward" "mid")
#echos=(single-echo multi-echo)
#confounds=(cnfds-fmriprep cnfds-tedana)
space="t1w"
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
    rm -f L3stats_sub-${sub}.o*
    rm -f L3stats_sub-${sub}.e*

    rm -f $logdir/re-runL3_sub-${sub}.log

    for task in "${tasks[@]}"; do
        for echo in "${echos[@]}"; do
            for confound in "${confounds[@]}"; do

    ITEMPLATE=${projectdir}/templates/L3_PairedTTest_n${nses}.fsf

    if [ "$sub" == "101" ]; then
        INPUT01=${MAINOUTPUT}/ses-01/L1_sub-${sub}_ses-01_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT02=${MAINOUTPUT}/ses-02/L1_sub-${sub}_ses-02_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT03=${MAINOUTPUT}/ses-03/L1_sub-${sub}_ses-03_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT04=${MAINOUTPUT}/ses-06/L1_sub-${sub}_ses-07_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT05=${MAINOUTPUT}/ses-07/L1_sub-${sub}_ses-07_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT06=${MAINOUTPUT}/ses-08/L1_sub-${sub}_ses-09_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT07=${MAINOUTPUT}/ses-09/L1_sub-${sub}_ses-09_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT08=${MAINOUTPUT}/ses-10/L1_sub-${sub}_ses-10_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT09=${MAINOUTPUT}/ses-11/L1_sub-${sub}_ses-11_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT10=${MAINOUTPUT}/ses-01/L1_sub-${sub}_ses-01_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT11=${MAINOUTPUT}/ses-02/L1_sub-${sub}_ses-02_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT12=${MAINOUTPUT}/ses-03/L1_sub-${sub}_ses-03_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT13=${MAINOUTPUT}/ses-06/L1_sub-${sub}_ses-06_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT14=${MAINOUTPUT}/ses-07/L1_sub-${sub}_ses-07_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT15=${MAINOUTPUT}/ses-08/L1_sub-${sub}_ses-08_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT16=${MAINOUTPUT}/ses-09/L1_sub-${sub}_ses-09_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT17=${MAINOUTPUT}/ses-10/L1_sub-${sub}_ses-10_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT18=${MAINOUTPUT}/ses-11/L1_sub-${sub}_ses-11_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
    elif [ "$sub" == "103" ]; then
        INPUT01=${MAINOUTPUT}/ses-01/L1_sub-${sub}_ses-01_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT02=${MAINOUTPUT}/ses-02/L1_sub-${sub}_ses-02_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT03=${MAINOUTPUT}/ses-03/L1_sub-${sub}_ses-03_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT04=${MAINOUTPUT}/ses-04/L1_sub-${sub}_ses-04_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT05=${MAINOUTPUT}/ses-05/L1_sub-${sub}_ses-05_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT06=${MAINOUTPUT}/ses-06/L1_sub-${sub}_ses-06_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT07=${MAINOUTPUT}/ses-07/L1_sub-${sub}_ses-07_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT08=${MAINOUTPUT}/ses-08/L1_sub-${sub}_ses-08_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT09=${MAINOUTPUT}/ses-09/L1_sub-${sub}_ses-09_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT10=${MAINOUTPUT}/ses-10/L1_sub-${sub}_ses-10_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT11=${MAINOUTPUT}/ses-11/L1_sub-${sub}_ses-11_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT12=${MAINOUTPUT}/ses-01/L1_sub-${sub}_ses-01_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT13=${MAINOUTPUT}/ses-02/L1_sub-${sub}_ses-02_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT14=${MAINOUTPUT}/ses-03/L1_sub-${sub}_ses-03_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT15=${MAINOUTPUT}/ses-04/L1_sub-${sub}_ses-04_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT16=${MAINOUTPUT}/ses-05/L1_sub-${sub}_ses-05_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT17=${MAINOUTPUT}/ses-06/L1_sub-${sub}_ses-06_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT18=${MAINOUTPUT}/ses-07/L1_sub-${sub}_ses-07_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT19=${MAINOUTPUT}/ses-08/L1_sub-${sub}_ses-08_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT20=${MAINOUTPUT}/ses-09/L1_sub-${sub}_ses-09_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT21=${MAINOUTPUT}/ses-10/L1_sub-${sub}_ses-10_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT22=${MAINOUTPUT}/ses-11/L1_sub-${sub}_ses-11_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
    else 
        INPUT01=${MAINOUTPUT}/ses-01/L1_sub-${sub}_ses-01_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT02=${MAINOUTPUT}/ses-02/L1_sub-${sub}_ses-02_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT03=${MAINOUTPUT}/ses-03/L1_sub-${sub}_ses-03_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT04=${MAINOUTPUT}/ses-04/L1_sub-${sub}_ses-04_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT05=${MAINOUTPUT}/ses-05/L1_sub-${sub}_ses-05_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT06=${MAINOUTPUT}/ses-06/L1_sub-${sub}_ses-06_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT07=${MAINOUTPUT}/ses-07/L1_sub-${sub}_ses-07_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT08=${MAINOUTPUT}/ses-08/L1_sub-${sub}_ses-08_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT09=${MAINOUTPUT}/ses-09/L1_sub-${sub}_ses-09_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT10=${MAINOUTPUT}/ses-10/L1_sub-${sub}_ses-10_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT11=${MAINOUTPUT}/ses-11/L1_sub-${sub}_ses-11_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT12=${MAINOUTPUT}/ses-12/L1_sub-${sub}_ses-12_task-${task}_model-${model}_type-${type}_run-1_space-${space}_${echo}_${confound}.feat
        INPUT13=${MAINOUTPUT}/ses-01/L1_sub-${sub}_ses-01_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT14=${MAINOUTPUT}/ses-02/L1_sub-${sub}_ses-02_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT15=${MAINOUTPUT}/ses-03/L1_sub-${sub}_ses-03_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT16=${MAINOUTPUT}/ses-04/L1_sub-${sub}_ses-04_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT17=${MAINOUTPUT}/ses-05/L1_sub-${sub}_ses-05_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT18=${MAINOUTPUT}/ses-06/L1_sub-${sub}_ses-06_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT19=${MAINOUTPUT}/ses-07/L1_sub-${sub}_ses-07_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT20=${MAINOUTPUT}/ses-08/L1_sub-${sub}_ses-08_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT21=${MAINOUTPUT}/ses-09/L1_sub-${sub}_ses-09_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT22=${MAINOUTPUT}/ses-10/L1_sub-${sub}_ses-10_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT23=${MAINOUTPUT}/ses-11/L1_sub-${sub}_ses-11_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
        INPUT24=${MAINOUTPUT}/ses-12/L1_sub-${sub}_ses-12_task-${task}_model-${model}_type-${type}_run-2_space-${space}_${echo}_${confound}.feat
    fi

    OUTPUT=${MAINOUTPUT}/L3_task-${task}_ttest_space-${space}_${echo}_${confound}
    #NCOPES=30

    NCOPES=$([ "$task" = "mid" ] && echo 10 || echo 15)

    # skip if output already exists
    if [ -e ${OUTPUT}.gfeat/cope${NCOPES}.feat/cluster_mask_zstat1.nii.gz ]; then
        echo "SKIP sub-${sub} ${task}: L3 already done" >> $logdir/re-runL3_sub-${sub}.log
        continue
    fi

    # build and run session FSF
    OTEMPLATE=${MAINOUTPUT}/L3_task-${task}_ttest_space-${space}_${echo}_${confound}.fsf

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
        -e 's@INPUT10@'$INPUT10'@g' \
        -e 's@INPUT11@'$INPUT11'@g' \
        -e 's@INPUT12@'$INPUT12'@g' \
        -e 's@INPUT13@'$INPUT13'@g' \
        -e 's@INPUT14@'$INPUT14'@g' \
        -e 's@INPUT15@'$INPUT15'@g' \
        -e 's@INPUT16@'$INPUT16'@g' \
        -e 's@INPUT17@'$INPUT17'@g' \
        -e 's@INPUT18@'$INPUT18'@g' \
        <$ITEMPLATE >$OTEMPLATE
        echo feat $OTEMPLATE >>$logdir/cmd_L3_sub-${sub}_${PBS_JOBID}.txt
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
        -e 's@INPUT12@'$INPUT12'@g' \
        -e 's@INPUT13@'$INPUT13'@g' \
        -e 's@INPUT14@'$INPUT14'@g' \
        -e 's@INPUT15@'$INPUT15'@g' \
        -e 's@INPUT16@'$INPUT16'@g' \
        -e 's@INPUT17@'$INPUT17'@g' \
        -e 's@INPUT18@'$INPUT18'@g' \
        -e 's@INPUT19@'$INPUT19'@g' \
        -e 's@INPUT20@'$INPUT20'@g' \
        -e 's@INPUT21@'$INPUT21'@g' \
        -e 's@INPUT22@'$INPUT22'@g' \
        <$ITEMPLATE >$OTEMPLATE
        echo feat $OTEMPLATE >>$logdir/cmd_L3_sub-${sub}_${PBS_JOBID}.txt
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
        -e 's@INPUT13@'$INPUT13'@g' \
        -e 's@INPUT14@'$INPUT14'@g' \
        -e 's@INPUT15@'$INPUT15'@g' \
        -e 's@INPUT16@'$INPUT16'@g' \
        -e 's@INPUT17@'$INPUT17'@g' \
        -e 's@INPUT18@'$INPUT18'@g' \
        -e 's@INPUT19@'$INPUT19'@g' \
        -e 's@INPUT20@'$INPUT20'@g' \
        -e 's@INPUT21@'$INPUT21'@g' \
        -e 's@INPUT22@'$INPUT22'@g' \
        -e 's@INPUT23@'$INPUT23'@g' \
        -e 's@INPUT24@'$INPUT24'@g' \
        <$ITEMPLATE >$OTEMPLATE
        echo feat $OTEMPLATE >>$logdir/cmd_L3_sub-${sub}_${PBS_JOBID}.txt
    fi

    done
    done
done
done

torque-launch -p "$logdir/chk_L3_${PBS_JOBID}.txt" "$logdir/cmd_L3_sub-${sub}_${PBS_JOBID}.txt"
