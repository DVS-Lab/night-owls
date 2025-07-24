
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

rm -f L2stats-SR.o*
rm -f L2stats-SR.e*

rm $logdir/re-runL2.log

type="act"               # "act" or "ppi" (or "nppi-dmn")
task=sharedreward       # edit if necessary
sm=5                    # smoothing kernel label
model=1                 # first-level model number

for sub in ${subjects[@]}; do

    MAINOUTPUT=${projectdir}/derivatives/fsl/space-MNI/sub-${sub}

    for ses in {01..12}; do
        SESDIR=${MAINOUTPUT}/ses-${ses}

        # skip if the session folder itself doesn’t exist
        if [ ! -d "${SESDIR}" ]; then
            echo "SKIP sub-${sub} ses-${ses}: session directory not found" >> $logdir/re-runL2.log
            continue
        fi

        # set template & cope count
        if [ "${type}" == "act" ]; then
            ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_2runs.fsf
            NCOPES=34
        else
            ITEMPLATE=${projectdir}/templates/L2_task-${task}_model-${model}_type-nppi-dmn.fsf
            NCOPES=$((34 + 1))
        fi

        INPUT1=${SESDIR}/L1_task-${task}_model-${model}_type-${type}_run-1_sm-${sm}.feat
        INPUT2=${SESDIR}/L1_task-${task}_model-${model}_type-${type}_run-2_sm-${sm}.feat
        OUTPUT=${SESDIR}/L2_task-${task}_model-${model}_type-${type}_ses-${ses}_sm-${sm}

        # skip if either run folder is missing
        missing=()
            [ ! -d "${INPUT1}" ] && missing+=(run-1)
            [ ! -d "${INPUT2}" ] && missing+=(run-2)
        if [ ${#missing[@]} -gt 0 ]; then
            echo "SKIP sub-${sub} ses-${ses}: missing ${missing[*]}" >> $logdir/re-runL2.log
            continue
        fi

        # skip if output already exists
        if [ -e ${OUTPUT}.gfeat/cope${NCOPES}.feat/cluster_mask_zstat1.nii.gz ]; then
            echo "SKIP sub-${sub} ses-${ses}: L2 already done" >> $logdir/re-runL2.log
            continue
        fi

        # build and run session FSF
        OTEMPLATE=${SESDIR}/L2_task-${task}_model-${model}_type-${type}.fsf
        sed -e 's@OUTPUT@'$OUTPUT'@g' \
        -e 's@INPUT1@'$INPUT1'@g' \
        -e 's@INPUT2@'$INPUT2'@g' \
        <$ITEMPLATE >$OTEMPLATE
        echo feat $OTEMPLATE >>$logdir/cmd_L2_${PBS_JOBID}.txt


    done
done

torque-launch -p "$logdir/chk_L2_${PBS_JOBID}.txt" "$logdir/cmd_L2_${PBS_JOBID}.txt"

# delete unused files
for sub in ${subjects[@]}; do
    for ses in {01..12}; do
        OUTPUT=${projectdir}/derivatives/fsl/space-MNI/sub-${sub}/ses-${ses}/L2_task-${task}_model-${model}_type-${type}_ses-${ses}_sm-${sm}
        rm -f ${OUTPUT}.gfeat/cope*.feat/stats/res4d.nii.gz
        rm -f ${OUTPUT}.gfeat/cope*.feat/stats/corrections.nii.gz
        rm -f ${OUTPUT}.gfeat/cope*.feat/stats/threshac1.nii.gz
        rm -f ${OUTPUT}.gfeat/cope*.feat/filtered_func_data.nii.gz
        rm -f ${OUTPUT}.gfeat/cope*.feat/var_filtered_func_data.nii.gz
    done
done