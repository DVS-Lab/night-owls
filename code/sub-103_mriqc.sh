#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -N mriqc103test
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -l mem=100gb
#PBS -o /home/tun47039/mriqc_103.o
#PBS -e /home/tun47039/mriqc_103.e

module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
module load singularity
cd $PBS_O_WORKDIR

umask 0000

dsroot=/gpfs/scratch/tug87422/smithlab-shared/night-owls
codedir=$dsroot/code
logdir=~/work/logs
mkdir -p $logdir

subjects=(103)

rm -f $logdir/cmd_mriqc_${PBS_JOBID}.txt
touch $logdir/cmd_mriqc_${PBS_JOBID}.txt

if [ ! -d $dsroot/derivatives/mriqc ]; then
        mkdir -p $dsroot/derivatives/mriqc
fi

scratch=$dsroot/scratch
if [ ! -d $scratch ]; then
	mkdir -p $scratch
fi

TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export APPTAINERENV_TEMPLATEFLOW_HOME=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
export APPTAINERENV_MPLCONFIGDIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir

export job_scratch=$scratch/mriqc_${PBS_JOBID}
mkdir -p $job_scratch

for sub in ${subjects[@]}; do
        echo singularity run --cleanenv \
        -B ${TEMPLATEFLOW_DIR}:/templateflow \
	-B ${MPLCONFIGDIR_DIR}:/mplconfigdir \
        -B $dsroot/bids:/data \
        -B $dsroot/derivatives/mriqc:/out \
        -B $job_scratch:/workdir \
        /gpfs/scratch/tug87422/smithlab-shared/tools/mriqc-24.0.2.simg \
        /data /out participant \
	--participant_label $sub \
	--no-datalad-get \
	--no-sub \
	--modalities T1w \
        -w /workdir >> $logdir/cmd_mriqc_${PBS_JOBID}.txt
done

torque-launch -p $logdir/chk_mriqc_${PBS_JOBID}.txt $logdir/cmd_mriqc_${PBS_JOBID}.txt