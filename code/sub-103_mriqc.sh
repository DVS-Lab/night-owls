#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -N mriqc103
#PBS -q normal
#PBS -l nodes=1:ppn=28
#PBS -l mem=100gb
#PBS -o /gpfs/scratch/tun47039/night-owls/logs/mriqc_103.o
#PBS -e /gpfs/scratch/tun47039/night-owls/logs/mriqc_103.e

module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
module load singularity
cd $PBS_O_WORKDIR

umask 0000

# New base directory
mydir=/gpfs/scratch/tun47039/night-owls
dsroot=/gpfs/scratch/tug87422/smithlab-shared/night-owls

# Create necessary directories
mkdir -p $mydir/logs
mkdir -p $mydir/scratch
mkdir -p $mydir/derivatives/mriqc

subjects=(103)

rm -f $mydir/logs/cmd_mriqc_${PBS_JOBID}.txt
touch $mydir/logs/cmd_mriqc_${PBS_JOBID}.txt

TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export APPTAINERENV_TEMPLATEFLOW_HOME=$TEMPLATEFLOW_DIR
export APPTAINERENV_MPLCONFIGDIR=$MPLCONFIGDIR_DIR

export job_scratch=$mydir/scratch/mriqc_${PBS_JOBID}
mkdir -p $job_scratch

for sub in ${subjects[@]}; do
        echo singularity run --cleanenv \
        -B ${TEMPLATEFLOW_DIR}:/templateflow \
	-B ${MPLCONFIGDIR_DIR}:/mplconfigdir \
        -B $dsroot/bids:/data \
        -B $mydir/derivatives/mriqc:/out \
        -B $job_scratch:/workdir \
        /gpfs/scratch/tug87422/smithlab-shared/tools/mriqc-24.0.2.simg \
        /data /out participant \
	--participant_label $sub \
	--no-datalad-get \
	--no-sub \
	--modalities T1w \
        -w /workdir >> $mydir/logs/cmd_mriqc_${PBS_JOBID}.txt
done

torque-launch -p $mydir/logs/chk_mriqc_${PBS_JOBID}.txt $mydir/logs/cmd_mriqc_${PBS_JOBID}.txt