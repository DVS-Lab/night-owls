#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -N mriqc_103
#PBS -q normal
#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb

module load singularity
cd $PBS_O_WORKDIR

dsroot=/gpfs/scratch/tug87422/smithlab-shared/night-owls
scratch=$dsroot/scratch
job_scratch=$scratch/mriqc_103_${PBS_JOBID}
mkdir -p $job_scratch

TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export APPTAINERENV_TEMPLATEFLOW_HOME=$TEMPLATEFLOW_DIR
export APPTAINERENV_MPLCONFIGDIR=$MPLCONFIGDIR_DIR

singularity run --cleanenv \
  -B ${TEMPLATEFLOW_DIR}:/templateflow \
  -B ${MPLCONFIGDIR_DIR}:/mplconfigdir \
  -B $dsroot/bids:/data \
  -B $dsroot/derivatives/mriqc:/out \
  -B $job_scratch:/workdir \
  /gpfs/scratch/tug87422/smithlab-shared/tools/mriqc-24.0.2.simg \
  /data /out participant \
  --participant_label 103 \
  --no-datalad-get \
  --no-sub \
  --modalities T1w \
  -w /workdir