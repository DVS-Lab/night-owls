#/usr/bin/env bash
#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -N mriqc
#PBS -q normal
#PBS -l nodes=1:ppn=28

# load modules and go to workdir
module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
module load singularity
cd $PBS_O_WORKDIR

umask 0000

# ensure paths are correct irrespective from where user runs the script
dsroot=/gpfs/scratch/tug87422/smithlab-shared/night-owls
codedir=$dsroot/code
logdir=~/work/logs
mkdir -p $logdir

#subjects=("${!1}")

rm -f $logdir/cmd_mriqc_${PBS_JOBID}.txt
touch $logdir/cmd_mriqc_${PBS_JOBID}.txt

# Delete old output files
rm -f $codedir/mriqc-test.o*
rm -f $codedir/mriqc-test.e*

## Run MRIQC on subject
# To-do: run through datalad with YODA principles

## make derivatives folder if it doesn't exist.
## let's keep this out of bids for now

if [ ! -d $dsroot/derivatives/mriqc ]; then
        mkdir -p $dsroot/derivatives/mriqc
fi

# make scratch
scratch=$dsroot/scratch
if [ ! -d $scratch ]; then
	mkdir -p $scratch
fi

TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export APPTAINERENV_TEMPLATEFLOW_HOME=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
export APPTAINERENV_MPLCONFIGDIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir

# no space left on device error for v0.15.2 and higher
# https://neurostars.org/t/mriqc-no-space-left-on-device-error/16187/1
# https://github.com/poldracklab/mriqc/issues/850

# Job writable scratch to deal with permissions issues
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
	--modalities bold \
        -w /workdir >> $logdir/cmd_mriqc_${PBS_JOBID}.txt
done

torque-launch -p $logdir/chk_mriqc_${PBS_JOBID}.txt $logdir/cmd_mriqc_${PBS_JOBID}.txt

# --bids-filter-file $codedir/fmriprep_config.json
