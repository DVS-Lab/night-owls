#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -N fmriprep-nightowls
#PBS -q large
#PBS -l nodes=1:ppn=14

# load modules and go to workdir.
module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
module load singularity
cd $PBS_O_WORKDIR

# ensure paths are correct
projectname=night-owls
maindir=/gpfs/scratch/tug87422/smithlab-shared/$projectname
scriptdir=$maindir/code
bidsdir=$maindir/bids
logdir=$maindir/logs
mkdir -p $logdir

#subjects=("${!1}")

rm -f $logdir/cmd_fmriprep_${PBS_JOBID}.txt
touch $logdir/cmd_fmriprep_${PBS_JOBID}.txt

# make derivatives folder if it doesn't exist.
# let's keep this out of bids for now
if [ ! -d $maindir/derivatives/anat-only ]; then
	mkdir -p $maindir/derivatives/anat-only
fi

scratchdir=~/scratch/$projectname/fmriprep-anat
if [ ! -d $scratchdir ]; then
	mkdir -p $scratchdir
fi

TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow
export SINGULARITYENV_MPLCONFIGDIR=/opt/mplconfigdir

for sub in ${subjects[@]}; do
		echo singularity run --cleanenv \
		-B ${TEMPLATEFLOW_DIR}:/opt/templateflow \
		-B ${MPLCONFIGDIR_DIR}:/opt/mplconfigdir \
		-B $maindir:/base \
		-B /gpfs/scratch/tug87422/smithlab-shared/tools/licenses:/opts \
		-B $scratchdir:/scratch \
		/gpfs/scratch/tug87422/smithlab-shared/tools/fmriprep-24.1.1.simg \
		/base/bids /base/derivatives/anat-only \
		participant --participant_label $sub \
		--stop-on-first-crash \
		--skip-bids-validation \
		--nthreads 14 \
		--me-output-echos \
		--output-spaces anat MNI152NLin6Asym \
        --anat-only \
        --longitudinal \
		--fs-no-reconall --fs-license-file /opts/fs_license.txt -w /scratch >> $logdir/cmd_fmriprep_${PBS_JOBID}.txt
done
torque-launch -p $logdir/chk_fmriprep_${PBS_JOBID}.txt $logdir/cmd_fmriprep_${PBS_JOBID}.txt

# --cifti-output 91k \
# --output-spaces fsLR fsaverage MNI152NLin6Asym \
