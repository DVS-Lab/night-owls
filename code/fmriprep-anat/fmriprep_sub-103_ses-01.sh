#!/bin/bash
#PBS -l walltime=4:00:00
#PBS -N fmriprep-sub-103-ses-01
#PBS -q normal
#PBS -l nodes=1:ppn=14

# load modules and go to workdir
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

TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow
export SINGULARITYENV_MPLCONFIGDIR=/opt/mplconfigdir

singularity run --cleanenv \
	-B ${TEMPLATEFLOW_DIR}:/opt/templateflow \
	-B ${MPLCONFIGDIR_DIR}:/opt/mplconfigdir \
	-B $maindir:/base \
	-B /gpfs/scratch/tug87422/smithlab-shared/tools/licenses:/opts \
	-B $maindir/scratch/anat-ses-01:/scratch \
	/gpfs/scratch/tug87422/smithlab-shared/tools/fmriprep-24.1.1.simg \
	/base/bids /base/derivatives/anat-only/ses-01 \
	participant --participant_label sub-103 \
	--stop-on-first-crash \
	--skip-bids-validation \
	--nthreads 14 \
	--me-output-echos \
	--output-spaces anat:res-2 MNI152NLin6Asym:res-2 \
  --derivatives /gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/anat-only \
	--bids-filter-file /base/code/fmriprep-anat/fmriprep_config_sub-103_ses-01.json \
	--fs-no-reconall --fs-license-file /opts/fs_license.txt \
	-w /scratch
