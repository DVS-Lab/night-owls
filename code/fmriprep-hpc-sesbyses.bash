#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -N fmriprep-nightowls
#PBS -q large
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

#Keep scratch and fmriprep dirs separate by sessions for now
for ses in {01..12}; do
  if [ ! -d "$maindir/scratch/ses-$ses" ]; then
    mkdir -p "$maindir/scratch/ses-$ses"
  fi
   if [ ! -d "$maindir/fmriprep/ses-$ses" ]; then
    mkdir -p "$maindir/fmriprep/ses-$ses"
  fi
done


rm -f $logdir/cmd_fmriprep_${PBS_JOBID}.txt
touch $logdir/cmd_fmriprep_${PBS_JOBID}.txt



TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow
export SINGULARITYENV_MPLCONFIGDIR=/opt/mplconfigdir

for sub in ${subjects[@]}; do
	
	sessions=($(find "${bidsdir}/${sub}" -type f -path "*/func/*.nii.gz" | grep -o "ses-[^/]*" | sort -u))
  
  for ses in "${sessions[@]}"; do
    sesid=${ses#ses-}
		
    # Create session-specific BIDS filter file
    configfile=$maindir/code/fmriprep_config_${sub}_${ses}.json
    cat <<EOF > "$configfile"
{
  "t1w": {"datatype": "anat", "suffix": "T1w", "session": "$sesid"},
  "sbref": {"datatype": "func", "suffix": "sbref", "part": [null, "mag"], "session": "$sesid"},
  "bold": {"datatype": "func", "suffix": "bold", "part": [null, "mag"], "session": "$sesid"}
}
EOF

	echo "singularity run --cleanenv \
		-B ${TEMPLATEFLOW_DIR}:/opt/templateflow \
		-B ${MPLCONFIGDIR_DIR}:/opt/mplconfigdir \
		-B $maindir:/base \
		-B /gpfs/scratch/tug87422/smithlab-shared/tools/licenses:/opts \
		-B $scratchdir:/scratch \ #update this
		/gpfs/scratch/tug87422/smithlab-shared/tools/fmriprep-24.1.1.simg \
		/base/bids /base/derivatives/fmriprep \
		participant --participant_label $sub \
		--stop-on-first-crash \
		--skip-bids-validation \
		--nthreads 14 \
		--me-output-echos \
		--output-spaces MNI152NLin6Asym \
		--bids-filter-file /base/code/fmriprep_config_${sub}_${ses}.json \
		--fs-no-reconall --fs-license-file /opts/fs_license.txt \
		-w /scratch \
		>> "$logdir/cmd_fmriprep_${PBS_JOBID}.txt"
	done
done