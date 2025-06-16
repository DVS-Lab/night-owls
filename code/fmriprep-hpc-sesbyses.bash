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


scratchdir=$maindir/scratch 
#status file
status_file=$scratchdir/status.txt
echo "# sub ses status" > "$status_file"


#subjects=("${!1}")

rm -f $logdir/cmd_fmriprep_${PBS_JOBID}.txt
touch $logdir/cmd_fmriprep_${PBS_JOBID}.txt

# make derivatives folder if it doesn't exist.
if [ ! -d $maindir/derivatives ]; then
	mkdir -p $maindir/derivatives
fi
# make derivatives/fmriprep parent folder if it doesn't exist
if [ ! -d $maindir/derivatives/fmriprep ]; then
    mkdir -p $maindir/derivatives/fmriprep
fi

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
		-B $scratchdir:/scratch \
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
		-w /scratch && \
	if ls $maindir/derivatives/fmriprep/$sub/$ses/func/${sub}_${ses}_task-*_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz 1> /dev/null 2>&1; then
		echo \"$sub $ses fmriprep successful\" >> \"$status_file\"; \
	else \
		echo \"$sub $ses fmriprep unsuccessful\" >> \"$status_file\"; \
		exit 1; \
	fi" >> "$logdir/cmd_fmriprep_${PBS_JOBID}.txt"
	
	done
done
torque-launch -p $logdir/chk_fmriprep_${PBS_JOBID}.txt $logdir/cmd_fmriprep_${PBS_JOBID}.txt

# --cifti-output 91k \
# --output-spaces fsLR fsaverage MNI152NLin6Asym \
