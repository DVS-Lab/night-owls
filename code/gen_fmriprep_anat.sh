#!/bin/bash

umask 0000
# ensure paths are correct
projectname=night-owls 
maindir=/gpfs/scratch/tug87422/smithlab-shared/$projectname
scriptdir=$maindir/code
bidsdir=$maindir/bids
logdir=$maindir/logs
mkdir -p $logdir

#Keep scratch and derivatives dirs separate by sessions
for ses in {01..12}; do
  if [ ! -d "$maindir/scratch/anat-ses-$ses" ]; then
    mkdir -p "$maindir/scratch/anat-ses-$ses"
  fi
   if [ ! -d "$maindir/derivatives/anat-only/ses-$ses" ]; then
    mkdir -p "$maindir/derivatives/anat-only/ses-$ses"
  fi
   if [ ! -d "$maindir/code/fmriprep-anat" ]; then
    mkdir -p "$maindir/code/fmriprep-anat"
  fi
done

# Generate configs and PBS scripts

mapfile -t subjects < "sublist.txt"

for sub in "${subjects[@]}"; do
  sub="sub-$sub"
    ## for sessions that have a functional image..
    sessions=( $(
        find "$bidsdir/$sub" \
            -type f \
            -path "*/ses-[0-9][0-9]/func/*.nii.gz" \
            2>/dev/null \
        | sed -n 's|.*/\(ses-[0-9][0-9]\)/.*|\1|p' \
        | sort -u
    ) )


  for ses in "${sessions[@]}"; do
    sesid=${ses#ses-}
		
    ## Create subject/session-specific BIDS filter file
    configfile=$maindir/code/fmriprep-anat/fmriprep_config_${sub}_${ses}.json
    cat <<EOF > "$configfile"
{
  "sbref": {"datatype": "func", "suffix": "sbref", "part": [null, "mag"], "session": ["$sesid"]},
  "bold": {"datatype": "func", "suffix": "bold", "part": [null, "mag"], "session": ["$sesid"]}
}
EOF
 # "t1w": {"datatype": "anat", "suffix": "T1w", "session": ["$sesid"]},

    # Generate subject/session-specific PBS script
    scriptfile=$maindir/code/fmriprep-anat/fmriprep_${sub}_${ses}.sh

    cat <<EOF > "$scriptfile"
#!/bin/bash
#PBS -l walltime=4:00:00
#PBS -N fmriprep-${sub}-${ses}
#PBS -q normal
#PBS -l nodes=1:ppn=14

# load modules and go to workdir
module load fsl/6.0.2
source \$FSLDIR/etc/fslconf/fsl.sh
module load singularity
cd \$PBS_O_WORKDIR

# ensure paths are correct
projectname=night-owls 
maindir=/gpfs/scratch/tug87422/smithlab-shared/\$projectname
scriptdir=\$maindir/code
bidsdir=\$maindir/bids
logdir=\$maindir/logs
mkdir -p \$logdir

TEMPLATEFLOW_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/templateflow
MPLCONFIGDIR_DIR=/gpfs/scratch/tug87422/smithlab-shared/tools/mplconfigdir
export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow
export SINGULARITYENV_MPLCONFIGDIR=/opt/mplconfigdir

singularity run --cleanenv \\
	-B \${TEMPLATEFLOW_DIR}:/opt/templateflow \\
	-B \${MPLCONFIGDIR_DIR}:/opt/mplconfigdir \\
	-B \$maindir:/base \\
	-B /gpfs/scratch/tug87422/smithlab-shared/tools/licenses:/opts \\
	-B \$maindir/scratch/anat-$ses:/scratch \\
	/gpfs/scratch/tug87422/smithlab-shared/tools/fmriprep-24.1.1.simg \\
	/base/bids /base/derivatives/anat-only/$ses \\
	participant --participant_label $sub \\
	--stop-on-first-crash \\
	--skip-bids-validation \\
	--nthreads 14 \\
	--me-output-echos \\
	--output-spaces anat:res-2 MNI152NLin6Asym:res-2 \\
  --derivatives $maindir/derivatives/anat-only \\
	--bids-filter-file /base/code/fmriprep-anat/fmriprep_config_${sub}_${ses}.json \\
	--fs-no-reconall --fs-license-file /opts/fs_license.txt \\
	-w /scratch
EOF

    chmod +x "$scriptfile"

  done
done