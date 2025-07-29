#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N tedana
#PBS -q normal
#PBS -l nodes=1:ppn=28

# load modules and go to workdir
module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
module load singularity/3.8.5
cd $PBS_O_WORKDIR

# ensure paths are correct
projectname=nightowls #this should be the only line that has to change if the rest of the script is set up correctly
maindir=/gpfs/scratch/tug87422/smithlab-shared/$projectname
scriptdir=$maindir/code
bidsdir=$maindir/bids
logdir=$maindir/logs

rm -f $logdir/cmd_tedana_${PBS_JOBID}.txt
touch $logdir/cmd_tedana_${PBS_JOBID}.txt

# estimated procs per subject: 4 tasks * 2 runs * 4 echoes = 32?
# wow... can't be right, but let's go for 9 subjects per job (84 procs per job). watch for memory issues.
for sub in ${subjects[@]}; do
	for ses in ${01..12};
		for task in mid sharedreward rest; do
			for run in 1 2; do

				indata=$maindir/derivatives/fmriprep/sub-${sub}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_desc-preproc_bold.nii.gz
				if [ -e $indata ]; then
					echo python $scriptdir/tedana-single.py \
					--fmriprepDir $maindir/derivatives/fmriprep \
					--bidsDir $bidsdir \
					--sub $sub \
					--task $task \
					--runnum $run >> $logdir/cmd_tedana_${PBS_JOBID}.txt
				fi

			done
		done
	done
done

torque-launch -p $logdir/chk_tedana_${PBS_JOBID}.txt $logdir/cmd_tedana_${PBS_JOBID}.txt
