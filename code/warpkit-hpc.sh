#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -N warpkit
#PBS -q normal
#PBS -l nodes=1:ppn=28

module load fsl/6.0.2
source $FSLDIR/etc/fslconf/fsl.sh
module load singularity
cd $PBS_O_WORKDIR


IFS=' ' read -r -a pairs <<< "$PAIRS"

maindir=/gpfs/scratch/tug87422/smithlab-shared/night-owls
toolsdir=/gpfs/scratch/tug87422/smithlab-shared/tools
scriptdir=$maindir/code
logdir=$maindir/logs
mkdir -p $logdir

umask 0000

rm -rf $scriptdir/warpkit.o*
rm -rf $scriptdir/warpkit.e*

TEMPLATEFLOW_DIR=$toolsdir/templateflow
MPLCONFIGDIR_DIR=$toolsdir/mplconfigdir
export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow
export SINGULARITYENV_MPLCONFIGDIR=/opt/mplconfigdir

echo "Processing pairs: ${pairs[@]}"

for pair in "${pairs[@]}"; do
	IFS=':' read -r sub ses <<< "$pair"
  	echo "Processing sub-${sub}, ses-${ses}"

	compiled_jobs="$logdir/cmd_warpkit_sub-${sub}_ses-${ses}_compiled.sh"
	rm -f $compiled_jobs

	echo "#!/bin/bash" > $compiled_jobs
	echo "" >> $compiled_jobs

	for task in "sharedreward" "mid" "rest"; do
		for run in 1 2; do

			outdir=$maindir/derivatives/warpkit/sub-${sub}/ses-${ses}
			if [ ! -d $outdir ]; then
				mkdir -p $outdir
			fi

			fmapdir=$maindir/bids/sub-${sub}/ses-${ses}/fmap
			if [ ! -d $fmapdir ]; then
				mkdir -p $fmapdir
			fi
	
			indir=${maindir}/bids/sub-${sub}/ses-${ses}/func
	
			if [ ! -e $indir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-mag_bold.json ]; then
				echo "NO DATA: sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-mag_bold.json"
				echo "NO DATA: sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-mag_bold.json" >> $scriptdir/missingFiles-warpkit.log
				continue
			fi
	
			# don't re-do existing output
			if [ -e $maindir/bids/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-${task}_run-${run}_fieldmap.nii.gz ]; then
				echo "EXISTS (skipping): sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-${task}_run-${run}_fieldmap.nii.gz"
				continue
			fi

			# delete default gre fmaps (phasediff)
			if [ -e $maindir/bids/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-bold_phasediff.nii.gz ]; then
				rm -rf $maindir/bids/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-bold*
			fi
	
			echo "singularity run --cleanenv \
-B $indir:/base \
-B $outdir:/out \
$toolsdir/warpkit.sif \
--magnitude /base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-mag_bold.nii.gz \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_bold.nii.gz \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-3_part-mag_bold.nii.gz \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-4_part-mag_bold.nii.gz \
--phase /base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-phase_bold.nii.gz \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-phase_bold.nii.gz \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-3_part-phase_bold.nii.gz \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-4_part-phase_bold.nii.gz \
--metadata /base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-phase_bold.json \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-phase_bold.json \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-3_part-phase_bold.json \
/base/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-4_part-phase_bold.json \
--out_prefix /out/sub-${sub}_ses-${ses}_task-${task}_run-${run}" >> $compiled_jobs
			echo "" >> $compiled_jobs

		done
	done
	
	chmod +x $compiled_jobs
	
	bash $logdir/cmd_warpkit_sub-${sub}_ses-${ses}_compiled.sh

	# Post-processing for this subject (after all jobs complete)
	for task in "sharedreward" "mid" "rest"; do
		for run in 1 2; do
			outdir=$maindir/derivatives/warpkit/sub-${sub}/ses-${ses}
			indir=${maindir}/bids/sub-${sub}/ses-${ses}/func
			
			# Check if the warpkit output exists before processing
			if [ -e $outdir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_fieldmaps.nii ]; then
				# extract first volume as fieldmap and copy to fmap dir
				fslroi $outdir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_fieldmaps.nii $maindir/bids/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-${task}_run-${run}_fieldmap 0 1
				fslroi $indir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-mag_bold.nii.gz $maindir/bids/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-${task}_run-${run}_magnitude 0 1

				# placeholders for json files
				cp $indir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-mag_bold.json $maindir/bids/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-${task}_run-${run}_magnitude.json
				cp $indir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-phase_bold.json $maindir/bids/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_acq-${task}_run-${run}_fieldmap.json

				# trash the rest
				rm -rf $outdir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_displacementmaps.nii
				rm -rf $outdir/sub-${sub}_ses-${ses}_task-${task}_run-${run}_fieldmaps_native.nii
			fi
    	done
	done
done

