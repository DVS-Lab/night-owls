#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N tedana
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=28
cd $PBS_O_WORKDIR

# ensure paths are correct irrespective from where user runs the script
projectname=night-owls
maindir=/gpfs/scratch/tug87422/smithlab-shared/$projectname
scriptdir=$maindir/code
logdir=$maindir/logs

mkdir -p $logdir

rm -f $logdir/cmd_tedana_${PBS_JOBID}.txt
touch $logdir/cmd_tedana_${PBS_JOBID}.txt

for sub in ${subjects[@]}; do
	for ses in {01..12}; do

        prepdir=$maindir/derivatives/anat-only/ses-${ses}/sub-${sub}/ses-${ses}/func
        [[ ! -d "$prepdir" ]] && continue

        for task in mid sharedreward rest; do
            for run in {1..2}; do

                # prepare inputs and outputs
                echo1=${prepdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-1_part-mag_desc-preproc_bold.nii.gz
                echo2=${prepdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_desc-preproc_bold.nii.gz
                echo3=${prepdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-3_part-mag_desc-preproc_bold.nii.gz
                echo4=${prepdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-4_part-mag_desc-preproc_bold.nii.gz
                
                outdir=${maindir}/derivatives/anat-only/tedana/sub-${sub}/ses-${ses}

                # Check for the presence of all echo files
                if [ ! -e $echo1 ] || [ ! -e $echo2 ] || [ ! -e $echo3 ] || [ ! -e $echo4 ]; then
                    echo "Missing one or more files for sub-${sub}, ses-${ses}, task-${task}, run-${run}" >> $scriptdir/missing-tedanaInput.log
                    echo "Skipping sub-${sub}, ses-${ses}, task-${task}, run-${run}" >> $logdir/cmd_tedana_${PBS_JOBID}.txt
                    continue
                fi

                mkdir -p $outdir
                
                    # Initialize echo time variables
                echotime1=""
                echotime2=""
                echotime3=""
                echotime4=""

                # Extract echo times from the first script output
                for echo in 1 2 3 4; do
                    json_file=$(find "$maindir/bids" -name "sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-${echo}_part-mag_bold.json")
                    if [ -n "$json_file" ]; then
                        echo_time=$(grep -o '"EchoTime": [0-9.]*' "$json_file" | cut -d' ' -f2 | tr -d '\r')
                        eval "echotime${echo}=${echo_time}"
                    else
                        echo "missing JSON for echo-${echo} for sub-${sub}, ses-${ses}, task-${task}, run-${run}"
                        echo "missing JSON for echo-${echo} for sub-${sub}, ses-${ses}, task-${task}, run-${run}" >> $scriptdir/missing-tedanaInput.log
                    fi
                done

                echo "tedana -d $echo1 $echo2 $echo3 $echo4 \
                -e $echotime1 $echotime2 $echotime3 $echotime4 \
                --out-dir $outdir \
                --prefix sub-${sub}_ses-${ses}_task-${task}_run-${run} \
                --convention bids \
                --fittype curvefit \
                --overwrite"  >> $logdir/cmd_tedana_${PBS_JOBID}.txt

                            # clean up and save space
                   # rm -rf ${outdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_*.nii.gz

            done
		done
	done
done

torque-launch -p $logdir/chk_tedana_${PBS_JOBID}.txt $logdir/cmd_tedana_${PBS_JOBID}.txt
