#!/bin/bash




# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#maindir="$(dirname "$scriptdir")"

# Define the data directory
datadir=/ZPOOL/data/projects/night-owls/derivatives/fmriprep
tasks=("mid" "sharedreward")

#echo -e "sub\tmean\tmax\tmeanVS_run1\tmeanVS_run2\t"
echo -e "sub\t mean_stan_run1\t mean_stan_run2\t mean_nat_run1\t mean_nat_run2\t max\t vsmean_run1\t vsmean_run2" 



# Process _stan data
for i in sub-${sub}_task-${task}_run-*_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz; do
    fslmaths "$i" -Tmean tmp_mean
    fslmaths "$i" -Tstd tmp_std
    fslmaths tmp_mean -div tmp_std tmp_tsnr
    fslmaths tmp_tsnr -thr 2 thr_tmp_tsnr
    max=$(fslstats thr_tmp_tsnr -R | awk '{ print $2 }')
    mean_stan_run1=$(fslstats thr_tmp_tsnr -k /ZPOOL/data/projects/rf1-sra-data/derivatives/fmriprep/sub-${sub}/func/sub-${sub}_task-${task}_run-1_space-MNI152NLin6Asym_desc-brain_mask.nii.gz -M)
    mean_stan_run2=$(fslstats thr_tmp_tsnr -k /ZPOOL/data/projects/rf1-sra-data/derivatives/fmriprep/sub-${sub}/func/sub-${sub}_task-${task}_run-2_space-MNI152NLin6Asym_desc-brain_mask.nii.gz -M)
    echo -e "$i\t $mean_stan_run1\t $mean_stan_run2\t -\t -\t $max\t -\t -"
done

# Process _nat data
for i in sub-${sub}_task-${task}_run-*_echo-*_desc-preproc_bold.nii.gz; do 
    fslmaths "$i" -Tmean tmp_mean
    fslmaths "$i" -Tstd tmp_std
    fslmaths tmp_mean -div tmp_std tmp_tsnr
    fslmaths tmp_tsnr -thr 2 thr_tmp_tsnr
    max=$(fslstats thr_tmp_tsnr -R | awk '{ print $2 }')
    mean_nat_run1=$(fslstats thr_tmp_tsnr -k /ZPOOL/data/projects/rf1-sra-data/derivatives/fmriprep/sub-${sub}/func/sub-${sub}_task-${task}_run-1_desc-brain_mask.nii.gz -M)
    mean_nat_run2=$(fslstats thr_tmp_tsnr -k /ZPOOL/data/projects/rf1-sra-data/derivatives/fmriprep/sub-${sub}/func/sub-${sub}_task-${task}_run-2_desc-brain_mask.nii.gz -M)
    vsmean_run1=$(fslstats thr_tmp_tsnr -k /ZPOOL/data/projects/rf1-sra-data/derivatives/fmriprep/sub-${sub}/sub-${sub}_task-${task}_run-1_space-native_roi-vs_thr_mask.nii.gz -M)
    vsmean_run2=$(fslstats thr_tmp_tsnr -k /ZPOOL/data/projects/rf1-sra-data/derivatives/fmriprep/sub-${sub}/sub-${sub}_task-${task}_run-2_space-native_roi-vs_thr_mask.nii.gz -M)
    echo -e "$i\t -\t -\t $mean_nat_run1\t $mean_nat_run2\t $max\t $vsmean_run1\t $vsmean_run2"
done
