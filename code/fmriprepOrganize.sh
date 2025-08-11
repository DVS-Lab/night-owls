#!/usr/bin/env bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 <subject_id>   e.g. $0 101"
  exit 1
fi

SUBJ="$1"

for ses in $(seq -w 1 12); do
  SRC="/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/anat-only/"
  DST="/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/fmriprep/sub-${SUBJ}/ses-${ses}"

  # Skip entire session if the main source folder is missing
  if [ ! -d "${SRC}/ses-${ses}/sub-${SUBJ}" ]; then
    echo "${ses} not found, skipping ses-${ses}" >&2
    continue
  fi
  mkdir -p "${DST}"

  #Move subject-level anat dir
  rsync -av "${SRC}/sub-${SUBJ}/anat/" "${DST}/"
  #Move ses-specific xfm and dseg files 
  rsync -av "${SRC}/sub-${SUBJ}/ses-${ses}/anat/sub-${SUBJ}_ses-${ses}_from-orig_to-T1w_mode-image_xfm.txt" "${DST}/anat/"
  rsync -av "${SRC}/ses-${ses}/sub-${SUBJ}/anat/sub-${SUBJ}_space-MNI152NLin6Asym_desc-preproc_dseg.nii.gz" "${DST}/anat/"

  #Move func and fmap 
  rsync -av "${SRC}/ses-${ses}/sub-${SUBJ}/ses-${ses}/func/" "${DST}/func"
  rsync -av "${SRC}/ses-${ses}/sub-${SUBJ}/ses-${ses}/fmap/" "${DST}/fmap"


  #Move figures
  rsync -av "${SRC}/sub-${SUBJ}_anat.html" "${DST}/"
  rsync -av "${SRC}/sub-${SUBJ}/figures/" "${DST}/"
  rsync -av "${SRC}/ses-${ses}/sub-${SUBJ}_ses-${ses}_func.html" "${DST}/"
  rsync -av "${SRC}/ses-${ses}/sub-${SUBJ}/figures" "${DST}/figures/"
done

#rm -rf /gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/anat-only/
#rm -rf /gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/ses-*
