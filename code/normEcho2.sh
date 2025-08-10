#!/usr/bin/env bash
set -euo pipefail

# --- Path bases (relative to where the script lives) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
fmriprepdir="${maindir}/derivatives/fmriprep"

# --- TemplateFlow res-2 reference (controls output grid/resolution) ---
TFLOW="/home/tug87422/work/tools/templateflow/tpl-MNI152NLin6Asym"
REF_MNI="${TFLOW}/tpl-MNI152NLin6Asym_res-02_T1w.nii.gz"
if [[ ! -f "$REF_MNI" ]]; then
  echo "ERROR: MNI reference not found at: $REF_MNI" >&2
  exit 1
fi
echo "Using MNI reference: $REF_MNI"

# --- Define entities to iterate (adjust as needed) ---
subs=(101 103)
#sessions=(01 02 03 04 05 06 07 08 09 10 11 12)
sessions=(01)
tasks=(mid sharedreward rest)
runs=(1 2)

# Optional: set FORCE=1 to overwrite existing outputs
: "${FORCE:=0}"

for sub in "${subs[@]}"; do
  for ses in "${sessions[@]}"; do
    for task in "${tasks[@]}"; do
      for run in "${runs[@]}"; do

        in4d="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_desc-preproc_bold.nii.gz"
        [[ -f "$in4d" ]] || { echo "Missing input, skip: $in4d"; continue; }

        # Transforms: boldref→T1w (per run) and MNI→T1w (subject- or session-level)
        boldref2t1="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"
        [[ -f "$boldref2t1" ]] || { echo "Missing boldref→T1w xfm, skip: $boldref2t1"; continue; }

        mni2t1="${fmriprepdir}/sub-${sub}/anat/sub-${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"
        if [[ ! -f "$mni2t1" ]]; then
          mni2t1="${fmriprepdir}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"
        fi
        [[ -f "$mni2t1" ]] || { echo "Missing MNI→T1w xfm, skip: $mni2t1"; continue; }

        # Output name: insert space/res before desc; keep echo/part
        out="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz"
        if [[ -f "$out" && "$FORCE" != "1" ]]; then
          echo "Exists, skipping: $out"
          continue
        fi

        echo "Warping: $in4d → $out"
        antsApplyTransforms \
          -d 3 -e 3 --float \
          -i "$in4d" \
          -r "$REF_MNI" \
          -o "$out" \
          -n Linear \
          -t "[$mni2t1,1]" \
          -t "$boldref2t1"

      done
    done
  done
done
