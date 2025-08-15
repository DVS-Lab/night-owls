#!/usr/bin/env bash
set -euo pipefail

# --- Path bases (relative to where the script lives) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
fmriprepdir="${maindir}/derivatives/fmriprep"

# --- Define entities (edit as needed) ---
subs=(101 103)
sessions=(01)                 # expand if needed
#sessions=(01 02 03 04 05 06 07 08 09 10 11 12)
tasks=(mid sharedreward rest)
runs=(1 2)                    # echo-2 only, per your original

# Optional: set FORCE=1 to overwrite existing outputs
: "${FORCE:=0}"

for sub in "${subs[@]}"; do
  for ses in "${sessions[@]}"; do
    for task in "${tasks[@]}"; do
      for run in "${runs[@]}"; do

        in4d="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_desc-preproc_bold.nii.gz"
        [[ -f "$in4d" ]] || { echo "Missing input, skip: $in4d"; continue; }

        # Per-run boldref→T1w transform (text affine from fMRIPrep)
        boldref2t1="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"
        [[ -f "$boldref2t1" ]] || { echo "Missing boldref→T1w xfm, skip: $boldref2t1"; continue; }

        # T1w↔MNI composite: prefer forward (T1w→MNI); otherwise invert the inverse (MNI→T1w)
        t1_to_mni=""
        try_forward="${fmriprepdir}/sub-${sub}/anat/sub-${sub}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5"
        try_forward_ses="${fmriprepdir}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5"
        try_inverse="${fmriprepdir}/sub-${sub}/anat/sub-${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"
        try_inverse_ses="${fmriprepdir}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"

        if   [[ -f "$try_forward"     ]]; then t1_to_mni="$try_forward"
        elif [[ -f "$try_forward_ses" ]]; then t1_to_mni="$try_forward_ses"
        elif [[ -f "$try_inverse"     ]]; then t1_to_mni="[$try_inverse,1]"
        elif [[ -f "$try_inverse_ses" ]]; then t1_to_mni="[$try_inverse_ses,1]"
        else
          echo "Missing T1w↔MNI xfm for sub-${sub} ses-${ses} — skipping."
          continue
        fi

        # 3-D references (use run-1 boldrefs to define output grids)
        ref_t1w="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-T1w_boldref.nii.gz"
        ref_mni="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-MNI152NLin6Asym_boldref.nii.gz"
        [[ -f "$ref_t1w" ]] || { echo "Missing T1w boldref: $ref_t1w"; continue; }
        [[ -f "$ref_mni" ]] || { echo "Missing MNI boldref: $ref_mni"; continue; }

        # --- Outputs (NO res-2 tags) ---
        out_t1w="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-T1w_desc-preproc_bold.nii.gz"
        out_mni="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"

        # Skip if exist (unless FORCE=1)
        if [[ -f "$out_t1w" && "$FORCE" != "1" ]]; then echo "Exists, skipping: $out_t1w"; else
          echo "Warping to T1w: $in4d → $out_t1w"
          antsApplyTransforms \
            -d 3 -e 3 --float \
            -i "$in4d" \
            -r "$ref_t1w" \
            -o "$out_t1w" \
            -n Linear \
            -t "$boldref2t1"
        fi

        if [[ -f "$out_mni" && "$FORCE" != "1" ]]; then echo "Exists, skipping: $out_mni"; else
          echo "Warping to MNI (regular): $in4d → $out_mni"
          antsApplyTransforms \
            -d 3 -e 3 --float \
            -i "$in4d" \
            -r "$ref_mni" \
            -o "$out_mni" \
            -n Linear \
            -t "$t1_to_mni" \
            -t "$boldref2t1"
        fi

      done
    done
  done
done
