#!/usr/bin/env bash
set -euo pipefail

# ---- paths relative to this script ----
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
fmriprepdir="${maindir}/derivatives/fmriprep"

# ---- edit as needed ----
subs=(103)              # e.g., (101 103)
sessions=(01)           # expand as needed
tasks=(mid sharedreward rest)
runs=(1 2)              # echo-2 only
: "${FORCE:=0}"

pick_first_existing() {
  # usage: pick_first_existing VAR path1 path2 ...
  local __outvar="$1"; shift
  local p
  for p in "$@"; do
    if [[ -f "$p" ]]; then printf -v "$__outvar" '%s' "$p"; return 0; fi
  done
  return 1
}

for sub in "${subs[@]}"; do
  for ses in "${sessions[@]}"; do
    for task in "${tasks[@]}"; do
      # 3-D reference grids (run-1)
      ref_t1w="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-T1w_boldref.nii.gz"
      ref_mni="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-MNI152NLin6Asym_boldref.nii.gz"
      [[ -f "$ref_t1w" ]] || { echo "Missing T1w boldref: $ref_t1w"; continue; }
      [[ -f "$ref_mni" ]] || { echo "Missing MNI boldref: $ref_mni"; continue; }

      # T1w↔MNI composite: check all plausible locations/names
      if ! pick_first_existing t1_to_mni \
        "${fmriprepdir}/sub-${sub}/anat/sub-${sub}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5" \
        "${fmriprepdir}/sub-${sub}/anat/sub-${sub}_ses-${ses}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5" \
        "${fmriprepdir}/sub-${sub}/ses-${ses}/anat/sub-${sub}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5" \
        "${fmriprepdir}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5" \
      ; then
        # fall back to inverse (invert on the fly)
        if ! pick_first_existing mni_to_t1 \
          "${fmriprepdir}/sub-${sub}/anat/sub-${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
          "${fmriprepdir}/sub-${sub}/anat/sub-${sub}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
          "${fmriprepdir}/sub-${sub}/ses-${ses}/anat/sub-${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
          "${fmriprepdir}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
        ; then
          echo "WARN sub-${sub} ses-${ses}: no T1w↔MNI .h5 found; will write T1w only."
          t1_to_mni=""
        else
          t1_to_mni="[${mni_to_t1},1]"
        fi
      fi

      for run in "${runs[@]}"; do
        in4d="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_desc-preproc_bold.nii.gz"
        [[ -f "$in4d" ]] || { echo "Missing input: $in4d"; continue; }

        # per-run boldref→T1w (exists in func/)
        boldref2t1="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"
        [[ -f "$boldref2t1" ]] || { echo "Missing boldref→T1w xfm: $boldref2t1"; continue; }

        # ---- space-T1w ----
        out_t1w="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-T1w_desc-preproc_bold.nii.gz"
        if [[ -f "$out_t1w" && "$FORCE" != "1" ]]; then
          echo "Exists, skipping: $out_t1w"
        else
          echo "T1w: $in4d → $out_t1w"
          antsApplyTransforms -d 3 -e 3 --float \
            -i "$in4d" \
            -r "$ref_t1w" \
            -o "$out_t1w" \
            -n Linear \
            -t "$boldref2t1"
        fi

        # ---- space-MNI152NLin6Asym (no res-2) ----
        if [[ -n "$t1_to_mni" ]]; then
          out_mni="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
          if [[ -f "$out_mni" && "$FORCE" != "1" ]]; then
            echo "Exists, skipping: $out_mni"
          else
            echo "MNI: $in4d → $out_mni"
            # order: last applied first → boldref→T1w then T1w→MNI
            antsApplyTransforms -d 3 -e 3 --float \
              -i "$in4d" \
              -r "$ref_mni" \
              -o "$out_mni" \
              -n Linear \
              -t "$t1_to_mni" \
              -t "$boldref2t1"
          fi
        fi

      done
    done
  done
done
