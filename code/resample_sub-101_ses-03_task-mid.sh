#!/usr/bin/env bash
set -euo pipefail

# --- standard header: paths relative to the code directory ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
fmriprepdir="${maindir}/derivatives/fmriprep"

# --- edit these if you want to target another subject/session/task/run ---
sub=101
ses=03
task=mid
run=2        # regrid this run to match run-1's reference grid

# --- small helpers ----------------------------------------------------------

need () { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found in PATH"; exit 1; }; }
need antsApplyTransforms
need fslval
need fslstats

# echo and check a file
must_exist () {
  local f="$1" msg="${2:-Missing required file}"
  [[ -f "$f" ]] || { echo "ERROR: $msg → $f"; exit 1; }
}

# regrid $src to $ref with identity transform.
# - auto-detects 3D vs 4D via dim4
# - interpolation passed as $3 (default Linear); masks should use NearestNeighbor
regrid_identity () {
  local ref="$1" src="$2" label="$3" interp="${4:-Linear}"

  must_exist "$ref" "Reference missing"
  must_exist "$src" "Source missing"

  local d4
  d4="$(fslval "$src" dim4 2>/dev/null || echo 1)"
  local ants_args=(-d 3 --float -i "$src" -r "$ref" -o "${src%.nii.gz}_tmp.nii.gz" -t identity)
  if [[ "${d4}" -gt 1 ]]; then
    ants_args=(-d 3 -e 3 --float -i "$src" -r "$ref" -o "${src%.nii.gz}_tmp.nii.gz" -t identity)
  fi

  # one-time backup
  local backup="${src%.nii.gz}.pre-resample.nii.gz"
  if [[ ! -f "$backup" ]]; then
    cp -n "$src" "$backup"
    echo "[$label] Backed up original → $backup"
  fi

  echo "[$label] Resampling (identity, interp=${interp})"
  antsApplyTransforms "${ants_args[@]}" -n "${interp}"

  # quick sanity check
  local tmp="${src%.nii.gz}_tmp.nii.gz"
  read -r _min _max < <(fslstats "$tmp" -R)
  if [[ "${_max:-0}" == "0" || "${_max:-0}" == "0.000000" ]]; then
    echo "[$label] ERROR: output looks empty (max=0). Aborting."
    rm -f "$tmp"
    exit 1
  fi

  mv -f "$tmp" "$src"
  echo "[$label] OK → $src"
}

# --- paths ------------------------------------------------------------------

funcdir="${fmriprepdir}/sub-${sub}/ses-${ses}/func"

# T1w reference is run-1 boldref
ref_t1w="${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-T1w_boldref.nii.gz"
# run-X BOLD & MASK (to be regridded to ref_t1w)
bold_t1w="${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-T1w_desc-preproc_bold.nii.gz"
mask_t1w="${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-T1w_desc-brain_mask.nii.gz"

# MNI reference is run-1 boldref
ref_mni="${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-MNI152NLin6Asym_boldref.nii.gz"
# run-X BOLD & MASK (to be regridded to ref_mni)
bold_mni="${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
mask_mni="${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-brain_mask.nii.gz"

# --- do it ------------------------------------------------------------------

echo "Subject: ${sub}  Session: ${ses}  Task: ${task}  Run: ${run}"
echo "Regridding run-${run} to match run-1 in both spaces (T1w & MNI152NLin6Asym)."

# T1w 4D BOLD
regrid_identity "$ref_t1w" "$bold_t1w" "BOLD T1w run-${run}" "Linear"
# T1w 3D MASK (nearest neighbor)
regrid_identity "$ref_t1w" "$mask_t1w" "MASK T1w run-${run}" "NearestNeighbor"

# MNI 4D BOLD
regrid_identity "$ref_mni" "$bold_mni" "BOLD MNI run-${run}" "Linear"
# MNI 3D MASK (nearest neighbor)
regrid_identity "$ref_mni" "$mask_mni" "MASK MNI run-${run}" "NearestNeighbor"

echo "Done."
