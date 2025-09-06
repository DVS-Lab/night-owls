#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Striatum atlas (MNI152NLin6Asym) -> subject native T1w space
# Reference grid = T1w-space BOLDREF (matches FEAT/EPI resolution)
# One transformed mask per subject to ./masks
# ------------------------------------------------------------

# All scripts run from ./code
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

fmriprep_dir="${maindir}/derivatives/fmriprep"
outdir="${maindir}/masks"

# Atlas (in MNI152NLin6Asym)
atlas="${maindir}/masks/StriatumMask_atlas.nii"
[[ -f "$atlas" ]] || atlas="${maindir}/masks/StriatumMask_atlas.nii.gz"
if [[ ! -f "$atlas" ]]; then
  echo "[ERROR] Could not find StriatumMask_atlas (.nii or .nii.gz) in ${maindir}/masks"
  exit 1
fi

# Need ANTs
if ! command -v antsApplyTransforms >/dev/null 2>&1 ; then
  echo "[ERROR] antsApplyTransforms not found in PATH."
  exit 1
fi

mkdir -p "$outdir"

# ---- helpers ------------------------------------------------

# Prefer a T1w-space BOLDREF (EPI in T1w space) so voxel sizes match FEAT outputs
pick_t1w_ref() {  # $1=subj_dir  $2=sub-XXX
  local subj_dir="$1" sub="$2" ref=""

  # 1) Any func boldref in T1w space (covers part-mag or not)
  ref="$(find "${subj_dir}"/ses-*/func -type f \
          -name "${sub}_*_space-T1w_*boldref.nii.gz" 2>/dev/null | head -n1 || true)"

  # 2) Fallback: native (non-space-*) T1w image
  if [[ -z "$ref" && -d "${subj_dir}/anat" ]]; then
    ref="$(find "${subj_dir}/anat" -maxdepth 1 -type f \
            -name "${sub}_desc-preproc_T1w.nii.gz" 2>/dev/null | head -n1 || true)"
  fi
  if [[ -z "$ref" ]]; then
    ref="$(find "${subj_dir}"/ses-*/anat -type f \
            -name "${sub}_desc-preproc_T1w.nii.gz" 2>/dev/null | grep -v '/space-' | head -n1 || true)"
  fi

  # 3) Last resort: native brain mask grid
  if [[ -z "$ref" && -d "${subj_dir}/anat" ]]; then
    ref="$(find "${subj_dir}/anat" -maxdepth 1 -type f \
            -name "${sub}_desc-brain_mask.nii.gz" 2>/dev/null | head -n1 || true)"
  fi
  if [[ -z "$ref" ]]; then
    ref="$(find "${subj_dir}"/ses-*/anat -type f \
            -name "${sub}_desc-brain_mask.nii.gz" 2>/dev/null | head -n1 || true)"
  fi

  echo "$ref"
}

pick_mni_to_t1w_xfm() {  # $1=subj_dir  $2=sub-XXX
  local subj_dir="$1" sub="$2" xfm=""
  if [[ -d "${subj_dir}/anat" ]]; then
    xfm="$(find "${subj_dir}/anat" -maxdepth 1 -type f \
            -name "${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" 2>/dev/null | head -n1 || true)"
  fi
  if [[ -z "$xfm" ]]; then
    xfm="$(find "${subj_dir}"/ses-*/anat -type f \
            -name "${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" 2>/dev/null | head -n1 || true)"
  fi
  echo "$xfm"
}

# ---- subjects ------------------------------------------------

mapfile -t subjects < <(find "$fmriprep_dir" -maxdepth 1 -type d -name "sub-*" -printf "%f\n" | sort)
echo "[INFO] Found ${#subjects[@]} subjects. Creating T1w-space striatum masks…"

for sub in "${subjects[@]}"; do
  subj_dir="${fmriprep_dir}/${sub}"

  t1w_ref="$(pick_t1w_ref "$subj_dir" "$sub")"
  if [[ -z "$t1w_ref" || ! -f "$t1w_ref" ]]; then
    echo "[WARN] ${sub}: no suitable T1w-space reference (boldref/T1w) found; skipping."
    continue
  fi
  echo "       ${sub}: REF -> ${t1w_ref}"

  xfm="$(pick_mni_to_t1w_xfm "$subj_dir" "$sub")"
  if [[ -z "$xfm" || ! -f "$xfm" ]]; then
    echo "[WARN] ${sub}: no MNI152NLin6Asym→T1w transform found; skipping."
    continue
  fi
  echo "       ${sub}: XFM -> ${xfm}"

  out="${outdir}/${sub}_space-T1w_desc-StriatumMask_atlas.nii.gz"

  antsApplyTransforms \
    -d 3 \
    -i "$atlas" \
    -r "$t1w_ref" \
    -t "$xfm" \
    -n NearestNeighbor \
    -o "$out"

  if [[ -f "$out" ]]; then
    echo "[OK]   ${sub}: wrote ${out}"
  else
    echo "[ERR]  ${sub}: failed to write ${out}"
  fi
done

echo "[INFO] Completed. Outputs in: ${outdir}"
