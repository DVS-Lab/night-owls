#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Striatum atlas (MNI152NLin6Asym) -> subject native T1w space
# Writes one transformed mask per subject into ./masks
# ------------------------------------------------------------

# Standard project roots (all scripts run from ./code)
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

fmriprep_dir="${maindir}/derivatives/fmriprep"
outdir="${maindir}/masks"

# Atlas in MNI152NLin6Asym space (accept .nii or .nii.gz)
atlas="${maindir}/masks/StriatumMask_atlas.nii"
[[ -f "${atlas}" ]] || atlas="${maindir}/masks/StriatumMask_atlas.nii.gz"

if [[ ! -f "${atlas}" ]]; then
  echo "[ERROR] Could not find StriatumMask_atlas (.nii or .nii.gz) in ${maindir}/masks"
  exit 1
fi

# Require ANTs
if ! command -v antsApplyTransforms >/dev/null 2>&1 ; then
  echo "[ERROR] antsApplyTransforms not found in PATH."
  exit 1
fi

mkdir -p "${outdir}"

# Helper: pick native T1w (exclude any space-MNI* variants). If not found,
# fall back to the native-space brain mask as the reference grid.
pick_t1w_ref() {
  local subj_dir="$1" sub="$2"
  local ref=""

  # Subject-level native T1w
  ref="$(find "${subj_dir}/anat" -maxdepth 1 -type f -name "${sub}_desc-preproc_T1w.nii.gz" | sort | head -n1)"

  # Session-level native T1w (no space- token)
  [[ -z "$ref" ]] && ref="$(find "${subj_dir}"/ses-*/anat -type f -name "${sub}_desc-preproc_T1w.nii.gz" | sort | head -n1)"
  [[ -z "$ref" ]] && ref="$(find "${subj_dir}"/ses-*/anat -type f -name "${sub}_ses-*_desc-preproc_T1w.nii.gz" | sort | head -n1)"

  # Fallback: native-space brain mask as reference grid
  [[ -z "$ref" ]] && ref="$(find "${subj_dir}/anat" -maxdepth 1 -type f -name "${sub}_desc-brain_mask.nii.gz" | sort | head -n1)"
  [[ -z "$ref" ]] && ref="$(find "${subj_dir}"/ses-*/anat -type f -name "${sub}_desc-brain_mask.nii.gz" | sort | head -n1)"
  [[ -z "$ref" ]] && ref="$(find "${subj_dir}"/ses-*/anat -type f -name "${sub}_ses-*_desc-brain_mask.nii.gz" | sort | head -n1)"

  echo "$ref"
}

# Collect subjects present in fmriprep
mapfile -t subjects < <(find "${fmriprep_dir}" -maxdepth 1 -type d -name "sub-*" -printf "%f\n" | sort)

echo "[INFO] Found ${#subjects[@]} subjects. Creating T1w-space striatum masks…"

idx=0
for sub in "${subjects[@]}"; do
  idx=$((idx+1))
  subj_dir="${fmriprep_dir}/${sub}"

  # Reference in native T1w space
  t1w_ref="$(pick_t1w_ref "${subj_dir}" "${sub}")"
  if [[ -z "${t1w_ref}" || ! -f "${t1w_ref}" ]]; then
    echo "[WARN] ${sub}: no native T1w (or T1w brain mask) found; skipping."
    continue
  fi
  echo "       ${sub}: REF -> ${t1w_ref}"

  # Transform: MNI152NLin6Asym -> T1w (prefer subject-level, else session-level)
  xfm="$(find "${subj_dir}/anat" -maxdepth 1 -type f -name "${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" | sort | head -n1)"
  [[ -z "${xfm}" ]] && xfm="$(find "${subj_dir}"/ses-*/anat -type f -name "${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" | sort | head -n1)"
  [[ -z "${xfm}" ]] && xfm="$(find "${subj_dir}"/ses-*/anat -type f -name "${sub}_ses-*_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" | sort | head -n1)"

  if [[ -z "${xfm}" || ! -f "${xfm}" ]]; then
    echo "[WARN] ${sub}: no MNI152NLin6Asym→T1w transform found; skipping."
    continue
  fi
  echo "       ${sub}: XFM -> ${xfm}"

  out="${outdir}/${sub}_space-T1w_desc-StriatumMask_atlas.nii.gz"

  # Apply transform (label mask -> use NN interpolation)
  antsApplyTransforms \
    -d 3 \
    -i "${atlas}" \
    -r "${t1w_ref}" \
    -t "${xfm}" \
    -n NearestNeighbor \
    -o "${out}"

  if [[ -f "${out}" ]]; then
    echo "[OK]   ${sub}: wrote ${out}"
  else
    echo "[ERR]  ${sub}: failed to write ${out}"
  fi
done

echo "[INFO] Completed. Outputs in: ${outdir}"
