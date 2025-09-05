#!/usr/bin/env bash
# Transform the group Striatum mask (MNI152NLin6Asym) into each subject's T1w space
# and write outputs into night-owls/masks as sub-XXX_space-T1w_desc-StriatumMask_atlas.nii.gz

set -euo pipefail

# --- canonical project roots: all scripts are run from code/ ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# --- inputs/outputs ---
fmriprep_dir="${maindir}/derivatives/fmriprep"
masks_dir="${maindir}/masks"
atlas_space="MNI152NLin6Asym"

# accept either .nii or .nii.gz for the source atlas mask
atlas_mask="${masks_dir}/StriatumMask_atlas.nii"
[[ -f "${atlas_mask}" ]] || atlas_mask="${masks_dir}/StriatumMask_atlas.nii.gz"

if [[ ! -f "${atlas_mask}" ]]; then
  echo "[FATAL] Missing atlas mask: ${masks_dir}/StriatumMask_atlas.nii[.gz]"
  exit 1
fi

# --- subject discovery ---
mapfile -t subs < <(find "${fmriprep_dir}" -maxdepth 1 -type d -name "sub-*" -printf "%f\n" | sort)

echo "[INFO] Found ${#subs[@]} subjects. Creating T1w-space striatum masks…"

total="${#subs[@]}"
done_count=0

for sub in "${subs[@]}"; do
  sid="${sub#sub-}"

  # --- choose subject-level T1w reference (or first session T1w) ---
  t1w_ref="${fmriprep_dir}/${sub}/anat/${sub}_desc-preproc_T1w.nii.gz"
  if [[ ! -f "${t1w_ref}" ]]; then
    # session folder, filename WITHOUT ses tag (your dataset case)
    t1w_ref=$(find "${fmriprep_dir}/${sub}" -type f -path "*/ses-*/anat/${sub}_desc-preproc_T1w.nii.gz" | sort | head -n1 || true)
  fi
  if [[ -z "${t1w_ref}" || ! -f "${t1w_ref}" ]]; then
    # session folder, filename WITH ses tag (fallback)
    t1w_ref=$(find "${fmriprep_dir}/${sub}" -type f -path "*/ses-*/anat/${sub}_ses-*_desc-preproc_T1w.nii.gz" | sort | head -n1 || true)
  fi
  if [[ -z "${t1w_ref}" || ! -f "${t1w_ref}" ]]; then
    echo "[WARN] ${sub}: no preproc T1w found; skipping."
    continue
  fi

  # --- prefer subject-level MNI->T1w transform, else session-level ---
  xfm="${fmriprep_dir}/${sub}/anat/${sub}_from-${atlas_space}_to-T1w_mode-image_xfm.h5"
  if [[ ! -f "${xfm}" ]]; then
    # session folder, filename WITHOUT ses tag (your dataset case)
    xfm=$(find "${fmriprep_dir}/${sub}" -type f -path "*/ses-*/anat/${sub}_from-${atlas_space}_to-T1w_mode-image_xfm.h5" | sort | head -n1 || true)
  fi
  if [[ -z "${xfm}" || ! -f "${xfm}" ]]; then
    # session folder, filename WITH ses tag (fallback)
    xfm=$(find "${fmriprep_dir}/${sub}" -type f -path "*/ses-*/anat/${sub}_ses-*_from-${atlas_space}_to-T1w_mode-image_xfm.h5" | sort | head -n1 || true)
  fi
  if [[ -z "${xfm}" || ! -f "${xfm}" ]]; then
    echo "[WARN] ${sub}: no ${atlas_space}->T1w transform found; skipping."
    continue
  fi

  out_mask="${masks_dir}/sub-${sid}_space-T1w_desc-StriatumMask_atlas.nii.gz"
  if [[ -f "${out_mask}" ]]; then
    echo "[INFO] ${sub}: output exists, skipping (${out_mask})."
  else
    # --- apply transform using ANTs (no FSL flirt; masks use nearest-neighbor) ---
    antsApplyTransforms \
      -d 3 \
      -i "${atlas_mask}" \
      -r "${t1w_ref}" \
      -o "${out_mask}" \
      -t "${xfm}" \
      -n NearestNeighbor

    echo "[OK]   ${sub}: wrote ${out_mask}"
  fi

  # progress
  ((done_count++))
  pct=$(( 100 * done_count / total ))
  echo "[PROGRESS] ${pct}% (${done_count}/${total})"
done

echo "[INFO] Completed. Outputs in: ${masks_dir}"
