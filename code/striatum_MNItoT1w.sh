#!/usr/bin/env bash
set -euo pipefail

# ----- directory anchors (as in your other scripts) -----
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# ----- inputs -----
atlas_space="MNI152NLin6Asym"
atlas_mask="${maindir}/masks/StriatumMask_atlas.nii"          # primary
[[ -f "${atlas_mask}" ]] || atlas_mask="${maindir}/masks/StriatumMask_atlas.nii.gz"

if [[ ! -f "${atlas_mask}" ]]; then
  echo "[ERROR] Could not find masks/StriatumMask_atlas.nii(.gz) under ${maindir}" >&2
  exit 1
fi

fmriprep_dir="${maindir}/derivatives/fmriprep"
outdir="${maindir}/masks"
mkdir -p "${outdir}"

# ----- tool checks -----
command -v antsApplyTransforms >/dev/null 2>&1 || { echo "[ERROR] antsApplyTransforms not found in PATH." >&2; exit 1; }

# ----- enumerate subjects from fMRIPrep -----
mapfile -t subs < <(find "${fmriprep_dir}" -maxdepth 1 -type d -name 'sub-*' -printf '%f\n' | sort)
if (( ${#subs[@]} == 0 )); then
  echo "[ERROR] No sub-* folders found under ${fmriprep_dir}" >&2
  exit 1
fi

echo "[INFO] Found ${#subs[@]} subjects. Creating T1w-space striatum masks…"

# progress
done_cnt=0
total=${#subs[@]}

for sub in "${subs[@]}"; do
  sid="${sub#sub-}"

  # --- choose subject-level T1w reference (average across sessions, if present) ---
  t1w_ref="${fmriprep_dir}/${sub}/anat/${sub}_desc-preproc_T1w.nii.gz"
  if [[ ! -f "${t1w_ref}" ]]; then
    # fallback: first session T1w (only if subject-level T1w not present)
    t1w_ref=$(find "${fmriprep_dir}/${sub}" -type f -path "*/ses-*/anat/${sub}_ses-*_desc-preproc_T1w.nii.gz" | sort | head -n1 || true)
  fi
  if [[ -z "${t1w_ref}" || ! -f "${t1w_ref}" ]]; then
    echo "[WARN] ${sub}: no preproc T1w found; skipping."
    continue
  fi

  # --- prefer subject-level MNI->T1w transform (no ses) ---
  xfm="${fmriprep_dir}/${sub}/anat/${sub}_from-${atlas_space}_to-T1w_mode-image_xfm.h5"
  if [[ ! -f "${xfm}" ]]; then
    # try session-specific transforms if subject-level H5 not present
    xfm=$(find "${fmriprep_dir}/${sub}" -type f -path "*/ses-*/anat/${sub}_ses-*_from-${atlas_space}_to-T1w_mode-image_xfm.h5" | sort | head -n1 || true)
  fi
  if [[ -z "${xfm}" || ! -f "${xfm}" ]]; then
    echo "[WARN] ${sub}: no ${atlas_space}->T1w transform found; skipping."
    continue
  fi

  outmask="${outdir}/sub-${sid}_space-T1w_desc-StriatumMask_atlas_mask.nii.gz"
  if [[ -f "${outmask}" ]]; then
    echo "[INFO] ${sub}: output exists, skipping (${outmask})"
  else
    echo "[INFO] ${sub}: applying ${atlas_space}->T1w transform"
    antsApplyTransforms \
      -d 3 \
      -i "${atlas_mask}" \
      -r "${t1w_ref}" \
      -t "${xfm}" \
      -n NearestNeighbor \
      -o "${outmask}"

    # quick sanity check: nonzero voxel count if FSL is available
    if command -v fslstats >/dev/null 2>&1; then
      nz=$(fslstats "${outmask}" -V | awk '{print $1}')
      echo "[INFO] ${sub}: wrote ${outmask} (nonzero voxels: ${nz})"
    else
      echo "[INFO] ${sub}: wrote ${outmask}"
    fi
  fi

  ((done_cnt++))
  printf '[%s] %d/%d (%.1f%%) subjects done\n' \
    "$(date '+%F %T')" "${done_cnt}" "${total}" "$(awk "BEGIN{print 100*${done_cnt}/${total}}")"
done

echo "[INFO] Completed. Outputs in: ${outdir}"
