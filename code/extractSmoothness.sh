#!/usr/bin/env bash
# Simplified extractor: compute AFNI smoothness (ACF) from smoothed and unsmoothed fMRIPrep preprocessed BOLD
# For each FEAT directory under derivatives/fsl, we:
#   1) locate the corresponding fMRIPrep preproc BOLD (unsmoothed) and the smoothed variant (e.g., *_bold_5mm.nii.gz)
#   2) run 3dFWHMx -detrend -ACF with the FEAT mask as -mask
#   3) write outputs into <featdir>/smoothness-<mm>mm.txt (smoothed) and <featdir>/smoothness-0mm.txt (unsmoothed)
# Notes:
#   - No zstat extraction is performed anymore.
#   - We do not assume a particular MB/ME label. We try both 'part-mag' and non-'part-mag' file patterns and pick what exists.
#   - We prefer sessioned fMRIPrep paths if ses-<id> is discoverable from FEAT path; otherwise we fall back to single-session layout.
#   - Requires AFNI (3dFWHMx) and FSL (for masks produced by FEAT).

set -euo pipefail

# -------- locations (relative to THIS script) --------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"
FSL_DERIV="${rootdir}/derivatives/fsl"
FMRIPREP_DERIV="${rootdir}/derivatives/fmriprep"

# -------- tools --------
command -v 3dFWHMx >/dev/null || { echo "ERROR: 3dFWHMx (AFNI) not found in PATH"; exit 1; }

OUT_DIR="${rootdir}/derivatives/extractions"
mkdir -p "${OUT_DIR}"
TSV="${OUT_DIR}/smoothness_acf.tsv"
if [[ ! -f "${TSV}" ]]; then
  echo -e "sub\tses\ttask\trun\tkernel_mm\tacf_a\tacf_b\tacf_c\tfwhm_x\tfwhm_y\tfwhm_z\tfeatdir\timg" > "${TSV}"
fi

# parse helpers
append_to_tsv() {
  local sub="$1" ses="$2" task="$3" run="$4" kernel="$5" featdir="$6" img="$7" txt="$8"
  # Default NA values
  local a="NA" b="NA" c="NA" fx="NA" fy="NA" fz="NA"

  # ACF parameters: try lines containing 'ACF', take last three numeric tokens
  if grep -qi "ACF" "$txt"; then
    read a b c < <(grep -i "ACF" "$txt" | tail -n1 | grep -Eo "[-+]?[0-9]*\.?[0-9]+" | tail -n3)
  fi

  # FWHM: try lines containing 'FWHM', take last three numeric tokens
  if grep -qi "FWHM" "$txt"; then
    read fx fy fz < <(grep -i "FWHM" "$txt" | tail -n1 | grep -Eo "[-+]?[0-9]*\.?[0-9]+" | tail -n3)
  fi

  echo -e "${sub}\t${ses}\t${task}\t${run}\t${kernel}\t${a}\t${b}\t${c}\t${fx}\t${fy}\t${fz}\t${featdir}\t${img}" >> "${TSV}"
}


# -------- main --------
# -------- main --------
echo "[$(date '+%F %T')] Starting smoothness extraction…"
shopt -s nullglob

# Iterate over FEAT directories (Level-1 results); adjust the find path to match your layout
while IFS= read -r -d '' featdir; do
  # Skip if mask is missing
  mask="${featdir}/mask.nii.gz"
  if [[ ! -f "$mask" ]]; then
    echo "WARN: No mask at ${mask}; skipping ${featdir}"
    continue
  fi

  # parse metadata from path
  IFS=$'\t' read -r sub ses task run <<<"$(parse_meta_from_path "$featdir")"
  if [[ -z "$sub" || -z "$task" ]]; then
    echo "WARN: Could not parse sub/task from ${featdir}; skipping"
    continue
  fi

  # find fMRIPrep raw+smoothed files
  IFS=$'\t' read -r raw_img smooth_img <<<"$(find_fmriprep_preproc "$sub" "$ses" "$task")"
  if [[ -z "$raw_img" ]]; then
    echo "WARN: No raw preproc image for sub-${sub} ses-${ses} task-${task}; skipping ${featdir}"
    continue
  fi
  if [[ -z "$smooth_img" ]]; then
    echo "WARN: No smoothed preproc image for sub-${sub} ses-${ses} task-${task}; will compute unsmoothed only"
  fi

  # AFNI writes 3dFWHMx.1D and .png in CWD. Ensure we run inside featdir to keep artifacts local, then delete them.
  pushd "$featdir" >/dev/null

  # UnsMoothed
  if [[ -f "$raw_img" ]]; then
    echo "   [${sub} ${ses} ${task} ${run}] 3dFWHMx (unsmoothed)…"
    rm -f 3dFWHMx.1D 3dFWHMx.1D.png
    3dFWHMx -detrend -ACF -mask "$mask" -input "$raw_img" > smoothness-0mm.txt
    append_to_tsv "$sub" "$ses" "$task" "$run" "0" "$featdir" "$raw_img" "smoothness-0mm.txt"
    rm -f 3dFWHMx.1D 3dFWHMx.1D.png
  fi

  # Smoothed
  if [[ -n "$smooth_img" && -f "$smooth_img" ]]; then
    smmm="$(kernel_from_name "$smooth_img")"
    echo "   [${sub} ${ses} ${task} ${run}] 3dFWHMx (smoothed ${smmm}mm)…"
    rm -f 3dFWHMx.1D 3dFWHMx.1D.png
    3dFWHMx -detrend -ACF -mask "$mask" -input "$smooth_img" > "smoothness-${smmm}mm.txt"
    append_to_tsv "$sub" "$ses" "$task" "$run" "${smmm}" "$featdir" "$smooth_img" "smoothness-${smmm}mm.txt"
    rm -f 3dFWHMx.1D 3dFWHMx.1D.png
  fi

  popd >/dev/null

done < <(find "$FSL_DERIV" -type d -name "*.feat" -print0)

echo "[$(date '+%F %T')] Done."
