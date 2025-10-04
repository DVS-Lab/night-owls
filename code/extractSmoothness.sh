#!/usr/bin/env bash
# AFNI 3dFWHMx -ACF smoothness for raw/smoothed fMRIPrep BOLD, per L1 FEAT, both acquisition types.

set -euo pipefail

# -------- locations (relative to THIS script) --------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"
FSL_DERIV="${rootdir}/derivatives/fsl"
FMRIPREP_DERIV="${rootdir}/derivatives/fmriprep"

# -------- tools --------
command -v 3dFWHMx >/dev/null || { echo "ERROR: 3dFWHMx (AFNI) not found in PATH"; exit 1; }

# -------- aggregation output --------
OUT_DIR="${rootdir}/derivatives/extractions"
mkdir -p "${OUT_DIR}"
TSV="${OUT_DIR}/smoothness_acf.tsv"
if [[ ! -f "${TSV}" ]]; then
  echo -e "sub\tses\ttask\trun\tacq\tkernel_mm\tacf_a\tacf_b\tacf_c\tfwhm_x\tfwhm_y\tfwhm_z\tfeatdir\timg" > "${TSV}"
fi

# -------- helpers --------

# Parse values from AFNI output and append row to TSV
append_to_tsv() {
  local sub="$1" ses="$2" task="$3" run="$4" kernel="$5" featdir="$6" img="$7" txt="$8" acq_in="${9:-}"
  local acq="${acq_in:-auto}"
  local a="NA" b="NA" c="NA" fx="NA" fy="NA" fz="NA"

  # ACF parameters: last three numeric tokens on any 'ACF' line
  if grep -qi "ACF" "$txt"; then
    read a b c < <(grep -i "ACF" "$txt" | tail -n1 | grep -Eo "[-+]?[0-9]*\.?[0-9]+" | tail -n3)
  fi

  # FWHM xyz: last three numeric tokens on any 'FWHM' line
  if grep -qi "FWHM" "$txt"; then
    read fx fy fz < <(grep -i "FWHM" "$txt" | tail -n1 | grep -Eo "[-+]?[0-9]*\.?[0-9]+" | tail -n3)
  fi

  echo -e "${sub}\t${ses}\t${task}\t${run}\t${acq}\t${kernel}\t${a}\t${b}\t${c}\t${fx}\t${fy}\t${fz}\t${featdir}\t${img}" >> "${TSV}"
}

# Pull sub/ses/task/run from FEAT path
parse_meta_from_path() {
  local p="$1"
  local sub="" ses="" task="" run=""
  if [[ "$p" =~ sub-([A-Za-z0-9]+) ]]; then sub="${BASH_REMATCH[1]}"; fi
  if [[ "$p" =~ ses-([A-Za-z0-9]+) ]]; then ses="${BASH_REMATCH[1]}"; fi
  if [[ "$p" =~ task-([^_/]+) ]]; then task="${BASH_REMATCH[1]}"; fi
  if [[ "$p" =~ run-([0-9]+) ]]; then run="${BASH_REMATCH[1]}"; fi
  printf "%s\t%s\t%s\t%s\n" "$sub" "$ses" "$task" "$run"
}

# For a given sub/ses/task/run, emit up to two lines:
#   "single-echo<TAB>RAW<TAB>SMOOTH"
#   "multiecho<TAB>RAW<TAB>SMOOTH"
find_fmriprep_for_run() {
  local sub="$1" ses="$2" task="$3" run="$4"
  local base="${FMRIPREP_DERIV}/sub-${sub}"
  local funcdir=""
  if [[ -n "$ses" && -d "${base}/ses-${ses}/func" ]]; then
    funcdir="${base}/ses-${ses}/func"
  elif [[ -d "${base}/func" ]]; then
    funcdir="${base}/func"
  else
    return 0
  fi
  shopt -s nullglob

  # single-echo (no part-mag)
  local se_raw=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
  )
  local se_sm=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
  )

  # multiecho (prefer combined part-mag without echo tag; else echo-2; else any echo)
  local me_raw=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
  )
  local me_sm=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-MNI152NLin6Asym_desc-preproc_bold_"*mm.nii.gz"
  )

  local pick_first=""
  pick_first_file() {
    pick_first=""
    for f in "$@"; do
      if [[ -f "$f" ]]; then pick_first="$f"; break; fi
    done
  }

  pick_first_file "${se_raw[@]}"; local se_r="$pick_first"
  pick_first_file "${se_sm[@]}";  local se_s="$pick_first"
  if [[ -n "$se_r" ]]; then
    echo -e "single-echo\t${se_r}\t${se_s}"
  fi

  pick_first_file "${me_raw[@]}"; local me_r="$pick_first"
  pick_first_file "${me_sm[@]}";  local me_s="$pick_first"
  if [[ -n "$me_r" ]]; then
    echo -e "multiecho\t${me_r}\t${me_s}"
  fi
}

# Extract kernel mm from filename
kernel_from_name() {
  local f="$1"
  if [[ "$f" =~ _([0-9]+)mm\.nii\.gz$ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "sm"
  fi
}

# -------- main --------
echo "[$(date '+%F %T')] Starting smoothness extraction under: ${FSL_DERIV}"
shopt -s nullglob

# Iterate only L1 FEAT dirs (FIXED pattern)
while IFS= read -r -d '' featdir; do
  # Skip non-L1 or aggregate dirs
  if [[ "$featdir" == *".gfeat"* || "$featdir" == *"/L2_"* || "$featdir" == *"/cope"* ]]; then
    continue
  fi

  mask="${featdir}/mask.nii.gz"
  if [[ ! -f "$mask" ]]; then
    echo "WARN: No mask at ${mask}; skipping ${featdir}"
    continue
  fi

  IFS=$'\t' read -r sub ses task run <<<"$(parse_meta_from_path "$featdir")"
  # PROGRESS echo for each loop:
  echo ">> FEAT: sub=${sub} ses=${ses} task=${task} run=${run} :: $(basename "$featdir")"

  if [[ -z "$sub" || -z "$task" || -z "$run" ]]; then
    echo "WARN: Could not parse sub/task/run from ${featdir}; skipping"
    continue
  fi

  # For this run, process both acquisition types discovered
  while IFS=$'\t' read -r acq raw_img smooth_img; do
    [[ -z "$raw_img" ]] && continue

    # PROGRESS echo for resolved inputs:
    echo "   >> acq=${acq} raw=$(basename "$raw_img") smooth=$(basename "${smooth_img:-NONE}")"

    pushd "$featdir" >/dev/null

    # unsmoothed
    rm -f 3dFWHMx.1D 3dFWHMx.1D.png
    3dFWHMx -detrend -ACF -mask "$mask" -input "$raw_img" > "smoothness-0mm_${acq}.txt"
    append_to_tsv "$sub" "$ses" "$task" "$run" "0" "$featdir" "$raw_img" "smoothness-0mm_${acq}.txt" "$acq"
    rm -f 3dFWHMx.1D 3dFWHMx.1D.png

    # smoothed
    if [[ -n "$smooth_img" && -f "$smooth_img" ]]; then
      smmm="$(kernel_from_name "$smooth_img")"
      rm -f 3dFWHMx.1D 3dFWHMx.1D.png
      3dFWHMx -detrend -ACF -mask "$mask" -input "$smooth_img" > "smoothness-${smmm}mm_${acq}.txt"
      append_to_tsv "$sub" "$ses" "$task" "$run" "${smmm}" "$featdir" "$smooth_img" "smoothness-${smmm}mm_${acq}.txt" "$acq"
      rm -f 3dFWHMx.1D 3dFWHMx.1D.png
    fi

    popd >/dev/null
  done < <(find_fmriprep_for_run "$sub" "$ses" "$task" "$run")

done < <(find "$FSL_DERIV" -type d -path "*/L1_*.feat" -print0)

echo "[$(date '+%F %T')] Done."
