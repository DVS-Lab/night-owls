#!/usr/bin/env bash
# AFNI 3dFWHMx -acf (ACF model) smoothness for raw/smoothed fMRIPrep BOLD, per L1 FEAT, both acquisition types.
# Progress is echoed for each FEAT and for each acquisition (single-echo, multiecho).

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
  # include acq and fwhm_eff (the gaussian_NEWmodel FWHM)
  echo -e "sub\tses\ttask\trun\tacq\tkernel_mm\tacf_a\tacf_b\tacf_c\tfwhm_eff\tfwhm_x\tfwhm_y\tfwhm_z\tfeatdir\timg" > "${TSV}"
fi

# -------- helpers --------

# Parse values from 3dFWHMx -acf NULL stdout file ($txt) and append to TSV
# - line 1: classic FWHM_x FWHM_y FWHM_z FWHM_combined
# - line 2: a b c FWHM_eff
append_to_tsv() {
  local sub="$1" ses="$2" task="$3" run="$4" kernel="$5" featdir="$6" img="$7" txt="$8" acq="$9"
  local a="NA" b="NA" c="NA" feff="NA" fx="NA" fy="NA" fz="NA"
  if [[ -s "$txt" ]]; then
    # classic FWHM line
    read fx fy fz _ < <(awk 'NR==1{print; exit}' "$txt" 2>/dev/null || echo)
    # ACF params + effective FWHM
    read a b c feff < <(awk 'NR==2{print; exit}' "$txt" 2>/dev/null || echo)
  fi
  echo -e "${sub}\t${ses}\t${task}\t${run}\t${acq}\t${kernel}\t${a}\t${b}\t${c}\t${feff}\t${fx}\t${fy}\t${fz}\t${featdir}\t${img}" >> "${TSV}"
}

# Parse sub/ses/task/run/space from FEAT path
parse_meta_from_path() {
  local p="$1"
  local sub="" ses="" task="" run="" space=""
  [[ "$p" =~ sub-([A-Za-z0-9]+) ]] && sub="${BASH_REMATCH[1]}"
  [[ "$p" =~ ses-([A-Za-z0-9]+) ]] && ses="${BASH_REMATCH[1]}"
  [[ "$p" =~ task-([^_/]+)      ]] && task="${BASH_REMATCH[1]}"
  [[ "$p" =~ run-([0-9]+)       ]] && run="${BASH_REMATCH[1]}"
  if [[ "$p" =~ space-([Tt]1[wW]) ]]; then
    space="T1w"
  else
    space="MNI152NLin6Asym"   # default if FEAT name says mni/space-mni
  fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$sub" "$ses" "$task" "$run" "$space"
}

# For a given sub/ses/task/run/space, emit up to two lines:
#   "single-echo<TAB>RAW<TAB>SMOOTH"
#   "multiecho<TAB>RAW<TAB>SMOOTH"
find_fmriprep_for_run() {
  local sub="$1" ses="$2" task="$3" run="$4" space="$5"
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

  # ------- single-echo (no part-mag) -------
  local se_raw=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_space-${space}_desc-preproc_bold.nii.gz"
  )
  local se_sm=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_space-${space}_desc-preproc_bold_"*mm.nii.gz"
  )

  # ------- multiecho (prefer combined part-mag; else echo-2; else any echo) -------
  local me_raw=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-"*"_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-"*"_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-${space}_desc-preproc_bold.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-${space}_desc-preproc_bold.nii.gz"
  )
  local me_sm=(
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-2_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-"*"_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_echo-"*"_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
    "${funcdir}/sub-${sub}_task-${task}_run-${run}_acq-"*"_echo-"*"_part-mag_space-${space}_desc-preproc_bold_"*mm.nii.gz"
  )

  local pick_first=""
  pick_first_file() { pick_first=""; for f in "$@"; do [[ -f "$f" ]] && { pick_first="$f"; break; }; done; }

  pick_first_file "${se_raw[@]}"; local se_r="$pick_first"
  pick_first_file "${se_sm[@]}";  local se_s="$pick_first"
  if [[ -n "$se_r" ]]; then echo -e "single-echo\t${se_r}\t${se_s}"; fi

  pick_first_file "${me_raw[@]}"; local me_r="$pick_first"
  pick_first_file "${me_sm[@]}";  local me_s="$pick_first"
  if [[ -n "$me_r" ]]; then echo -e "multiecho\t${me_r}\t${me_s}"; fi
}

# Kernel from filename suffix
kernel_from_name() {
  local f="$1"
  if [[ "$f" =~ _([0-9]+)mm\.nii\.gz$ ]]; then echo "${BASH_REMATCH[1]}"; else echo "sm"; fi
}

# -------- main --------
echo "[$(date '+%F %T')] Starting smoothness extraction under: ${FSL_DERIV}"
shopt -s nullglob

# Iterate only L1 FEAT dirs (each FEAT supplies the mask and naming context)
while IFS= read -r -d '' featdir; do
  # Skip aggregate dirs just in case
  [[ "$featdir" == *".gfeat"* || "$featdir" == *"/L2_"* || "$featdir" == *"/cope"* ]] && continue

  mask="${featdir}/mask.nii.gz"
  [[ -f "$mask" ]] || { echo "WARN: No mask at ${mask}; skipping ${featdir}"; continue; }

  IFS=$'\t' read -r sub ses task run space <<<"$(parse_meta_from_path "$featdir")"
  echo ">> FEAT: sub=${sub} ses=${ses} task=${task} run=${run} space=${space} :: $(basename "$featdir")"
  [[ -n "$sub" && -n "$task" && -n "$run" ]] || { echo "WARN: Could not parse sub/task/run; skipping"; continue; }

  # Resolve raw/smoothed inputs for BOTH acquisition types, space-matched
  while IFS=$'\t' read -r acq raw_img smooth_img; do
    [[ -n "$raw_img" ]] || continue
    echo "   >> acq=${acq} raw=$(basename "$raw_img") smooth=$(basename "${smooth_img:-NONE}")"

    pushd "$featdir" >/dev/null

    # UNSMOOTHED: use -acf NULL to suppress 1D/PNG files; stdout is pure numeric (2 lines)
    3dFWHMx -detrend -acf NULL -mask "$mask" -input "$raw_img" > "smoothness-0mm_${acq}.txt"
    append_to_tsv "$sub" "$ses" "$task" "$run" "0" "$featdir" "$raw_img" "smoothness-0mm_${acq}.txt" "$acq"

    # SMOOTHED
    if [[ -n "${smooth_img:-}" && -f "$smooth_img" ]]; then
      smmm="$(kernel_from_name "$smooth_img")"
      3dFWHMx -detrend -acf NULL -mask "$mask" -input "$smooth_img" > "smoothness-${smmm}mm_${acq}.txt"
      append_to_tsv "$sub" "$ses" "$task" "$run" "${smmm}" "$featdir" "$smooth_img" "smoothness-${smmm}mm_${acq}.txt" "$acq"
    fi

    popd >/dev/null
  done < <(find_fmriprep_for_run "$sub" "$ses" "$task" "$run" "$space")

done < <(find "$FSL_DERIV" -type d -path "*/L1_*.feat" -print0)

echo "[$(date '+%F %T')] Done."
