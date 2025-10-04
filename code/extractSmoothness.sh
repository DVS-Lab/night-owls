#!/usr/bin/env bash
# AFNI 3dFWHMx -acf smoothness for raw/smoothed fMRIPrep BOLD, per L1 FEAT, both acquisition types.
# Minimal progress output; TSV only includes FWHM_eff (gaussian_NEWmodel).

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
  echo -e "sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff" > "${TSV}"
fi

# -------- helpers --------

# Parse sub/ses/task/run/space from FEAT path (space = T1w or MNI152NLin6Asym)
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
    space="MNI152NLin6Asym"
  fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$sub" "$ses" "$task" "$run" "$space"
}

# Extract kernel size from smoothed filename suffix
kernel_from_name() {
  local f="$1"
  if [[ "$f" =~ _([0-9]+)mm\.nii\.gz$ ]]; then echo "${BASH_REMATCH[1]}"; else echo "sm"; fi
}

# Append a single row (only FWHM_eff) to TSV
append_to_tsv() {
  local sub="$1" ses="$2" task="$3" run="$4" acq="$5" kernel="$6" txt="$7"
  local feff="NA"
  if [[ -s "$txt" ]]; then
    feff="$(awk 'NR==2{print $4}' "$txt" 2>/dev/null || echo NA)"
  fi
  echo -e "${sub}\t${ses}\t${task}\t${run}\t${acq}\t${kernel}\t${feff}" >> "${TSV}"
}

# For a given sub/ses/task/run/space, emit up to two lines:
# "single-echo<TAB>RAW<TAB>SMOOTH" and "multiecho<TAB>RAW<TAB>SMOOTH"
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

  # single-echo (no part-mag)
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

  # multiecho (prefer combined part-mag; else echo-2; else any echo)
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

# -------- main --------
echo "[$(date '+%F %T')] Starting smoothness extraction under: ${FSL_DERIV}"
shopt -s nullglob

# Only L1 FEATs; explicitly SKIP LSS trial FEATs (model-LSS)
while IFS= read -r -d '' featdir; do
  [[ "$featdir" == *".gfeat"* || "$featdir" == *"/L2_"* || "$featdir" == *"/cope"* ]] && continue
  [[ "$featdir" == *"model-LSS"* ]] && continue

  mask="${featdir}/mask.nii.gz"
  [[ -f "$mask" ]] || { echo "WARN: No mask at ${mask}; skipping ${featdir}"; continue; }

  IFS=$'\t' read -r sub ses task run space <<<"$(parse_meta_from_path "$featdir")"
  echo "FEAT: sub=${sub} ses=${ses} task=${task} run=${run} space=${space}"
  [[ -n "$sub" && -n "$task" && -n "$run" ]] || { echo "WARN: Could not parse sub/task/run; skipping"; continue; }

  while IFS=$'\t' read -r acq raw_img smooth_img; do
    [[ -n "$raw_img" ]] || continue

    # UNSMOOTHED (quiet numeric; no plots), suppress AFNI chatter to stderr
    out_raw="${featdir}/smoothness-0mm_${acq}.txt"
    3dFWHMx -detrend -acf NULL -mask "$mask" -input "$raw_img" > "$out_raw" 2>/dev/null
    append_to_tsv "$sub" "$ses" "$task" "$run" "$acq" "0" "$out_raw"

    # SMOOTHED (if present)
    if [[ -n "${smooth_img:-}" && -f "$smooth_img" ]]; then
      smmm="$(kernel_from_name "$smooth_img")"
      out_sm="${featdir}/smoothness-${smmm}mm_${acq}.txt"
      3dFWHMx -detrend -acf NULL -mask "$mask" -input "$smooth_img" > "$out_sm" 2>/dev/null
      append_to_tsv "$sub" "$ses" "$task" "$run" "$acq" "$smmm" "$out_sm"
    fi
  done < <(find_fmriprep_for_run "$sub" "$ses" "$task" "$run" "$space")

done < <(find "$FSL_DERIV" -type d -path "*/L1_*.feat" -not -name "*model-LSS*" -print0)

echo "[$(date '+%F %T')] Done."
