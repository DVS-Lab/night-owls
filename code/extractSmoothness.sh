#!/usr/bin/env bash
# Extract accurate smoothness (FWHM_eff from 3dFWHMx -acf) for raw/smoothed fMRIPrep BOLD
# from L1 FEAT directories, for tasks mid and sharedreward only.

set -u  # no -e so a single failure won’t stop everything
set -o pipefail

# ---------------- paths ----------------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"
FSL_DERIV="${rootdir}/derivatives/fsl"
FMRIPREP_DERIV="${rootdir}/derivatives/fmriprep"

# Output TSV lives at <project>/extractions/ (to match your workflow)
OUT_DIR="${rootdir}/derivatives/extractions"
mkdir -p "${OUT_DIR}"
TSV="${OUT_DIR}/smoothness_acf.tsv"

# ---------------- tools ----------------
command -v 3dFWHMx >/dev/null || { echo "ERROR: 3dFWHMx (AFNI) not in PATH"; exit 1; }

# ---------------- env toggles ----------
VERBOSE="${VERBOSE:-0}"  # set VERBOSE=1 for progress lines

# ---------------- helpers --------------
log() { [[ "$VERBOSE" == "1" ]] && echo "$*"; }

parse_meta_from_path() {
  # prints: sub \t ses \t task \t run \t space
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

kernel_from_name() {
  local f="$1"
  if [[ "$f" =~ _([0-9]+)mm\.nii\.gz$ ]]; then echo "${BASH_REMATCH[1]}"; else echo "sm"; fi
}

append_row() {
  # only FWHM_eff to TSV
  local sub="$1" ses="$2" task="$3" run="$4" acq="$5" kernel="$6" txt="$7"
  local feff="NA"
  if [[ -s "$txt" ]]; then
    # 3dFWHMx -acf NULL prints 2 numeric lines; FWHM_eff is field 4 on line 2
    feff="$(awk 'NR==2{print $4}' "$txt" 2>/dev/null || echo NA)"
  fi
  echo -e "${sub}\t${ses}\t${task}\t${run}\t${acq}\t${kernel}\t${feff}" >> "${TSV}"
}

find_inputs_for_run() {
  # emits up to two lines:
  #   "single-echo<TAB>RAW<TAB>SMOOTH"
  #   "multiecho<TAB>RAW<TAB>SMOOTH"
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

# header
if [[ ! -f "${TSV}" ]]; then
  echo -e "sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff" > "${TSV}"
fi

echo "[$(date +'%F %T')] Scanning L1 FEATs; writing to ${TSV}"

shopt -s nullglob
rows=0 feats=0 skipped=0

# Use shell globs to avoid find() touching L2 or gfeat trees (and spamming Permission denied)
# Common on your tree: sub-*/ses-*/L1_*.feat plus some L1_*.feat right under sub-*/ses-*
for pat in \
  "${FSL_DERIV}"/sub-*/ses-*/L1_*.feat \
  "${FSL_DERIV}"/sub-*/L1_*.feat
do
  for featdir in $pat; do
    # hard skips
    [[ -d "$featdir" ]] || continue
    [[ "$featdir" == *"model-LSS"* ]] && continue
    [[ "$featdir" == *"/subject-level/"* ]] && continue
    [[ "$featdir" == *".gfeat"* ]] && continue
    [[ "$featdir" == *"/L2_"* ]] && continue

    # task filter: only mid and sharedreward
    [[ "$featdir" =~ task-(mid|sharedreward) ]] || continue

    mask="${featdir}/mask.nii.gz"
    [[ -f "$mask" ]] || { ((skipped++)); log "skip (no mask): $featdir"; continue; }

    IFS=$'\t' read -r sub ses task run space <<<"$(parse_meta_from_path "$featdir")"
    [[ -n "$sub" && -n "$task" && -n "$run" ]] || { ((skipped++)); log "skip (parse): $featdir"; continue; }
    ((feats++))
    log "FEAT: sub=$sub ses=$ses task=$task run=$run space=$space"

    # find raw/smoothed for each acquisition (single-echo, multiecho)
    while IFS=$'\t' read -r acq raw_img smooth_img; do
      [[ -n "$raw_img" ]] || continue

      # UNSMOOTHED
      out_raw="${featdir}/smoothness-0mm_${acq}.txt"
      if 3dFWHMx -detrend -acf NULL -mask "$mask" -input "$raw_img" > "$out_raw" 2>/dev/null; then
        append_row "$sub" "$ses" "$task" "$run" "$acq" "0" "$out_raw"
        ((rows++))
      else
        log "warn: 3dFWHMx failed (raw) for $featdir [$acq]"
      fi

      # SMOOTHED (if present)
      if [[ -n "${smooth_img:-}" && -f "$smooth_img" ]]; then
        smmm="$(kernel_from_name "$smooth_img")"
        out_sm="${featdir}/smoothness-${smmm}mm_${acq}.txt"
        if 3dFWHMx -detrend -acf NULL -mask "$mask" -input "$smooth_img" > "$out_sm" 2>/dev/null; then
          append_row "$sub" "$ses" "$task" "$run" "$acq" "$smmm" "$out_sm"
          ((rows++))
        else
          log "warn: 3dFWHMx failed (sm) for $featdir [$acq ${smmm}mm]"
        fi
      fi
    done < <(find_inputs_for_run "$sub" "$ses" "$task" "$run" "$space")
  done
done

echo "[$(date +'%F %T')] Done. FEATs: ${feats}, rows written: ${rows}, skipped: ${skipped}"