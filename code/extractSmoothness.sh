#!/usr/bin/env bash
# Extract accurate smoothness (FWHM_eff from 3dFWHMx -acf) for raw/smoothed fMRIPrep BOLD
# from L1 FEAT directories, tasks mid/sharedreward only.
# De-duplicates across FEATs and image paths; MNI-only by default.

set -u
set -o pipefail

# ---------------- paths ----------------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"
FSL_DERIV="${rootdir}/derivatives/fsl"
FMRIPREP_DERIV="${rootdir}/derivatives/fmriprep"

# Output TSV at <project>/extractions/
OUT_DIR="${rootdir}/derivatives/extractions"
mkdir -p "${OUT_DIR}"
TSV="${OUT_DIR}/smoothness_acf.tsv"

# ---------------- tools ----------------
command -v 3dFWHMx >/dev/null || { echo "ERROR: 3dFWHMx (AFNI) not in PATH"; exit 1; }

# ---------------- env toggles ----------
VERBOSE="${VERBOSE:-0}"           # VERBOSE=1 for progress lines
ONLY_SPACE="${ONLY_SPACE:-MNI}"   # MNI | T1w | ALL

log() { [[ "$VERBOSE" == "1" ]] && echo "$*"; }

# ---------------- helpers --------------
parse_meta_from_path() {
  local p="$1"
  local sub="" ses="" task="" run="" space=""
  [[ "$p" =~ sub-([A-Za-z0-9]+) ]] && sub="${BASH_REMATCH[1]}"
  [[ "$p" =~ ses-([A-Za-z0-9]+) ]] && ses="${BASH_REMATCH[1]}"
  [[ "$p" =~ task-([^_/]+)      ]] && task="${BASH_REMATCH[1]}"
  [[ "$p" =~ run-([0-9]+)       ]] && run="${BASH_REMATCH[1]}"
  if [[ "$p" =~ space-([Tt]1[wW]) ]]; then space="T1w"; else space="MNI152NLin6Asym"; fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$sub" "$ses" "$task" "$run" "$space"
}

kernel_from_name() {
  local f="$1"
  if [[ "$f" =~ _([0-9]+)mm\.nii\.gz$ ]]; then echo "${BASH_REMATCH[1]}"; else echo "sm"; fi
}

append_row() {
  local sub="$1" ses="$2" task="$3" run="$4" acq="$5" kernel="$6" txt="$7"
  local feff="NA"
  if [[ -s "$txt" ]]; then feff="$(awk 'NR==2{print $4}' "$txt" 2>/dev/null || echo NA)"; fi
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

process_img() {
  # De-duplicate by image path + acq + kernel across the whole run
  local sub="$1" ses="$2" task="$3" run="$4" acq="$5" kernel="$6" img="$7" mask="$8" out_txt="$9"
  local key="${acq}|${kernel}|${img}"
  if [[ -n "${DONE_IMG[$key]:-}" ]]; then
    log "skip duplicate img: ${key}"
    return 0
  fi
  if 3dFWHMx -detrend -acf NULL -mask "$mask" -input "$img" > "$out_txt" 2>/dev/null; then
    append_row "$sub" "$ses" "$task" "$run" "$acq" "$kernel" "$out_txt"
    DONE_IMG[$key]=1
    return 0
  else
    log "warn: 3dFWHMx failed for $img"
    return 1
  fi
}

# header
if [[ ! -f "${TSV}" ]]; then
  echo -e "sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff" > "${TSV}"
fi

echo "[$(date +'%F %T')] Scanning L1 FEATs; writing to ${TSV}"

shopt -s nullglob
declare -A SEEN_RUN
declare -A DONE_IMG

rows=0 feats=0 skipped=0

for pat in \
  "${FSL_DERIV}"/sub-*/ses-*/L1_*.feat \
  "${FSL_DERIV}"/sub-*/L1_*.feat
do
  for featdir in $pat; do
    [[ -d "$featdir" ]] || continue

    # hard skips
    [[ "$featdir" == *"model-LSS"* ]] && continue
    [[ "$featdir" == *"/subject-level/"* ]] && continue
    [[ "$featdir" == *".gfeat"* ]] && continue
    [[ "$featdir" == *"/L2_"* ]] && continue
    [[ "$featdir" =~ task-(mid|sharedreward) ]] || continue

    # space selection
    if [[ "$ONLY_SPACE" == "MNI" ]]; then
      [[ "$featdir" =~ space-[Mm][Nn][Ii] ]] || continue
    elif [[ "$ONLY_SPACE" == "T1w" ]]; then
      [[ "$featdir" =~ space-[Tt]1[wW] ]] || continue
    fi

    mask="${featdir}/mask.nii.gz"
    [[ -f "$mask" ]] || { ((skipped++)); log "skip (no mask): $featdir"; continue; }

    IFS=$'\t' read -r sub ses task run space <<<"$(parse_meta_from_path "$featdir")"
    [[ -n "$sub" && -n "$task" && -n "$run" ]] || { ((skipped++)); log "skip (parse): $featdir"; continue; }

    # run-level de-duplication: only handle each (sub|ses|task|run|space) once
    run_key="${sub}|${ses}|${task}|${run}|${space}"
    if [[ -n "${SEEN_RUN[$run_key]:-}" ]]; then
      log "skip duplicate FEAT for run: $run_key"
      continue
    fi
    SEEN_RUN[$run_key]=1
    ((feats++))
    log "FEAT: sub=$sub ses=$ses task=$task run=$run space=$space"

    # reset per-run image memo
    DONE_IMG=()

    # find raw/smoothed per acq and process
    while IFS=$'\t' read -r acq raw_img smooth_img; do
      [[ -n "$raw_img" ]] || continue

      out_raw="${featdir}/smoothness-0mm_${acq}.txt"
      if process_img "$sub" "$ses" "$task" "$run" "$acq" "0" "$raw_img" "$mask" "$out_raw"; then
        ((rows++))
      fi

      if [[ -n "${smooth_img:-}" && -f "$smooth_img" ]]; then
        smmm="$(kernel_from_name "$smooth_img")"
        out_sm="${featdir}/smoothness-${smmm}mm_${acq}.txt"
        if process_img "$sub" "$ses" "$task" "$run" "$acq" "$smmm" "$smooth_img" "$mask" "$out_sm"; then
          ((rows++))
        fi
      fi
    done < <(find_inputs_for_run "$sub" "$ses" "$task" "$run" "$space")
  done
done

echo "[$(date +'%F %T')] Done. FEATs considered: ${feats}, rows written: ${rows}, skipped: ${skipped}"
