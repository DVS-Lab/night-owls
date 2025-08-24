#!/usr/bin/env bash
set -euo pipefail

# Where this script lives: .../code
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"

# Inputs (relative to repo layout you've been using)
FSL_DERIV="${rootdir}/derivatives/fsl"
FMRIPREP_DERIV="${rootdir}/derivatives/fmriprep"
MASKS_DIR="${rootdir}/masks"
OUT_DIR="${rootdir}/derivatives/extractions"
mkdir -p "$OUT_DIR"

# Masks/maps (in MNI space by construction)
BRS_MNI="${MASKS_DIR}/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"
VS_MNI="${MASKS_DIR}/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

# Check required tools
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH." >&2; exit 1; }; }
for c in fslstats fslcc flirt; do need_cmd "$c"; done
have_ants=1
command -v antsApplyTransforms >/dev/null 2>&1 || have_ants=0

# Map zstat number -> human-readable label for each task
contrast_label() {
  local task="$1" z="$2"
  case "$task" in
    mid)
      case "$z" in
        7)  echo "anticipation_reward>neutral" ;;
        8)  echo "positive>negative" ;;
        9)  echo "reward:pos>neg" ;;
        10) echo "neutral:pos>neg" ;;
        *)  echo "" ;;
      esac
      ;;
    sharedreward)
      case "$z" in
        9)  echo "stranger>comp" ;;
        10) echo "neu>pun" ;;
        11) echo "rew>pun" ;;
        12) echo "S_rew>pun" ;;
        13) echo "C_rew>pun" ;;
        14) echo "S-C_rew>pun" ;;
        15) echo "rew>neu" ;;
        *)  echo "" ;;
      esac
      ;;
    *) echo "" ;;
  esac
}

# Resample/warp a source (MNI) image into the *exact* grid of a target image.
# If target is MNI grid, do identity apply with flirt. If target is T1w grid, use ANTs MNI->T1w xfm if available.
# Usage: resample_mni_to_target <sub> <ses> <target_img> <out_img> <src_mni_img>
resample_mni_to_target() {
  local sub="$1" ses="$2" target="$3" out="$4" src="$5"
  local space
  if [[ "$target" == *"_space-mni_"* ]] || [[ "$target" == *"_space-MNI"* ]]; then
    space="MNI"
  else
    space="T1w"
  fi

  if [[ "$space" == "MNI" ]]; then
    flirt -in "$src" -ref "$target" -applyxfm -usesqform -interp trilinear -out "$out"
  else
    local anatdir="${FMRIPREP_DERIV}/sub-${sub}/ses-${ses}/anat"
    local h5=
    for cand in \
      "${anatdir}/sub-${sub}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
      "${anatdir}/sub-${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
      "${anatdir}/sub-${sub}_ses-${ses}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5" \
      "${anatdir}/sub-${sub}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
    do
      [[ -f "$cand" ]] && { h5="$cand"; break; }
    done

    if [[ $have_ants -eq 1 && -n "${h5}" ]]; then
      antsApplyTransforms -d 3 \
        -i "$src" \
        -r "$target" \
        -t "$h5" \
        -n Linear \
        -o "$out"
    else
        echo "can't do ANTS, so something is wrong..."
        exit
    fi
  fi
}

# Compute VS mean and BRS correlation for a single zstat image
# Usage: process_one <zimg> <featdir> <sub> <ses> <run> <task> <space> <acq> <confounds> <znum> <label>
process_one() {
  local zimg="$1" featdir="$2" sub="$3" ses="$4" run="$5" task="$6" space="$7" acq="$8" confounds="$9" znum="${10}" label="${11}"

  local tmp
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' RETURN

  local brainmask="${featdir}/mask.nii.gz"
  if [[ ! -f "$brainmask" ]]; then
    echo "WARN: ${brainmask} missing; skipping ${zimg}" >&2
    return
  fi

  local brs_res="${tmp}/brs_res.nii.gz"
  local vs_res="${tmp}/vs_res.nii.gz"

  resample_mni_to_target "$sub" "$ses" "$zimg" "$brs_res" "$BRS_MNI"
  resample_mni_to_target "$sub" "$ses" "$zimg" "$vs_res" "$VS_MNI"

  local vs_mean="NA"
  if [[ -f "$vs_res" ]]; then
    vs_mean="$(fslstats "$zimg" -k "$vs_res" -M 2>/dev/null || echo NA)"
  fi

  local brs_corr="NA"
  if [[ -f "$brs_res" ]]; then
    # Whole-brain mask from FEAT; keep sign; include full range (-1..1)
    brs_corr="$(fslcc -m "$brainmask" --noabs -t -1 "$zimg" "$brs_res" 2>/dev/null | awk '{print $NF}' || echo NA)"
  fi

  echo -e "${sub}\t${ses}\t${run}\t${task}\t${space}\t${acq}\t${confounds}\t${znum}\t${label}\t${vs_mean}\t${brs_corr}"
}

main() {
  local out="${OUT_DIR}/extractions_L1stats.tsv"
  echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tVS_mean\tBRS_corr" > "$out"

  local glob="${GLOB:-}"
  shopt -s nullglob
  local arr=()
  if [[ -n "$glob" ]]; then
    while IFS= read -r -d '' f; do arr+=("$f"); done < <(find "$FSL_DERIV" -type f -path "*/L1_*/*.feat/stats/zstat*.nii.gz" -path "*${glob}*" -print0)
  else
    while IFS= read -r -d '' f; do arr+=("$f"); done < <(find "$FSL_DERIV" -type f -path "*/L1_*/*.feat/stats/zstat*.nii.gz" -print0)
  fi

  for zimg in "${arr[@]}"; do
    local featdir; featdir="$(dirname "$(dirname "$zimg")")"
    local subdir sesdir; subdir="$(basename "$(dirname "$(dirname "$featdir")")")"; sesdir="$(basename "$(dirname "$featdir")")"
    local sub="${subdir#sub-}"; local ses="${sesdir#ses-}"

    local featbase; featbase="$(basename "$featdir")"; featbase="${featbase%.feat}"

    local task run space_raw acq_raw confounds znum
    task="$(sed -E 's/^.*_task-([^_]+).*$/\1/' <<<"$featbase")"
    run="$(sed -E 's/^.*_run-([0-9]+).*$/\1/' <<<"$featbase")"
    space_raw="$(sed -E 's/^.*_space-([^_]+).*$/\1/' <<<"$featbase")"
    acq_raw="$(sed -E 's/^.*_(multi-echo|single-echo)_.*$/\1/' <<<"$featbase")"
    confounds="$(sed -E 's/^.*_cnfds-([^_]+).*$/\1/' <<<"$featbase")"
    znum="$(basename "$zimg" | sed -E 's/^zstat([0-9]+).*$/\1/')"

    local space acq
    if [[ "$space_raw" =~ ^(mni|MNI) ]]; then space="MNI152NLin6Asym"; else space="T1w"; fi
    if [[ "$acq_raw" == "multi-echo" ]]; then acq="multiecho"; else acq="single"; fi

    local label; label="$(contrast_label "$task" "$znum")"
    [[ -z "$label" ]] && continue

    process_one "$zimg" "$featdir" "$sub" "$ses" "$run" "$task" "$space" "$acq" "$confounds" "$znum" "$label" >> "$out"
  done

  echo "[${EPOCHREALTIME%.*}] Done. Wrote: $out"
}

main "$@"
