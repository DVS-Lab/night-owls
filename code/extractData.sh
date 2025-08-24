#!/usr/bin/env bash
set -euo pipefail

# -------- fixed locations (relative to THIS script), and required tools --------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"                 # project root (…/night-owls)
FSL_DERIV="${rootdir}/derivatives/fsl"
FMRIPREP_DERIV="${rootdir}/derivatives/fmriprep"
MASKS_DIR="${rootdir}/masks"
OUT_DIR="${rootdir}/derivatives/extractions"

BRS_MNI="${MASKS_DIR}/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"
VS_MNI="${MASKS_DIR}/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

command -v antsApplyTransforms >/dev/null || { echo "ERROR: antsApplyTransforms not found in PATH"; exit 1; }
command -v fslstats >/dev/null || { echo "ERROR: fslstats not found in PATH"; exit 1; }
command -v fslcc >/dev/null    || { echo "ERROR: fslcc not found in PATH"; exit 1; }

[[ -d "$FSL_DERIV" ]] || { echo "ERROR: Can't find ${FSL_DERIV}"; exit 1; }
mkdir -p "$OUT_DIR"

# -------- contrast label maps (from your screenshots) --------
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

# -------- warp/resample helper (NO FLIRT) --------
# Resample a MNI-space mask/map to the exact grid of a target zstat image.
# - If zstat is in MNI space: identity into zstat grid (no transform).
# - If zstat is in T1w space: apply fMRIPrep MNI→T1w H5 (ses-specific if present).
mni_to_target() {
  local sub="$1" ses="$2" target_img="$3" src_mni="$4" out_img="$5"

  # detect target space from FEAT name or header path (we pass space separately to avoid header reads)
  local space_hint="$6"  # "mni" or "t1w"
  if [[ "$space_hint" == "mni" ]]; then
    antsApplyTransforms -d 3 -i "$src_mni" -r "$target_img" -n Linear -o "$out_img"
  else
    local anat="${FMRIPREP_DERIV}/sub-${sub}/ses-${ses}/anat"
    local h5=""
    for cand in \
      "${anat}/sub-${sub}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
      "${anat}/sub-${sub}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5" \
      "${anat}/sub-${sub}_ses-${ses}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5" \
      "${anat}/sub-${sub}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
    do
      [[ -f "$cand" ]] && { h5="$cand"; break; }
    done
    [[ -z "$h5" ]] && { echo "WARN: no MNI→T1w transform for sub-${sub} ses-${ses}; skipping $target_img" >&2; return 1; }
    antsApplyTransforms -d 3 -i "$src_mni" -r "$target_img" -t "$h5" -n Linear -o "$out_img"
  fi
}

# -------- per-image processing --------
process_one() {
  local zimg="$1" featdir="$2" sub="$3" ses="$4" run="$5" task="$6" space_tag="$7" acq="$8" confounds="$9" znum="${10}" label="${11}"

  local brainmask="${featdir}/mask.nii.gz"
  [[ -f "$brainmask" ]] || { echo "WARN: missing FEAT mask: $brainmask — skipping"; return; }

  local tmp; tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' RETURN


if [[ "$space_tag" == "mni" ]]; then
  brs_res="$BRS_MNI"   # use MNI map as-is
  vs_res="$VS_MNI"
else
  brs_res="$tmp/brs_in_T1w.nii.gz"   # define outputs first
  vs_res="$tmp/vs_in_T1w.nii.gz"
  mni_to_target "$sub" "$ses" "$zimg" "$BRS_MNI" "$brs_res" "t1w" || return
  mni_to_target "$sub" "$ses" "$zimg" "$VS_MNI"  "$vs_res"  "t1w" || return
fi



  # VS mean
  local vs_mean="NA"
  [[ -f "$vs_res" ]] && vs_mean="$(fslstats "$zimg" -k "$vs_res" -M 2>/dev/null || echo NA)"

  # Signed whole-brain spatial corr with BRS (mask to FEAT brainmask)
  local brs_corr="NA"
  if [[ -f "$brs_res" ]]; then
    brs_corr="$(fslcc -m "$brainmask" --noabs -t -1 "$zimg" "$brs_res" 2>/dev/null | awk '{print $NF}' || echo NA)"
  fi

  echo -e "${sub}\t${ses}\t${run}\t${task}\t${space_tag}\t${acq}\t${confounds}\t${znum}\t${label}\t${vs_mean}\t${brs_corr}"
}

# -------- main --------
out="${OUT_DIR}/extractions_L1stats.tsv"
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tVS_mean\tBRS_corr" > "$out"

shopt -s nullglob
# ONLY L1_task-… FEATs under derivatives/fsl; ignore everything else
while IFS= read -r -d '' zimg; do
  featdir="$(dirname "$(dirname "$zimg")")"                 # …/L1_task-….feat
  sesdir="$(basename "$(dirname "$featdir")")"              # ses-XX
  subdir="$(basename "$(dirname "$(dirname "$featdir")")")" # sub-XXX
  sub="${subdir#sub-}"; ses="${sesdir#ses-}"

  fbase="$(basename "$featdir")"; fbase="${fbase%.feat}"

  task="$(sed -E 's/^.*_task-([^_]+).*$/\1/' <<<"$fbase")"
  run="$(sed -E 's/^.*_run-([0-9]+).*$/\1/' <<<"$fbase")"
  space_raw="$(sed -E 's/^.*_space-([^_]+).*$/\1/' <<<"$fbase")"
  acq_raw="$(sed -E 's/^.*_(multi-echo|single-echo)_.*$/\1/' <<<"$fbase")"
  confounds="$(sed -E 's/^.*_cnfds-([^_]+).*$/\1/' <<<"$fbase")"
  znum="$(basename "$zimg" | sed -E 's/^zstat([0-9]+).*$/\1/')"

  # normalize tags used in output and in mni_to_target's space hint
  space_tag="t1w"; [[ "$space_raw" =~ ^(mni|MNI) ]] && space_tag="mni"
  acq="single"; [[ "$acq_raw" == "multi-echo" ]] && acq="multiecho"

  label="$(contrast_label "$task" "$znum")"
  [[ -z "$label" ]] && continue  # only extract requested contrasts

  process_one "$zimg" "$featdir" "$sub" "$ses" "$run" "$task" "$space_tag" "$acq" "$confounds" "$znum" "$label" >> "$out"
done < <(find "$FSL_DERIV" -type f -path "*/L1_task-*/stats/zstat*.nii.gz" -print0)

echo "[$(date '+%F %T')] Done. Wrote: $out"
