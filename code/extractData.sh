#!/usr/bin/env bash
set -euo pipefail

# -------------------- fixed project-relative paths --------------------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"                 # project root
FSL_DERIV="${rootdir}/derivatives/fsl"
maskdir="${rootdir}/masks"
OUT_DIR="${rootdir}/derivatives/extractions"

# -------------------- MNI masks/maps (as provided) --------------------
NAcc_MNI="$maskdir/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
CORT_MNI="$maskdir/BRS_Cortical_3pt1.nii.gz"
BRS_MNI="$maskdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

# -------------------- required tools --------------------
command -v fslstats >/dev/null || { echo "ERROR: fslstats not found in PATH"; exit 1; }
command -v fslcc    >/dev/null || { echo "ERROR: fslcc not found in PATH"; exit 1; }

# -------------------- validate fixed paths --------------------
[[ -d "$FSL_DERIV" ]] || { echo "ERROR: Can't find FSL derivatives dir: $FSL_DERIV"; exit 1; }
[[ -f "$NAcc_MNI" ]]  || { echo "ERROR: Missing NAcc mask: $NAcc_MNI"; exit 1; }
[[ -f "$CORT_MNI" ]]  || { echo "ERROR: Missing cortical mask: $CORT_MNI"; exit 1; }
[[ -f "$BRS_MNI" ]]   || { echo "ERROR: Missing BRS map: $BRS_MNI"; exit 1; }

mkdir -p "$OUT_DIR"

# -------------------- contrast label maps --------------------
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
    *)
      echo ""
      ;;
  esac
}

# -------------------- tag extraction from FEAT basename --------------------
# Example basename (no .feat):
# L1_sub-105_ses-12_task-sharedreward_model-1_type-act_run-1_space-mni_single-echo_cnfds-fmriprep_unsmoothed
get_tag() {
  local key="$1" str="$2"
  if [[ "$str" =~ _${key}-([^_]+) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "NA"
  fi
}

# -------------------- per-contrast processing --------------------
process_one() {
  local zimg="$1" copeimg="$2" varcopeimg="$3" featdir="$4" \
        sub="$5" ses="$6" run="$7" task="$8" acq="$9" confounds="${10}" \
        znum="${11}" label="${12}"

  local brainmask="${featdir}/mask.nii.gz"
  [[ -f "$brainmask" ]] || { echo "WARN: missing FEAT mask: $brainmask — skipping" >&2; return; }

  local NAcc_z="NA" NAcc_c="NA" NAcc_v="NA"
  local Cort_z="NA" Cort_c="NA" Cort_v="NA"

  [[ -f "$zimg" ]]       && NAcc_z="$(fslstats "$zimg"       -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$copeimg" ]]    && NAcc_c="$(fslstats "$copeimg"    -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$varcopeimg" ]] && NAcc_v="$(fslstats "$varcopeimg" -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"

  [[ -f "$zimg" ]]       && Cort_z="$(fslstats "$zimg"       -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$copeimg" ]]    && Cort_c="$(fslstats "$copeimg"    -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$varcopeimg" ]] && Cort_v="$(fslstats "$varcopeimg" -k "$CORT_MNI" -M 2>/dev/null || echo NA)"

  local brs_corr="NA"
  if [[ -f "$zimg" ]]; then
    brs_corr="$(fslcc -m "$brainmask" --noabs -t -1 -p 6 "$zimg" "$BRS_MNI" 2>/dev/null | awk '{print $NF}' || echo NA)"
  fi

  echo -e "${sub}\t${ses}\t${run}\t${task}\tmni\t${acq}\t${confounds}\t${znum}\t${label}\t${NAcc_z}\t${NAcc_c}\t${NAcc_v}\t${Cort_z}\t${Cort_c}\t${Cort_v}\t${brs_corr}"
}

# -------------------- main --------------------
out="${OUT_DIR}/extractions_L1stats-revised.tsv"
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tNAcc_zstat_mean\tNAcc_cope_mean\tNAcc_varcope_mean\tBRS_Cort_zstat_mean\tBRS_Cort_cope_mean\tBRS_Cort_varcope_mean\tBRS_corr" \
  > "$out"

# Strictly target UNSMOOTHED L1 FEATs
ZPATH="*/L1_*unsmoothed.feat/stats/zstat*.nii.gz"

total=$(
  find "$FSL_DERIV" \
    -type d \( -name '*.gfeat' -o -name 'subject-level' \) -prune -o \
    -type f -path "$ZPATH" -print0 \
  | tr -cd '\0' | wc -c
)

if (( total == 0 )); then
  echo "[WARN] No zstat files matched." >&2
  echo "[WARN] FSL_DERIV = $FSL_DERIV" >&2
  echo "[WARN] Pattern   = $ZPATH" >&2
fi

done_cnt=0

while IFS= read -r -d '' zimg; do
  ((done_cnt++))

  featdir="$(dirname "$(dirname "$zimg")")"                 # …/L1_...unsmoothed.feat
  sesdir="$(basename "$(dirname "$featdir")")"              # ses-XX
  subdir="$(basename "$(dirname "$(dirname "$featdir")")")" # sub-XXX

  sub="${subdir#sub-}"
  ses="${sesdir#ses-}"

  fbase="$(basename "$featdir")"
  fbase="${fbase%.feat}"

  task="$(get_tag task "$fbase")"
  run="$(get_tag run "$fbase")"
  space_raw="$(get_tag space "$fbase")"

  # acquisition + confounds tags from your naming convention
  # (these are not key-value tags in the same way, so we parse lightly)
  acq_raw="NA"
  if [[ "$fbase" == *_multi-echo_* ]]; then
    acq_raw="multi-echo"
  elif [[ "$fbase" == *_single-echo_* ]]; then
    acq_raw="single-echo"
  fi

  confounds="NA"
  if [[ "$fbase" =~ _cnfds-([^_]+) ]]; then
    confounds="${BASH_REMATCH[1]}"
  fi

  znum="$(basename "$zimg" | sed -E 's/^zstat([0-9]+).*$/\1/')"

  # normalize tags used in output
  space_tag="t1w"
  [[ "$space_raw" =~ ^(mni|MNI)$ ]] && space_tag="mni"

  acq="single"
  [[ "$acq_raw" == "multi-echo" ]] && acq="multiecho"

  # ignore non-MNI
  [[ "$space_tag" == "mni" ]] || continue

  label="$(contrast_label "$task" "$znum")"
  [[ -n "$label" ]] || continue

  statsdir="$(dirname "$zimg")"
  copeimg="${statsdir}/cope${znum}.nii.gz"
  varcopeimg="${statsdir}/varcope${znum}.nii.gz"

  if [[ ! -f "$copeimg" || ! -f "$varcopeimg" ]]; then
    echo "WARN: missing cope/varcope for $zimg (expected $copeimg and $varcopeimg)" >&2
  fi

  process_one "$zimg" "$copeimg" "$varcopeimg" "$featdir" \
    "$sub" "$ses" "$run" "$task" "$acq" "$confounds" "$znum" "$label" \
    >> "$out"

done < <(
  find "$FSL_DERIV" \
    -type d \( -name '*.gfeat' -o -name 'subject-level' \) -prune -o \
    -type f -path "$ZPATH" -print0
)

echo "[$(date '+%F %T')] Done. Wrote: $out"
