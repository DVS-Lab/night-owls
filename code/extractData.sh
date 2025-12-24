#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# -------- fixed locations (relative to THIS script) --------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"                 # project root (…/night-owls)

FSL_DERIV="${rootdir}/derivatives/fsl"
maskdir="${rootdir}/masks"
OUT_DIR="${rootdir}/derivatives/extractions"

# MNI masks/maps (use exactly what you provided)
NAcc_MNI="$maskdir/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
CORT_MNI="$maskdir/BRS_Cortical_3pt1.nii.gz"
BRS_MNI="$maskdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

# -------- required tools --------
command -v fslstats >/dev/null || { echo "ERROR: fslstats not found in PATH"; exit 1; }
command -v fslcc    >/dev/null || { echo "ERROR: fslcc not found in PATH"; exit 1; }

# -------- sanity checks --------
[[ -d "$FSL_DERIV" ]] || { echo "ERROR: Can't find ${FSL_DERIV}"; exit 1; }
mkdir -p "$OUT_DIR"

[[ -f "$BRS_MNI"  ]] || { echo "ERROR: Missing BRS map: $BRS_MNI"; exit 1; }
[[ -f "$NAcc_MNI" ]] || { echo "ERROR: Missing NAcc mask: $NAcc_MNI"; exit 1; }
[[ -f "$CORT_MNI" ]] || { echo "ERROR: Missing cortical mask: $CORT_MNI"; exit 1; }

# -------- contrast label maps --------
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

# -------- helpers --------
extract_tag() {
  # Extracts value for patterns like _task-XXXX_ from a string.
  # Usage: extract_tag task "$basename"
  local key="$1" str="$2"
  if [[ "$str" =~ _${key}-([^_]+) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "NA"
  fi
}

extract_acq() {
  local str="$1"
  if [[ "$str" =~ _multi-echo_ ]]; then
    echo "multiecho"
  elif [[ "$str" =~ _single-echo_ ]]; then
    echo "single"
  else
    # fallback to your earlier convention if naming drifts slightly
    echo "single"
  fi
}

# -------- per-contrast processing (MNI-only) --------
process_one() {
  local zimg="$1" copeimg="$2" varcopeimg="$3" featdir="$4" \
        sub="$5" ses="$6" run="$7" task="$8" acq="$9" confounds="${10}" \
        znum="${11}" label="${12}"

  local brainmask="${featdir}/mask.nii.gz"
  [[ -f "$brainmask" ]] || { echo "WARN: missing FEAT mask: $brainmask — skipping" >&2; return; }

  # ROI means
  local NAcc_z="NA" NAcc_c="NA" NAcc_v="NA"
  local Cort_z="NA" Cort_c="NA" Cort_v="NA"

  [[ -f "$zimg" ]]       && NAcc_z="$(fslstats "$zimg"       -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$copeimg" ]]    && NAcc_c="$(fslstats "$copeimg"    -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$varcopeimg" ]] && NAcc_v="$(fslstats "$varcopeimg" -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"

  [[ -f "$zimg" ]]       && Cort_z="$(fslstats "$zimg"       -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$copeimg" ]]    && Cort_c="$(fslstats "$copeimg"    -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$varcopeimg" ]] && Cort_v="$(fslstats "$varcopeimg" -k "$CORT_MNI" -M 2>/dev/null || echo NA)"

  # Signed whole-brain spatial corr with BRS (mask to FEAT brainmask) using zstat
  local brs_corr="NA"
  if [[ -f "$zimg" ]]; then
    brs_corr="$(fslcc -m "$brainmask" --noabs -t -1 -p 6 "$zimg" "$BRS_MNI" 2>/dev/null | awk '{print $NF}' || echo NA)"
  fi

  echo -e "${sub}\t${ses}\t${run}\t${task}\tmni\t${acq}\t${confounds}\t${znum}\t${label}\t${NAcc_z}\t${NAcc_c}\t${NAcc_v}\t${Cort_z}\t${Cort_c}\t${Cort_v}\t${brs_corr}"
}

# -------- main --------
out="${OUT_DIR}/extractions_L1stats-revised_smoothed.tsv"
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tNAcc_zstat_mean\tNAcc_cope_mean\tNAcc_varcope_mean\tBRS_Cort_zstat_mean\tBRS_Cort_cope_mean\tBRS_Cort_varcope_mean\tBRS_corr" \
  > "$out"

# Tight, permission-safe discovery: ONLY L1_*{fmriprep,tedana}.feat under sub-*/ses-*
zfiles=(
  "$FSL_DERIV"/sub-*/ses-*/L1_*fmriprep.feat/stats/zstat*.nii.gz
  "$FSL_DERIV"/sub-*/ses-*/L1_*tedana.feat/stats/zstat*.nii.gz
)

if (( ${#zfiles[@]} == 0 )); then
  echo "WARN: No files found." >&2
  echo "Tried patterns:" >&2
  echo "  $FSL_DERIV/sub-*/ses-*/L1_*fmriprep.feat/stats/zstat*.nii.gz" >&2
  echo "  $FSL_DERIV/sub-*/ses-*/L1_*tedana.feat/stats/zstat*.nii.gz" >&2
  echo "Check that you're running from the expected project tree and that L1 FEATs exist." >&2
  echo "Wrote header only to: $out" >&2
  exit 0
fi

for zimg in "${zfiles[@]}"; do
  featdir="$(dirname "$(dirname "$zimg")")"  # …/L1_...{fmriprep|tedana}.feat

  # Guard: ensure we only process desired FEAT dirs
  case "$featdir" in
    *fmriprep.feat|*tedana.feat) : ;;
    *) continue ;;
  esac

  sesdir="$(basename "$(dirname "$featdir")")"
  subdir="$(basename "$(dirname "$(dirname "$featdir")")")"

  sub="${subdir#sub-}"
  ses="${sesdir#ses-}"

  fbase="$(basename "$featdir")"
  fbase="${fbase%.feat}"

  task="$(extract_tag task "$fbase")"
  run="$(extract_tag run "$fbase")"
  space_raw="$(extract_tag space "$fbase")"
  confounds="$(extract_tag cnfds "$fbase")"
  acq="$(extract_acq "$fbase")"

  # Only MNI space in this script
  [[ "$space_raw" =~ ^(mni|MNI) ]] || continue

  znum="$(basename "$zimg" | sed -E 's/^zstat([0-9]+).*$/\1/')"
  [[ "$znum" =~ ^[0-9]+$ ]] || continue

  label="$(contrast_label "$task" "$znum")"
  [[ -z "$label" ]] && continue

  statsdir="$(dirname "$zimg")"
  copeimg="${statsdir}/cope${znum}.nii.gz"
  varcopeimg="${statsdir}/varcope${znum}.nii.gz"

  if [[ ! -f "$copeimg" || ! -f "$varcopeimg" ]]; then
    echo "WARN: missing cope/varcope for $zimg" >&2
    echo "      expected: $copeimg" >&2
    echo "                $varcopeimg" >&2
  fi

  process_one "$zimg" "$copeimg" "$varcopeimg" "$featdir" \
    "$sub" "$ses" "$run" "$task" "$acq" "$confounds" "$znum" "$label" \
    >> "$out"
done

echo "[$(date '+%F %T')] Done. Wrote: $out"
