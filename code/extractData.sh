#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ------------------------------------------------------------
# Fixed locations relative to this script
# Assumes: /.../night-owls/code/extractData.sh
# ------------------------------------------------------------
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rootdir="$(dirname "$scriptdir")"

FSL_DERIV="${rootdir}/derivatives/fsl"
maskdir="${rootdir}/masks"
OUT_DIR="${rootdir}/derivatives/extractions"

# ------------------------------------------------------------
# MNI masks/maps (as provided)
# ------------------------------------------------------------
NAcc_MNI="$maskdir/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
CORT_MNI="$maskdir/BRS_Cortical_3pt1.nii.gz"
BRS_MNI="$maskdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

# ------------------------------------------------------------
# Required tools
# ------------------------------------------------------------
command -v fslstats >/dev/null || { echo "ERROR: fslstats not found in PATH" >&2; exit 1; }
command -v fslcc    >/dev/null || { echo "ERROR: fslcc not found in PATH" >&2; exit 1; }

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------
[[ -d "$FSL_DERIV" ]] || { echo "ERROR: Can't find FSL derivatives dir: $FSL_DERIV" >&2; exit 1; }
[[ -f "$NAcc_MNI"  ]] || { echo "ERROR: Missing NAcc MNI mask: $NAcc_MNI" >&2; exit 1; }
[[ -f "$CORT_MNI"  ]] || { echo "ERROR: Missing cortical mask: $CORT_MNI" >&2; exit 1; }
[[ -f "$BRS_MNI"   ]] || { echo "ERROR: Missing BRS MNI map: $BRS_MNI" >&2; exit 1; }

mkdir -p "$OUT_DIR"

# ------------------------------------------------------------
# Contrast label maps
# (Keep these aligned with your established numbering)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Find pattern: ONLY unsmoothed L1 FEAT zstats
# ------------------------------------------------------------
ZPATH="*/L1_*_unsmoothed.feat/stats/zstat*.nii.gz"

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------
out="${OUT_DIR}/extractions_L1stats_unsmoothed.tsv"
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tNAcc_zstat_mean\tNAcc_cope_mean\tNAcc_varcope_mean\tBRS_Cort_zstat_mean\tBRS_Cort_cope_mean\tBRS_Cort_varcope_mean\tBRS_corr" > "$out"

# ------------------------------------------------------------
# Count matches for a quick early warning
# ------------------------------------------------------------
total=$(
  find "$FSL_DERIV" -type f -path "$ZPATH" -print0 \
    | tr -cd '\0' | wc -c
)

if [[ "${total//[[:space:]]/}" == "0" ]]; then
  echo "WARN: No files matched:" >&2
  echo "      FSL_DERIV = $FSL_DERIV" >&2
  echo "      ZPATH     = $ZPATH" >&2
  echo "      Example expected FEAT dir like:" >&2
  echo "      .../sub-105/ses-12/L1_..._unsmoothed.feat/stats/zstat*.nii.gz" >&2
fi

done_cnt=0

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------
while IFS= read -r -d '' zimg; do
  ((done_cnt++))

  statsdir="$(dirname "$zimg")"
  featdir="$(dirname "$statsdir")"   # .../L1_..._unsmoothed.feat

  # Extract subject/session from path
  # Robust to deep paths as long as sub-*/ses-* exist somewhere above
  subdir="$(awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^sub-/){print $i; exit}}' <<<"$zimg")"
  sesdir="$(awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^ses-/){print $i; exit}}' <<<"$zimg")"

  sub="${subdir#sub-}"
  ses="${sesdir#ses-}"

  # FEAT basename without .feat
  fbase="$(basename "$featdir")"
  fbase="${fbase%.feat}"

  # Parse tags from FEAT name
  task="$(sed -E 's/^.*_task-([^_]+).*$/\1/' <<<"$fbase" 2>/dev/null || echo "")"
  run="$(sed -E 's/^.*_run-([0-9]+).*$/\1/' <<<"$fbase" 2>/dev/null || echo "NA")"
  space_raw="$(sed -E 's/^.*_space-([^_]+).*$/\1/' <<<"$fbase" 2>/dev/null || echo "")"
  confounds="$(sed -E 's/^.*_cnfds-([^_]+).*$/\1/' <<<"$fbase" 2>/dev/null || echo "unknown")"

  # Normalize/derive acquisition tag
  acq="unknown"
  [[ "$fbase" == *"_multi-echo_"*  ]] && acq="multiecho"
  [[ "$fbase" == *"_single-echo_"* ]] && acq="single"

  # Normalize space output tag
  space="t1w"
  [[ "${space_raw,,}" == mni* ]] && space="mni"

  # Only MNI
  [[ "$space" == "mni" ]] || continue

  # zstat index
  znum="$(basename "$zimg" | sed -E 's/^zstat([0-9]+).*$/\1/')"

  # Keep only requested contrasts
  label="$(contrast_label "$task" "$znum")"
  [[ -n "$label" ]] || continue

  copeimg="${statsdir}/cope${znum}.nii.gz"
  varcopeimg="${statsdir}/varcope${znum}.nii.gz"

  brainmask="${featdir}/mask.nii.gz"
  if [[ ! -f "$brainmask" ]]; then
    echo "WARN: missing FEAT mask (skipping): $brainmask" >&2
    continue
  fi

  # ROI means (default NA if missing)
  NAcc_z="NA"; NAcc_c="NA"; NAcc_v="NA"
  Cort_z="NA"; Cort_c="NA"; Cort_v="NA"

  [[ -f "$zimg"       ]] && NAcc_z="$(fslstats "$zimg"       -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$copeimg"    ]] && NAcc_c="$(fslstats "$copeimg"    -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$varcopeimg" ]] && NAcc_v="$(fslstats "$varcopeimg" -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"

  [[ -f "$zimg"       ]] && Cort_z="$(fslstats "$zimg"       -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$copeimg"    ]] && Cort_c="$(fslstats "$copeimg"    -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  [[ -f "$varcopeimg" ]] && Cort_v="$(fslstats "$varcopeimg" -k "$CORT_MNI" -M 2>/dev/null || echo NA)"

  if [[ ! -f "$copeimg" || ! -f "$varcopeimg" ]]; then
    echo "WARN: missing cope/varcope for zstat${znum} in $statsdir" >&2
  fi

  # Signed whole-brain spatial correlation with BRS using FEAT brain mask
  brs_corr="NA"
  if [[ -f "$zimg" ]]; then
    brs_corr="$(fslcc -m "$brainmask" --noabs -t -1 -p 6 "$zimg" "$BRS_MNI" 2>/dev/null | awk '{print $NF}' || echo NA)"
  fi

  echo -e "${sub}\t${ses}\t${run}\t${task}\t${space}\t${acq}\t${confounds}\t${znum}\t${label}\t${NAcc_z}\t${NAcc_c}\t${NAcc_v}\t${Cort_z}\t${Cort_c}\t${Cort_v}\t${brs_corr}" >> "$out"

  # Light-touch progress
  if (( done_cnt % 200 == 0 )); then
    echo "[$(date '+%F %T')] matched ${done_cnt} zstats so far..." >&2
  fi

done < <(find "$FSL_DERIV" -type f -path "$ZPATH" -print0)

echo "[$(date '+%F %T')] Done. Wrote: $out" >&2
