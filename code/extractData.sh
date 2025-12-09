#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# --------------------------------------------
# extractData.sh
#
# Extract ROI means from L1 UNSMOOTHED zstat/cope/varcope images
# and compute whole-brain BRS spatial correlation (zstat only),
# MNI-space only.
#
# Assumes:
#   project_root/
#     code/              (this script lives here)
#     masks/
#     derivatives/fsl/
#     derivatives/extractions/
# --------------------------------------------

# -------- fixed locations --------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"

FSL_DERIV="${rootdir}/derivatives/fsl"
maskdir="${rootdir}/masks"
outdir="${rootdir}/derivatives/extractions"

# -------- MNI masks/maps (as provided) --------
NAcc_MNI="$maskdir/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
CORT_MNI="$maskdir/BRS_Cortical_3pt1.nii.gz"
BRS_MNI="$maskdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

# -------- required tools --------
command -v fslstats >/dev/null || { echo "ERROR: fslstats not found in PATH" >&2; exit 1; }
command -v fslcc    >/dev/null || { echo "ERROR: fslcc not found in PATH" >&2; exit 1; }

# -------- sanity checks --------
[[ -d "$FSL_DERIV" ]] || { echo "ERROR: Can't find $FSL_DERIV" >&2; exit 1; }
[[ -f "$NAcc_MNI" ]]  || { echo "ERROR: Missing NAcc mask: $NAcc_MNI" >&2; exit 1; }
[[ -f "$CORT_MNI" ]]  || { echo "ERROR: Missing cortical mask: $CORT_MNI" >&2; exit 1; }
[[ -f "$BRS_MNI" ]]   || { echo "ERROR: Missing BRS map: $BRS_MNI" >&2; exit 1; }

mkdir -p "$outdir"

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

# -------- find pattern (tight to avoid noise) --------
ZPATH="*/sub-*/ses-*/L1_*_unsmoothed.feat/stats/zstat*.nii.gz"

find_zstats() {
  # Prune group/subject-level outputs and suppress permission chatter.
  find "$FSL_DERIV" \
    -type d \( -name '*.gfeat' -o -name 'subject-level' \) -prune -o \
    -type f -path "$ZPATH" -print0 2>/dev/null || true
}

# -------- per-image processing (MNI-only) --------
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

  # Signed whole-brain spatial corr with BRS (mask to FEAT brainmask), zstat only
  local brs_corr="NA"
  if [[ -f "$zimg" ]]; then
    brs_corr="$(fslcc -m "$brainmask" --noabs -t -1 -p 6 "$zimg" "$BRS_MNI" 2>/dev/null | awk '{print $NF}' || echo NA)"
  fi

  echo -e "${sub}\t${ses}\t${run}\t${task}\tmni\t${acq}\t${confounds}\t${znum}\t${label}\t${NAcc_z}\t${NAcc_c}\t${NAcc_v}\t${Cort_z}\t${Cort_c}\t${Cort_v}\t${brs_corr}"
}

# -------- main --------
out="${outdir}/extractions_L1stats_unsmoothed.tsv"
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tNAcc_zstat_mean\tNAcc_cope_mean\tNAcc_varcope_mean\tBRS_Cort_zstat_mean\tBRS_Cort_cope_mean\tBRS_Cort_varcope_mean\tBRS_corr" \
  > "$out"

total="$(find_zstats | tr -cd '\0' | wc -c | tr -d ' ')"
done_cnt=0

if [[ "${total:-0}" -eq 0 ]]; then
  echo "WARN: No zstat files found with pattern:" >&2
  echo "      $ZPATH" >&2
  echo "      under $FSL_DERIV" >&2
fi

while IFS= read -r -d '' zimg; do
  ((done_cnt++))

  if (( total > 0 )) && (( done_cnt % 200 == 0 || done_cnt == total )); then
    pct="$(awk "BEGIN{print 100*$done_cnt/$total}")"
    printf '[%s] %d/%d (%.1f%%) processed\n' "$(date '+%F %T')" "$done_cnt" "$total" "$pct" >&2
  fi

  # Derive FEAT dir and BIDS-ish identifiers
  featdir="$(dirname "$(dirname "$zimg")")"                 # .../L1_...feat
  sesdir="$(basename "$(dirname "$featdir")")"              # ses-XX
  subdir="$(basename "$(dirname "$(dirname "$featdir")")")" # sub-XXX

  sub="${subdir#sub-}"
  ses="${sesdir#ses-}"

  fbase="$(basename "$featdir")"
  fbase="${fbase%.feat}"

  # Parse from FEAT name
  task="$(sed -E 's/^.*_task-([^_]+).*$/\1/' <<<"$fbase")"
  run="$(sed -E 's/^.*_run-([0-9]+).*$/\1/' <<<"$fbase")"
  space_raw="$(sed -E 's/^.*_space-([^_]+).*$/\1/' <<<"$fbase")"
  acq_raw="$(sed -E 's/^.*_(multi-echo|single-echo)_.*$/\1/' <<<"$fbase")"
  confounds="$(sed -E 's/^.*_cnfds-([^_]+).*$/\1/' <<<"$fbase")"

  znum="$(basename "$zimg" | sed -E 's/^zstat([0-9]+).*$/\1/')"

  # Normalize tags
  space_tag="t1w"
  [[ "$space_raw" =~ ^(mni|MNI) ]] && space_tag="mni"

  acq="single"
  [[ "$acq_raw" == "multi-echo" ]] && acq="multiecho"

  # MNI-only
  [[ "$space_tag" == "mni" ]] || continue

  label="$(contrast_label "$task" "$znum")"
  # Keep all matching zstats even if label is empty
  # (label can be blank for unlisted contrasts)

  statsdir="$(dirname "$zimg")"
  copeimg="${statsdir}/cope${znum}.nii.gz"
  varcopeimg="${statsdir}/varcope${znum}.nii.gz"

  if [[ ! -f "$copeimg" || ! -f "$varcopeimg" ]]; then
    echo "WARN: missing cope/varcope for $zimg" >&2
    echo "      expected: $copeimg" >&2
    echo "                $varcopeimg" >&2
  fi

  process_one "$zimg" "$copeimg" "$varcopeimg" "$featdir" \
    "$sub" "$ses" "$run" "$task" "$acq" "$confounds" "$znum" "${label}" \
    >> "$out"

done < <(find_zstats)

echo "[$(date '+%F %T')] Done. Wrote: $out" >&2
