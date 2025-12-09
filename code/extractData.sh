#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# -------------------- helpers --------------------
die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARN:  $*" >&2; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found in PATH: $1"
}

# -------------------- fixed locations (relative to THIS script) --------------------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rootdir="$(dirname "$scriptdir")"                 # project root (…/night-owls)

FSL_DERIV="${rootdir}/derivatives/fsl"
maskdir="${rootdir}/masks"
OUT_DIR="${rootdir}/derivatives/extractions"

# -------------------- your exact MNI masks/maps --------------------
NAcc_MNI="$maskdir/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
CORT_MNI="$maskdir/BRS_Cortical_3pt1.nii.gz"
BRS_MNI="$maskdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

# -------------------- tool checks --------------------
need_cmd fslstats
need_cmd fslcc
need_cmd find
need_cmd awk
need_cmd sed
need_cmd date
need_cmd wc
need_cmd tr

# -------------------- sanity checks --------------------
[[ -d "$FSL_DERIV" ]] || die "Can't find FSL derivatives directory: $FSL_DERIV"
[[ -d "$maskdir"   ]] || die "Can't find masks directory: $maskdir"

[[ -f "$NAcc_MNI" ]] || die "NAcc MNI mask not found: $NAcc_MNI"
[[ -f "$CORT_MNI" ]] || die "Cortical MNI mask not found: $CORT_MNI"
[[ -f "$BRS_MNI"  ]] || die "BRS MNI map not found: $BRS_MNI"

mkdir -p "$OUT_DIR"

# -------------------- contrast label maps --------------------
# Keeps your task-aware labels, but NEVER filters rows to empty.
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

# -------------------- tag parsing --------------------
extract_tag() {
  local key="$1" base="$2"
  if [[ "$base" =~ _${key}-([^_]+) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "NA"
  fi
}

# -------------------- per-contrast extraction --------------------
process_one() {
  local zimg="$1"

  local statsdir featdir sesdir subdir sub ses fbase
  statsdir="$(dirname "$zimg")"
  featdir="$(dirname "$statsdir")"
  sesdir="$(basename "$(dirname "$featdir")")"
  subdir="$(basename "$(dirname "$(dirname "$featdir")")")"

  sub="${subdir#sub-}"
  ses="${sesdir#ses-}"

  fbase="$(basename "$featdir")"
  fbase="${fbase%.feat}"

  local task run space_raw acq_raw confounds znum
  task="$(extract_tag "task" "$fbase")"
  run="$(extract_tag "run" "$fbase")"
  space_raw="$(extract_tag "space" "$fbase")"

  if [[ "$fbase" =~ _(multi-echo|single-echo)_ ]]; then
    acq_raw="${BASH_REMATCH[1]}"
  else
    acq_raw="NA"
  fi

  confounds="$(extract_tag "cnfds" "$fbase")"
  znum="$(basename "$zimg" | sed -E 's/^zstat([0-9]+).*$/\1/')"

  # normalize tags
  local space_tag="t1w"
  [[ "$space_raw" =~ ^(mni|MNI) ]] && space_tag="mni"

  local acq="NA"
  [[ "$acq_raw" == "single-echo" ]] && acq="single"
  [[ "$acq_raw" == "multi-echo"  ]] && acq="multiecho"

  # If everything is truly MNI, we should not be seeing non-MNI here.
  if [[ "$space_tag" != "mni" ]]; then
    warn "Skipping non-MNI image (no resampling requested): $zimg"
    return 0
  fi

  local label
  label="$(contrast_label "$task" "$znum")"
  [[ -n "$label" ]] || label="z${znum}"

  # corresponding cope/varcope
  local copeimg="${statsdir}/cope${znum}.nii.gz"
  local varcopeimg="${statsdir}/varcope${znum}.nii.gz"

  [[ -f "$copeimg"    ]] || warn "Missing cope for $zimg (expected $copeimg)"
  [[ -f "$varcopeimg" ]] || warn "Missing varcope for $zimg (expected $varcopeimg)"

  # --- zstat metrics ---
  local nacc_z="NA" cort_z="NA" brs_corr="NA"
  nacc_z="$(fslstats "$zimg" -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
  cort_z="$(fslstats "$zimg" -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  brs_corr="$(fslcc -m "$CORT_MNI" --noabs -t -1 -p 6 "$zimg" "$BRS_MNI" 2>/dev/null | awk '{print $NF}' || echo NA)"

  # --- cope metrics ---
  local nacc_c="NA" cort_c="NA"
  if [[ -f "$copeimg" ]]; then
    nacc_c="$(fslstats "$copeimg" -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
    cort_c="$(fslstats "$copeimg" -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  fi

  # --- varcope metrics ---
  local nacc_v="NA" cort_v="NA"
  if [[ -f "$varcopeimg" ]]; then
    nacc_v="$(fslstats "$varcopeimg" -k "$NAcc_MNI" -M 2>/dev/null || echo NA)"
    cort_v="$(fslstats "$varcopeimg" -k "$CORT_MNI" -M 2>/dev/null || echo NA)"
  fi

  echo -e "${sub}\t${ses}\t${run}\t${task}\t${space_tag}\t${acq}\t${confounds}\t${znum}\t${label}\t${nacc_z}\t${cort_z}\t${brs_corr}\t${nacc_c}\t${cort_c}\t${nacc_v}\t${cort_v}"
}

# -------------------- main --------------------
out="${OUT_DIR}/extractions_L1_unsmoothed_z_cope_varcope.tsv"
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tNAcc_mean_z\tCORT_mean_z\tBRS_corr_z\tNAcc_mean_cope\tCORT_mean_cope\tNAcc_mean_varcope\tCORT_mean_varcope" > "$out"

# Precision path: your specific unsmoothed L1 FEATs
ZPATH="${FSL_DERIV}/sub-*/ses-*/L1_*_unsmoothed.feat/stats/zstat*.nii.gz"

# quick count
total=$(
  find "$FSL_DERIV" \
    -type d \( -name '*.gfeat' -o -name 'subject-level' \) -prune -o \
    -type f -path "$ZPATH" -print0 \
  | tr -cd '\0' | wc -c
)

if (( total == 0 )); then
  warn "No zstat files matched:"
  warn "  $ZPATH"
  warn "Example expected pattern:"
  warn "  ${FSL_DERIV}/sub-*/ses-*/L1_*_unsmoothed.feat/stats/zstat*.nii.gz"
  exit 0
fi

done_cnt=0

while IFS= read -r -d '' zimg; do
  ((done_cnt++))

  if (( done_cnt % 200 == 0 || done_cnt == total )); then
    printf '[%s] %d/%d (%.1f%%) processed\n' \
      "$(date '+%F %T')" "$done_cnt" "$total" \
      "$(awk "BEGIN{print 100*$done_cnt/$total}")" >&2
  fi

  process_one "$zimg" >> "$out"
done < <(
  find "$FSL_DERIV" \
    -type d \( -name '*.gfeat' -o -name 'subject-level' \) -prune -o \
    -type f -path "$ZPATH" -print0
)

echo "[$(date '+%F %T')] Done. Wrote: $out"
