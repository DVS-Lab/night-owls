#!/usr/bin/env bash
set -euo pipefail
umask 0000

# Run from: night-owls/code
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
projectdir="$(dirname "$scriptdir")"

DERIV_FSL="$projectdir/derivatives/fsl"
MASKS_DIR="$projectdir/masks"
OUTDIR="$projectdir/derivatives/extractions"
mkdir -p "$OUTDIR"
OUTTSV="$OUTDIR/extractions_L1stats.tsv"

# Tools
for cmd in fslstats fslcc flirt ; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found in PATH"; exit 1; }
done

# Masks in MNI space
VS_MASK_MNI="$MASKS_DIR/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"
BRS_MAP_MNI="$MASKS_DIR/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"
[[ -f "$VS_MASK_MNI" && -f "$BRS_MAP_MNI" ]] || { echo "ERROR: expected masks not found in $MASKS_DIR"; exit 1; }

# --- contrast label maps (from your screenshots) ---
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

# Helpers
get_token () { echo "$1" | sed -n "s/.*_${2}-\([^_]*\)\(_\|$\).*/\1/p"; }

# Correlation with BRS using FEAT brainmask; keep sign; full range
corr_with_brs () { # $1=zimg  $2=featdir  $3=map
  local zimg="$1" feat="$2" map="$3" brainmask="$feat/mask.nii.gz"
  [[ -f "$brainmask" ]] || { echo "NA"; return; }
  fslcc -m "$brainmask" --noabs -t -1 "$zimg" "$map" 2>/dev/null | awk '{print $1}'
}

# Mean within mask
mean_in_mask () { fslstats "$1" -k "$2" -M 2>/dev/null || echo "NA"; }

# Prepare masks in the zstat grid for T1w space; for MNI just use originals
prep_masks_for_feat () { # $1=space(mni|t1w) $2=zimg $3=featdir
  local space="$1" zimg="$2" feat="$3"
  if [[ "$space" == "mni" ]]; then
    echo "$VS_MASK_MNI $BRS_MAP_MNI"
  else
    local vs="$feat/aux_VS_inT1w.nii.gz"
    local brs="$feat/aux_BRS_inT1w.nii.gz"
    [[ -f "$vs" ]] || flirt -in "$VS_MASK_MNI" -ref "$zimg" -out "$vs" -applyxfm -usesqform -interp nearestneighbour >/dev/null 2>&1
    [[ -f "$brs" ]] || flirt -in "$BRS_MAP_MNI" -ref "$zimg" -out "$brs" -applyxfm -usesqform -interp trilinear >/dev/null 2>&1
    echo "$vs $brs"
  fi
}

# Output header
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tVS_mean\tBRS_corr\tfeatdir" > "$OUTTSV"

# Find ONLY L1 FEAT zstats under derivatives/fsl (never scratch, never L2, never gfeat)
mapfile -t ZLIST < <(
  find "$DERIV_FSL" \
    \( -path '*/L2_*' -o -path '*/*.gfeat' \) -prune -o \
    -type f -path '*/L1_*.feat/stats/zstat*.nii.gz' -print | sort
)

[[ "${#ZLIST[@]}" -gt 0 ]] || { echo "No L1 zstats found under $DERIV_FSL"; exit 1; }

# Main
for zimg in "${ZLIST[@]}"; do
  feat="$(dirname "$(dirname "$zimg")")"
  [[ "$(basename "$feat")" == L1_* ]] || continue

  rel="${feat#"$DERIV_FSL/"}"            # sub-XXX/ses-YY/...
  sub="$(echo "$rel" | cut -d/ -f1)"     # sub-101
  ses="$(echo "$rel" | cut -d/ -f2)"     # ses-01
  featbase="$(basename "$feat" .feat)"

  task="$(get_token "$featbase" task)"                   # mid | sharedreward | ...
  [[ "$task" == "mid" || "$task" == "sharedreward" ]] || continue

  run="$(get_token "$featbase" run)"
  space="$(get_token "$featbase" space)"                 # mni | t1w
  acq_token="$(echo "$featbase" | sed -n 's/.*_\(multi-echo\|single-echo\)\(_\|$\).*/\1/p')"
  acq="${acq_token//-echo/}"; acq="${acq//-/}"           # multiecho | single
  confounds="$(get_token "$featbase" cnfds)"             # fmriprep | tedana

  zname="$(basename "$zimg" .nii.gz)"                    # zstat7
  znum="${zname#zstat}"
  label="$(contrast_label "$task" "$znum")"
  [[ -n "$label" ]] || continue  # only keep zstats we explicitly mapped

  read VS_MASK BRS_MAP < <(prep_masks_for_feat "$space" "$zimg" "$feat")
  vs_mean="$(mean_in_mask "$zimg" "$VS_MASK")"
  brs_corr="$(corr_with_brs "$zimg" "$feat" "$BRS_MAP")"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sub" "$ses" "$run" "$task" "$space" "$acq" "$confounds" \
    "$zname" "$label" "$vs_mean" "$brs_corr" "$feat" >> "$OUTTSV"
done

echo "[${HOSTNAME:-node}] $(date +'%F %T') Done. Wrote: $OUTTSV"
