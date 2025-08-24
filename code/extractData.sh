#!/usr/bin/env bash
set -euo pipefail

umask 0000

# --- locate project relative to this script ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
projectdir="$(dirname "$scriptdir")"

DERIV_FSL="$projectdir/derivatives/fsl"
MASKS_DIR="$projectdir/masks"
OUTDIR="$projectdir/extracts"
mkdir -p "$OUTDIR"
OUTTSV="$OUTDIR/extract_mid_sharedreward.tsv"

# --- required tools ---
for cmd in fslstats fslcc fslmaths fslnvols ; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found in PATH"; exit 1; }
done

# --- masks in standard (MNI152NLin6Asym) space ---
VS_MASK_MNI="$MASKS_DIR/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"
BRS_MAP_MNI="$MASKS_DIR/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"
[[ -f "$VS_MASK_MNI" && -f "$BRS_MAP_MNI" ]] || { echo "ERROR: expected masks not found in $MASKS_DIR"; exit 1; }

# --- small helper: extract token value from FEAT basename like: L1_task-..._run-2_space-mni_multi-echo_cnfds-tedana.feat
get_token () {  # $1=featbase  $2=key (e.g., task, run, space)
  local b="$1" k="$2"
  echo "$b" | sed -n "s/.*_${k}-\([^_]*\)\(_\|$\).*/\1/p"
}

# optional mapping stub: default to "zstat<N>" if no nicer label known
label_for () { # $1=task  $2=znum
  local task="$1" z="$2"
  # If you want human-readable names, hard-code your map here.
  # For now we keep the zstat number so you can proceed immediately.
  echo "zstat${z}"
}

# --- fslcc with subject-wise brainmask; keep sign; no thresholding of input (t=-1 => include full range) ---
corr_with_brs () { # $1=zstat_img  $2=featdir  $3=map_in_same_space
  local zimg="$1" feat="$2" map="$3"
  local brainmask="$feat/mask.nii.gz"
  [[ -f "$brainmask" ]] || { echo "NA"; return; }
  # fslcc prints: "<corr>    <path1>    <path2>"
  local cc
  if ! cc=$(fslcc -m "$brainmask" --noabs -t -1 "$zimg" "$map" 2>/dev/null | awk '{print $1}'); then
    echo "NA"; return
  fi
  [[ -z "$cc" ]] && echo "NA" || echo "$cc"
}

# --- Mean within VS mask (assumes mask already in same space as zimg) ---
mean_in_mask () { # $1=zstat_img  $2=mask_same_space
  local zimg="$1" mask="$2"
  local val
  if ! val=$(fslstats "$zimg" -k "$mask" -M 2>/dev/null); then
    echo "NA"; return
  fi
  [[ -z "$val" ]] && echo "NA" || echo "$val"
}

# --- ensure masks are in the same space as the zstat ---
# If zstat is MNI, use the MNI masks directly.
# If zstat is T1w, reuse the FEAT's own mask as a crude brain mask; for VS/BRS we’ll align to T1w by copying MNI masks’ voxel grid via flirt (nearest) to the FEAT space if needed.
# (This avoids poking into fmriprep internals; we stay inside derivatives/.)
prep_masks_for_feat () { # $1=space (mni|t1w)  $2=zimg  $3=featdir  -> prints "VS_PATH  BRS_PATH"
  local space="$1" zimg="$2" feat="$3"
  if [[ "$space" == "mni" ]]; then
    echo "$VS_MASK_MNI $BRS_MAP_MNI"
  else
    # Bring MNI masks to the zstat grid using nearest-neighbor resampling
    local vs="$feat/aux_VS_inT1w.nii.gz"
    local brs="$feat/aux_BRS_inT1w.nii.gz"
    if [[ ! -f "$vs" ]]; then
      flirt -in "$VS_MASK_MNI"  -ref "$zimg" -out "$vs"  -applyxfm -usesqform -interp nearestneighbour >/dev/null 2>&1 || true
    fi
    if [[ ! -f "$brs" ]]; then
      flirt -in "$BRS_MAP_MNI"  -ref "$zimg" -out "$brs" -applyxfm -usesqform -interp nearestneighbour >/dev/null 2>&1 || true
    fi
    echo "$vs $brs"
  fi
}

# --- write header ---
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tVS_mean\tBRS_corr\tfeatdir" > "$OUTTSV"

# --- collect only L1_* FEAT zstats; prune L2_* and *.gfeat completely ---
mapfile -t ZLIST < <(
  find "$DERIV_FSL" \
    \( -path '*/L2_*' -o -path '*/*.gfeat' \) -prune -o \
    -type f -path '*/L1_*.feat/stats/zstat*.nii.gz' -print \
  | sort
)

# bail out early if nothing found (helps debugging wrong roots)
if [[ "${#ZLIST[@]}" -eq 0 ]]; then
  echo "No zstat images found under $DERIV_FSL (L1_* only). Check paths."
  exit 1
fi

# --- main loop ---
for zimg in "${ZLIST[@]}"; do
  feat="$(dirname "$(dirname "$zimg")")"                 # .../L1_....feat
  [[ "$(basename "$feat")" == L1_* ]] || continue        # belt-and-suspenders

  rel="${feat#"$DERIV_FSL/"}"                            # sub-XXX/ses-YY/L1_...feat
  sub="$(echo "$rel" | cut -d/ -f1)"                     # sub-101
  ses="$(echo "$rel" | cut -d/ -f2)"                     # ses-01
  featbase="$(basename "$feat" .feat)"                   # L1_task-..._...

  task="$(get_token "$featbase" task)"
  run="$(get_token "$featbase" run)"
  rawspace="$(get_token "$featbase" space)"               # mni or t1w (based on your names)
  space="$rawspace"                                       # keep as mni|t1w for transparency

  # acquisition: token is "multi-echo" or "single-echo" in your names
  acq_token="$(echo "$featbase" | sed -n 's/.*_\(multi-echo\|single-echo\)\(_\|$\).*/\1/p')"
  acq="${acq_token//-/}"                                  # multiecho|singleecho
  confounds="$(get_token "$featbase" cnfds)"              # fmriprep|tedana

  zname="$(basename "$zimg" .nii.gz)"                     # zstat7
  znum="${zname#zstat}"                                   # 7
  label="$(label_for "$task" "$znum")"

  # prep VS/BRS masks in this zstat’s space
  read VS_MASK BRS_MAP < <(prep_masks_for_feat "$space" "$zimg" "$feat")

  # compute metrics
  vs_mean="$(mean_in_mask "$zimg" "$VS_MASK")"
  brs_corr="$(corr_with_brs "$zimg" "$feat" "$BRS_MAP")"

  # append row
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sub" "$ses" "$run" "$task" "$space" "$acq" "$confounds" \
    "$zname" "$label" "$vs_mean" "$brs_corr" "$feat" >> "$OUTTSV"
done

echo "[${HOSTNAME:-node}] $(date +'%F %T') Done. Wrote: $OUTTSV"
