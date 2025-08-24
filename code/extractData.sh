#!/usr/bin/env bash
# extract_mid_sharedreward.sh
# Run from: derivatives/fsl/code

set -euo pipefail

# ---- layout (all relative to this "code" directory) ----
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
fslroot="$(dirname "$scriptdir")"              # .../derivatives/fsl
derivdir="$(dirname "$fslroot")"               # .../derivatives
projroot="$(dirname "$derivdir}")"             # .../night-owls
masksdir="${projroot}/masks"

VS_MNI="${masksdir}/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"
BRS_MNI="${masksdir}/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

outdir="${fslroot}/extracts"
tmpdir="${scriptdir}/tmp_warps"
mkdir -p "$outdir" "$tmpdir"

outtsv="${outdir}/extract_mid_sharedreward.tsv"

log(){ echo "[$(date +'%F %T')] $*"; }
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH" >&2; exit 1; }; }

# ---- label map: zstat# -> human-readable label ----
declare -A LABELS

# MID (from your screenshot)
LABELS["mid:7"]="ant_rew>neu"
LABELS["mid:8"]="pos>neg"
LABELS["mid:9"]="rew:pos>neg"
LABELS["mid:10"]="neu:pos>neg"

# SharedReward (from your screenshot)
LABELS["sharedreward:9"]="stranger>comp"
LABELS["sharedreward:10"]="neu>pun"
LABELS["sharedreward:11"]="rew>pun"
LABELS["sharedreward:12"]="S rew>pun"
LABELS["sharedreward:13"]="C rew>pun"
LABELS["sharedreward:14"]="S-C rew>pun"
LABELS["sharedreward:15"]="rew>neu"

label_for() {
  local task="$1" zidx="$2"
  local k="${task}:${zidx}"
  if [[ -n "${LABELS[$k]:-}" ]]; then
    printf "%s" "${LABELS[$k]}"
  else
    printf "%s:zstat%s" "$task" "$zidx"
  fi
}

# ---- warping helpers ----
warp_to_zspace_if_needed() {
  local in_img="$1" ref="$2" out_img="$3" sub="$4" ses="$5" interp="$6"
  local tag="MNI"
  [[ "$ref" == *"/space-t1w"* ]] && tag="T1w"
  if [[ "$tag" == "T1w" ]]; then
    # session-specific fMRIPrep transform MNI->T1w
    local xfm="${derivdir}/fmriprep/${sub}/${ses}/anat/${sub}_${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"
    if [[ ! -f "$xfm" ]]; then
      echo "WARN: missing transform: $xfm ; copying input mask/map without warp" >&2
      cp -f "$in_img" "$out_img"
      return 0
    fi
    need antsApplyTransforms
    antsApplyTransforms -d 3 \
      -i "$in_img" \
      -r "$ref" \
      -o "$out_img" \
      -t "$xfm" \
      -n "$interp" \
      --float
  else
    cp -f "$in_img" "$out_img"
  fi
}

mean_in_mask() {
  local img="$1" mask="$2"
  fslstats "$img" -k "$mask" -M
}

corr_with_map() {
  local img="$1" map="$2" mask="$3"
  # keep sign; use full range by setting threshold to -1
  fslcc --noabs -t -1 -m "$mask" "$img" "$map" | awk 'NF>=2{print $2; exit}'
}

# ---- header ----
if [[ ! -f "$outtsv" ]]; then
  printf "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tVS_mean\tBRS_corr\tfeatdir\n" > "$outtsv"
fi

need fslstats
need fslcc

# ---- main scan ----
while IFS= read -r -d '' zfile; do
  statsdir="$(dirname "$zfile")"
  featdir="$(dirname "$statsdir")"
  featbase="$(basename "$featdir")"           # L1_task-... .feat
  featstem="${featbase%.feat}"
  sesdir="$(dirname "$featdir")"              # .../ses-XX
  subdir="$(dirname "$sesdir")"               # .../sub-XXX

  sub="$(basename "$subdir")"                 # sub-101
  ses="$(basename "$sesdir")"                 # ses-01

  task="$(sed -nE 's/.*task-([^_]+).*/\1/p' <<<"$featstem")"
  run="$(sed -nE 's/.*_run-([0-9]+).*/\1/p' <<<"$featstem")"
  space_tag="$(sed -nE 's/.*_space-([^_]+).*/\1/p' <<<"$featstem")"       # mni | t1w
  acq="$(sed -nE 's/.*_(multi-echo|single-echo).*/\1/p' <<<"$featstem")"  # multi-echo | single-echo
  confounds="$(sed -nE 's/.*_cnfds-([^_]+).*/\1/p' <<<"$featstem")"       # fmriprep | tedana

  zidx="$(basename "$zfile" | sed -nE 's/^zstat([0-9]+).*/\1/p')"
  label="$(label_for "$task" "$zidx")"

  wb_mask="${featdir}/mask.nii.gz"
  if [[ ! -f "$wb_mask" ]]; then
    log "WARN: missing FEAT mask: $wb_mask ; skipping ${zfile}"
    continue
  fi

  # cache per (sub,ses,space)
  cache_tag="${sub}_${ses}_${space_tag}"
  vs_cache="${tmpdir}/VS_${cache_tag}.nii.gz"
  brs_cache="${tmpdir}/BRS_${cache_tag}.nii.gz"

  [[ -f "$vs_cache" ]] || warp_to_zspace_if_needed "$VS_MNI" "$zfile" "$vs_cache" "$sub" "$ses" "NearestNeighbor"
  [[ -f "$brs_cache" ]] || warp_to_zspace_if_needed "$BRS_MNI" "$zfile" "$brs_cache" "$sub" "$ses" "Linear"

  vs_mean="NA"
  brs_cc="NA"
  [[ -f "$vs_cache" ]] && vs_mean="$(mean_in_mask "$zfile" "$vs_cache" 2>/dev/null || echo NA)"
  [[ -f "$brs_cache" ]] && brs_cc="$(corr_with_map "$zfile" "$brs_cache" "$wb_mask" 2>/dev/null || echo NA)"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tzstat%s\t%s\t%s\t%s\t%s\n" \
    "$sub" "$ses" "$run" "$task" "$space_tag" "$acq" "$confounds" \
    "$zidx" "$label" "$vs_mean" "$brs_cc" "$featdir" >> "$outtsv"

done < <(cd "$fslroot" && find . -type f -path "./sub-*/ses-*/L1_task-*.feat/stats/zstat*.nii.gz" -print0)

log "Done. Wrote: $outtsv"
