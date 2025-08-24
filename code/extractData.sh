#!/usr/bin/env bash
# extract_mid_sharedreward.sh
# Run from: derivatives/fsl/code
# Requires: FSL (fslstats, fslcc), ANTs (antsApplyTransforms) for warps when zstats are in T1w space.

set -euo pipefail

# --- paths (all relative to this "code" directory) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
fslroot="$(dirname "$scriptdir")"                         # .../derivatives/fsl
derivdir="$(dirname "$fslroot")"                          # .../derivatives
projroot="$(dirname "$derivdir")"                         # .../night-owls
masksdir="${projroot}/masks"

VS_MNI="${masksdir}/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"
BRS_MNI="${masksdir}/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

outdir="${fslroot}/extracts"
tmpdir="${scriptdir}/tmp_warps"
mkdir -p "$outdir" "$tmpdir"

outtsv="${outdir}/extract_mid_sharedreward.tsv"

# --- helpers ---
log(){ echo "[$(date +'%F %T')] $*"; }

need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH" >&2; exit 1; }; }

# If zstat is T1w space, warp MNI-space mask/map into that zstat space using ANTs transform from fMRIPrep
# $1=in_img (MNI mask or map), $2=ref_img (the zstat), $3=out_img, $4=sub-XXX, $5=ses-YY, $6=interp (NearestNeighbor|Linear)
warp_to_zspace_if_needed() {
  local in_img="$1" ref="$2" out_img="$3" sub="$4" ses="$5" interp="$6"
  local space_tag
  if [[ $(basename "$ref") == *".nii.gz" ]]; then
    # heuristic: if FEAT path contained 'space-t1w' treat as T1w, else MNI
    if [[ "$ref" == *"/space-t1w_"* || "$ref" == *"/space-t1w."* || "$ref" == *"/space-t1w/"* ]]; then
      space_tag="T1w"
    else
      space_tag="MNI152NLin6Asym"
    fi
  fi

  if [[ "$space_tag" == "T1w" ]]; then
    local xfm="${derivdir}/fmriprep/${sub}/${ses}/anat/${sub}_${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"
    if [[ ! -f "$xfm" ]]; then
      echo "WARN: missing transform: $xfm ; skipping warp and returning source image" >&2
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
    # already MNI
    cp -f "$in_img" "$out_img"
  fi
}

# mean extraction within mask
mean_in_mask() {
  local img="$1" mask="$2"
  fslstats "$img" -k "$mask" -M
}

# spatial correlation with BRS map, using FEAT whole-brain mask, keep sign, full range
corr_with_map() {
  local img="$1" map="$2" mask="$3"
  # -t -1 disables thresholding so we keep the full [-1,1] range
  fslcc --noabs -t -1 -m "$mask" "$img" "$map" | awk 'NF>=2{print $2; exit}'
}

# optional: label map by task + z-index if you want human-readable names
contrast_label() {
  local task="$1" zidx="$2"
  # Leave conservative defaults; fill in precise names as desired.
  case "$task" in
    mid)
      case "$zidx" in
        7) echo "MID:zstat7" ;;   # TODO: replace with your exact label for MID z7
        *) echo "MID:zstat${zidx}" ;;
      esac
      ;;
    sharedreward|shared-reward|shared_reward)
      case "$zidx" in
        7) echo "SharedReward:zstat7" ;;  # TODO: replace with your exact label for SharedReward z7
        *) echo "SharedReward:zstat${zidx}" ;;
      esac
      ;;
    *)
      echo "${task}:zstat${zidx}"
      ;;
  esac
}

# --- header ---
if [[ ! -f "$outtsv" ]]; then
  printf "sub\tses\trun\ttask\tspace\tacq\tconfounds\tzstat\tlabel\tVS_mean\tBRS_corr\tfeatdir\n" > "$outtsv"
fi

need fslstats
need fslcc

# --- main loop: scan all FEAT zstat7s (adjust glob to include others if desired) ---
# From derivatives/fsl root
while IFS= read -r -d '' zfile; do
  # featdir, run context
  statsdir="$(dirname "$zfile")"                   # .../stats
  featdir="$(dirname "$statsdir")"                 # .../L1_task-... .feat
  featbase="$(basename "$featdir")"                # L1_task-... .feat
  featstem="${featbase%.feat}"
  sesdir="$(dirname "$featdir")"                   # .../ses-XX
  subdir="$(dirname "$sesdir")"                    # .../sub-XXX

  sub="$(basename "$subdir")"                      # sub-101
  ses="$(basename "$sesdir")"                      # ses-01

  # parse fields from featstem
  # examples: L1_task-mid_model-1_type-act_run-1_space-mni_multi-echo_cnfds-tedana
  task="$(sed -nE 's/.*task-([^_]+).*/\1/p' <<<"$featstem")"
  run="$(sed -nE 's/.*_run-([0-9]+).*/\1/p' <<<"$featstem")"
  space_tag="$(sed -nE 's/.*_space-([^_]+).*/\1/p' <<<"$featstem")"      # mni | t1w
  acq="$(sed -nE 's/.*_(multi-echo|single-echo).*/\1/p' <<<"$featstem")" # multi-echo | single-echo
  confounds="$(sed -nE 's/.*_cnfds-([^_]+).*/\1/p' <<<"$featstem")"      # fmriprep | tedana

  zidx="$(basename "$zfile" | sed -nE 's/^zstat([0-9]+).*/\1/p')"
  label="$(contrast_label "$task" "$zidx")"

  wb_mask="${featdir}/mask.nii.gz"
  if [[ ! -f "$wb_mask" ]]; then
    log "WARN: missing FEAT mask: $wb_mask ; skipping ${zfile}"
    continue
  fi

  # Prepare masks/maps in the zstat space (cache per (sub,ses,space) to avoid recomputation)
  cache_tag="${sub}_${ses}_${space_tag}"
  vs_cache="${tmpdir}/VS_${cache_tag}.nii.gz"
  brs_cache="${tmpdir}/BRS_${cache_tag}.nii.gz"

  if [[ ! -f "$vs_cache" ]]; then
    warp_to_zspace_if_needed "$VS_MNI" "$zfile" "$vs_cache" "$sub" "$ses" "NearestNeighbor"
  fi
  if [[ ! -f "$brs_cache" ]]; then
    warp_to_zspace_if_needed "$BRS_MNI" "$zfile" "$brs_cache" "$sub" "$ses" "Linear"
  fi

  # Extract metrics
  vs_mean="NA"
  brs_cc="NA"

  if [[ -f "$vs_cache" ]]; then
    vs_mean="$(mean_in_mask "$zfile" "$vs_cache" 2>/dev/null || echo NA)"
  fi
  if [[ -f "$brs_cache" ]]; then
    brs_cc="$(corr_with_map "$zfile" "$brs_cache" "$wb_mask" 2>/dev/null || echo NA)"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tzstat%s\t%s\t%s\t%s\t%s\n" \
    "$sub" "$ses" "$run" "$task" "$space_tag" "$acq" "$confounds" \
    "$zidx" "$label" "$vs_mean" "$brs_cc" "$featdir" >> "$outtsv"

done < <(cd "$fslroot" && find . -type f -path "./sub-*/ses-*/L1_task-*.feat/stats/zstat7.nii.gz" -print0)

log "Done. Wrote: $outtsv"
