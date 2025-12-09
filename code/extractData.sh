#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Extract ROI means and BRS spatial correlations from subject-level
# L1 *unsmoothed* FEAT outputs.

# Resolve project root relative to this script
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rootdir="$(cd "${scriptdir}/.." && pwd)"

# Fixed project paths
FSL_DERIV="${rootdir}/derivatives/fsl"
FMRIPREP_DERIV="${rootdir}/derivatives/fmriprep"
maskdir="${rootdir}/masks"
outdir="${rootdir}/derivatives/extractions"

mkdir -p "$outdir"

# MNI masks/maps (provided)
NAcc_MNI="$maskdir/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
CORT_MNI="$maskdir/BRS_Cortical_3pt1.nii.gz"
BRS_MNI="$maskdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

# Output
outtsv="${outdir}/extractions_L1stats_unsmoothed.tsv"

die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARN:  $*" >&2; }

require_tool() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required tool: $1"
}

require_tool antsApplyTransforms
require_tool fslstats
require_tool fslcc
require_tool fslmaths
require_tool grep
require_tool awk
require_tool find
require_tool sort

[[ -d "$FSL_DERIV" ]] || die "Missing FSL derivatives directory: $FSL_DERIV"
[[ -d "$FMRIPREP_DERIV" ]] || die "Missing fMRIPrep derivatives directory: $FMRIPREP_DERIV"

[[ -f "$NAcc_MNI" ]] || die "NAcc mask not found: $NAcc_MNI"
[[ -f "$CORT_MNI" ]] || die "Cortical mask not found: $CORT_MNI"
[[ -f "$BRS_MNI"  ]] || die "BRS map not found: $BRS_MNI"

# Map zstat index to a readable contrast label for known tasks
contrast_label() {
  local task="$1" znum="$2"
  case "$task" in
    mid)
      case "$znum" in
        1) echo "gain>baseline" ;;
        2) echo "loss>baseline" ;;
        3) echo "gain>loss" ;;
        *) echo "z${znum}" ;;
      esac
      ;;
    sharedreward)
      case "$znum" in
        1) echo "friend>computer" ;;
        2) echo "stranger>computer" ;;
        3) echo "friend>stranger" ;;
        *) echo "z${znum}" ;;
      esac
      ;;
    *)
      echo "z${znum}"
      ;;
  esac
}

# Resample/transform an MNI image into the target grid.
# space_hint: "mni" (stats already in MNI) or "t1w"
# interp: defaults to Linear
mni_to_target() {
  local sub="$1" ses="$2" target_img="$3" src_mni="$4" out_img="$5"
  local space_hint="$6"
  local interp="${7:-Linear}"

  if [[ "$space_hint" == "mni" ]]; then
    antsApplyTransforms -d 3 -i "$src_mni" -r "$target_img" -n "$interp" -o "$out_img" >/dev/null 2>&1
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
    [[ -z "$h5" ]] && { warn "no MNI→T1w transform for sub-${sub} ses-${ses}; skipping"; return 1; }
    antsApplyTransforms -d 3 -i "$src_mni" -r "$target_img" -t "$h5" -n "$interp" -o "$out_img" >/dev/null 2>&1
  fi
}

mean_in_mask_or_na() {
  local img="$1" mask="$2"
  if [[ -f "$img" ]]; then
    fslstats "$img" -k "$mask" -M 2>/dev/null || echo "NA"
  else
    echo "NA"
  fi
}

brs_corr_or_na() {
  local img="$1" brs_img="$2" mask="$3"
  if [[ -f "$img" && -f "$brs_img" && -f "$mask" ]]; then
    fslcc -m "$mask" --noabs -t -1 -p 6 "$img" "$brs_img" 2>/dev/null | awk '{print $NF}' || echo "NA"
  else
    echo "NA"
  fi
}

# -------- per-image processing --------
process_one() {
  local zimg="$1" featdir="$2" sub="$3" ses="$4" run="$5" task="$6"
  local space_tag="$7" acq="$8" confounds="$9" znum="${10}" label="${11}"

  local brainmask="${featdir}/mask.nii.gz"
  [[ -f "$brainmask" ]] || { warn "missing FEAT mask: $brainmask — skipping"; return 0; }

  local tmpbase="${TMPROOT}/sub-${sub}_ses-${ses}_run-${run}_task-${task}_z${znum}"
  local nacc_res="${tmpbase}_NAcc_res.nii.gz"
  local cort_res="${tmpbase}_CORT_res.nii.gz"
  local brs_res="${tmpbase}_BRS_res.nii.gz"
  local brs_mask="${tmpbase}_BRSmask_res.nii.gz"

  # Bring masks/maps into the zstat grid
  mni_to_target "$sub" "$ses" "$zimg" "$NAcc_MNI" "$nacc_res" "$space_tag" "NearestNeighbor" || return 0
  mni_to_target "$sub" "$ses" "$zimg" "$CORT_MNI" "$cort_res" "$space_tag" "NearestNeighbor" || return 0
  mni_to_target "$sub" "$ses" "$zimg" "$BRS_MNI"  "$brs_res"  "$space_tag" "Linear"          || return 0

  # Intersection mask for BRS correlation: FEAT brainmask × BRS cortical mask
  fslmaths "$brainmask" -mul "$cort_res" "$brs_mask" >/dev/null 2>&1 || true

  local statsdir
  statsdir="$(dirname "$zimg")"

  local cope_img="${statsdir}/cope${znum}.nii.gz"
  local varcope_img="${statsdir}/varcope${znum}.nii.gz"

  local nacc_z_mean
  local nacc_cope_mean
  local nacc_varcope_mean
  local brs_corr

  nacc_z_mean="$(mean_in_mask_or_na "$zimg" "$nacc_res")"
  nacc_cope_mean="$(mean_in_mask_or_na "$cope_img" "$nacc_res")"
  nacc_varcope_mean="$(mean_in_mask_or_na "$varcope_img" "$nacc_res")"
  brs_corr="$(brs_corr_or_na "$zimg" "$brs_res" "$brs_mask")"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "sub-${sub}" "ses-${ses}" "run-${run}" "$task" "$acq" "${confounds:-NA}" "$space_tag" \
    "$znum" "$label" "$nacc_z_mean" "$nacc_cope_mean" "$nacc_varcope_mean" "$brs_corr" "$zimg"
}

# Temporary workspace
TMPROOT="$(mktemp -d -t extractData.XXXXXX)"
cleanup() { rm -rf "$TMPROOT"; }
trap cleanup EXIT

# Header
printf "sub\tses\trun\ttask\tacq\tconfounds\tspace\tznum\tcontrast\tNAcc_z_mean\tNAcc_cope_mean\tNAcc_varcope_mean\tBRS_corr\tzstat_path\n" > "$outtsv"

# Count expected files for progress reporting
total=$(
  find "$FSL_DERIV" -type f \
    -path "*/L1_*unsmoothed.feat/stats/zstat*.nii.gz" \
    ! -path "*/gfeat/*" -print0 2>/dev/null | tr -cd '\0' | wc -c
)
processed=0

# Main loop
while IFS= read -r -d '' zimg; do
  ((processed++)) || true

  featdir="$(cd "$(dirname "$zimg")/.." && pwd)"
  featbase="$(basename "$featdir")"

  # Parse metadata from FEAT directory name
  sub_raw="$(echo "$featbase" | grep -oE 'sub-[0-9]+' || true)"
  ses_raw="$(echo "$featbase" | grep -oE 'ses-[0-9]+' || true)"
  run_raw="$(echo "$featbase" | grep -oE 'run-[0-9]+' || true)"
  task_raw="$(echo "$featbase" | grep -oE 'task-[^_]+' || true)"
  space_raw="$(echo "$featbase" | grep -oE 'space-[^_]+' || true)"
  acq_raw="$(echo "$featbase" | grep -oE '(multi-echo|single-echo)' || true)"
  confounds_raw="$(echo "$featbase" | grep -oE 'cnfds-[^_]+' || true)"

  sub="${sub_raw#sub-}"
  ses="${ses_raw#ses-}"
  run="${run_raw#run-}"
  task="${task_raw#task-}"
  space="${space_raw#space-}"
  acq="${acq_raw:-NA}"
  confounds="${confounds_raw#cnfds-}"

  [[ -n "$sub" && -n "$ses" && -n "$task" ]] || { warn "Missing metadata in FEAT name: $featbase — skipping"; continue; }

  # Determine space tag
  space_tag="t1w"
  [[ "$space" == mni* || "$space" == MNI* ]] && space_tag="mni"
  [[ "$space" == T1w* || "$space" == t1w* ]] && space_tag="t1w"

  # Determine zstat index
  zbase="$(basename "$zimg")"
  znum="$(echo "$zbase" | grep -oE '[0-9]+' | head -n 1 || true)"
  [[ -n "$znum" ]] || { warn "Could not parse znum from $zimg — skipping"; continue; }

  label="$(contrast_label "$task" "$znum")"

  process_one "$zimg" "$featdir" "$sub" "$ses" "${run:-NA}" "$task" "$space_tag" \
    "$acq" "${confounds:-NA}" "$znum" "$label" >> "$outtsv"

  if (( processed % 250 == 0 )); then
    if (( total > 0 )); then
      echo "Processed ${processed}/${total} zstats..." >&2
    else
      echo "Processed ${processed} zstats..." >&2
    fi
  fi
done < <(
  find "$FSL_DERIV" -type f \
    -path "*/L1_*unsmoothed.feat/stats/zstat*.nii.gz" \
    ! -path "*/gfeat/*" -print0 2>/dev/null | sort -z
)

echo "Wrote: $outtsv" >&2
