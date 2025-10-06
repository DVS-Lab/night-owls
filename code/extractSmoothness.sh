#!/usr/bin/env bash
# Extract ACF-based smoothness (effective FWHM) from fMRIPrep BOLD series
# for both unsmoothed (0 mm) and smoothed (e.g., 5 mm) images, across
# multi-echo aggregate and single-echo (echo-2) acquisitions.
#
# Usage (env or args):
#   FSL_ROOT=/path/to/derivatives/fsl \
#   FMRIPREP_ROOT=/path/to/derivatives/fmriprep \
#   TSV=./smoothness_acf.tsv \
#   KERNELS="0 5" \
#   VERBOSE=0 \
#   bash extractSmoothness.sh
#
# Or provide roots and tsv as args:
#   bash extractSmoothness.sh /path/to/derivatives/fsl /path/to/derivatives/fmriprep ./smoothness_acf.tsv
#
# Output TSV columns:
#   sub  ses  task  run  acq  kernel_mm  fwhm_eff
#
# Notes:
# - Only L1, MNI-space FEATs are scanned (no L2, no T1w, no gfeat).
# - acq="multi-echo" selects files WITHOUT an "echo-" token.
# - acq="single-echo" selects files WITH "echo-2".
# - We read the accurate combined FWHM (field 4 of the second numeric line from 3dFWHMx -acf).
# - Duplicates are suppressed by (sub|ses|task|run|acq|kernel_mm).

set -euo pipefail

# ----------------------------
# Config / inputs
# ----------------------------
FSL_ROOT="${FSL_ROOT:-${1:-}}"
FMRIPREP_ROOT="${FMRIPREP_ROOT:-${2:-}}"
TSV="${TSV:-${3:-smoothness_acf.tsv}}"

# Kernels to extract. Default assumes 0 (raw) and 5 mm (smoothed).
KERNELS="${KERNELS:-0 5}"

# Verbose progress (0 = quiet, 1 = print per-row OK)
VERBOSE="${VERBOSE:-0}"

# If roots not set, try to infer ../derivatives/* from current dir
if [[ -z "${FSL_ROOT}" || -z "${FMRIPREP_ROOT}" ]]; then
  PARENT="$(cd .. 2>/dev/null && pwd || pwd)"
  [[ -z "${FSL_ROOT}" ]] && [[ -d "${PARENT}/derivatives/fsl" ]] && FSL_ROOT="${PARENT}/derivatives/fsl"
  [[ -z "${FMRIPREP_ROOT}" ]] && [[ -d "${PARENT}/derivatives/fmriprep" ]] && FMRIPREP_ROOT="${PARENT}/derivatives/fmriprep"
fi

# Validate roots
if [[ ! -d "${FSL_ROOT}" ]]; then
  echo "ERROR: FSL_ROOT does not exist: ${FSL_ROOT:-<unset>}" >&2
  exit 1
fi
if [[ ! -d "${FMRIPREP_ROOT}" ]]; then
  echo "ERROR: FMRIPREP_ROOT does not exist: ${FMRIPREP_ROOT:-<unset>}" >&2
  exit 1
fi

# Check AFNI tool
if ! command -v 3dFWHMx >/dev/null 2>&1; then
  echo "ERROR: AFNI's 3dFWHMx not found on PATH." >&2
  exit 1
fi

ts() { date +'%F %T'; }

# ----------------------------
# FEAT discovery: L1 + MNI only
# ----------------------------
# Examples we want:
#   .../L1_sub-104_ses-03_task-mid_model-1_type-act_run-1_space-mni_single-echo_cnfds-fmriprep.feat
# We exclude:
#   - L2_* and any gfeat
#   - space-*t1w*
#   - subject-level trees (often cause permission noise and not per-run FEAT)
mapfile -t FEATS < <(
  find "${FSL_ROOT}" -type d -name "L1_*" \
    -path "*/space-*mni*_*" \
    -not -path "*/L2_*" \
    -not -path "*/gfeat/*" \
    -not -path "*/space-*t1w*/*" \
    -not -path "*/subject-level/*" \
    2>/dev/null \
  | sort -u
)

echo "[${ts}] Scanning L1/MNI FEATs under: ${FSL_ROOT}"
echo "[${ts}] fMRIPrep root: ${FMRIPREP_ROOT}"
echo "[${ts}] Writing to: ${TSV}"

# ----------------------------
# Helpers
# ----------------------------

parse_meta_from_feat() {
  # Print: sub ses task run acq
  local featdir="$1" base sub ses task run acq
  base="$(basename "$featdir")"

  # Extract IDs from anywhere in the name
  sub="$(sed -n 's/.*sub-\([0-9][0-9]*\).*/\1/p' <<<"$base")"
  ses="$(sed -n 's/.*ses-\([0-9][0-9]*\).*/\1/p' <<<"$base")"
  task="$(sed -n 's/.*task-\([^_]*\).*/\1/p' <<<"$base")"
  run="$(sed -n 's/.*run-\([0-9][0-9]*\).*/\1/p' <<<"$base")"

  # Acquisition encoded in FEAT name
  if [[ "$base" == *"single-echo"* ]]; then
    acq="single-echo"
  elif [[ "$base" == *"multi-echo"* ]]; then
    acq="multi-echo"
  else
    acq=""
  fi

  echo "$sub" "$ses" "$task" "$run" "$acq"
}

# Return fMRIPrep BOLD file path for given IDs and kernel (0 = raw, k>0 smoothed)
fmriprep_func_file() {
  local sub="$1" ses="$2" task="$3" run="$4" acq="$5" kernel="$6"
  local space="MNI152NLin6Asym"
  local funcdir="${FMRIPREP_ROOT}/sub-${sub}/ses-${ses}/func"
  local tail="_space-${space}_desc-preproc_bold"
  [[ "$kernel" -gt 0 ]] && tail="${tail}_${kernel}mm"
  tail="${tail}.nii.gz"

  shopt -s nullglob
  if [[ "$acq" == "single-echo" ]]; then
    # Require echo-2
    for f in "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_*"echo-2"*${tail}"; do
      [[ -f "$f" ]] && { echo "$f"; shopt -u nullglob; return 0; }
    done
  else
    # multi-echo aggregate → pick a file without any "_echo-" segment
    for f in "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_*"${tail}"; do
      [[ "$f" == *"_echo-"* ]] && continue
      [[ -f "$f" ]] && { echo "$f"; shopt -u nullglob; return 0; }
    done
  fi
  shopt -u nullglob

  echo ""   # not found
  return 1
}

# Compute effective FWHM from ACF model (field 4 on 2nd numeric line)
acf_fwhm_eff() {
  local mask="$1" img="$2"
  # -acf NULL avoids writing 3dFWHMx.1D and .png; we parse stdout only.
  3dFWHMx -acf NULL -detrend -mask "${mask}" -input "${img}" 2>/dev/null \
    | awk 'NR==2 {print $4}'
}

# ----------------------------
# Output header
# ----------------------------
if [[ ! -s "${TSV}" ]]; then
  printf "sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff\n" > "${TSV}"
fi

# De-dup set
declare -A seen=()

# ----------------------------
# Main
# ----------------------------
rows=0
for feat in "${FEATS[@]}"; do
  # Parse identifiers
  read -r sub ses task run acq < <(parse_meta_from_feat "$feat")
  [[ -z "$sub" || -z "$ses" || -z "$task" || -z "$run" || -z "$acq" ]] && continue

  mask="${feat}/mask.nii.gz"
  [[ ! -f "$mask" ]] && continue

  # For each kernel in KERNELS (e.g., 0 and 5)
  for k in ${KERNELS}; do
    key="${sub}|${ses}|${task}|${run}|${acq}|${k}"
    [[ -n "${seen[$key]:-}" ]] && continue

    img="$(fmriprep_func_file "$sub" "$ses" "$task" "$run" "$acq" "$k")"
    [[ -z "$img" || ! -f "$img" ]] && continue

    fwhm="$(acf_fwhm_eff "$mask" "$img")"
    [[ -z "$fwhm" ]] && continue

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sub" "$ses" "$task" "$run" "$acq" "$k" "$fwhm" >> "${TSV}"

    seen[$key]=1
    ((rows++))
    if [[ "${VERBOSE}" == "1" ]]; then
      echo "OK: sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} k=${k} -> ${fwhm}"
    fi
  done
done

echo "[${ts}] Done. Wrote ${rows} rows to ${TSV}"
