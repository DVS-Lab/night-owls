#!/usr/bin/env bash
# Extract ACF-based smoothness (effective FWHM) from fMRIPrep BOLD series.
# Handles both multiecho aggregate and single-echo (echo-2).
# Writes: sub  ses  task  run  acq  kernel_mm  fwhm_eff
#
# Usage:
#   FSL_ROOT=/path/to/derivatives/fsl \
#   FMRIPREP_ROOT=/path/to/derivatives/fmriprep \
#   TSV=./smoothness_acf.tsv \
#   KERNELS="0 5" \
#   VERBOSE=0 \
#   bash extractSmoothness.sh
#
# Or as args:
#   bash extractSmoothness.sh /path/to/derivatives/fsl /path/to/derivatives/fmriprep ./smoothness_acf.tsv

set -euo pipefail

# Guard: must be bash
if [ -z "${BASH_VERSION:-}" ]; then
  echo "ERROR: Run with bash (not sh)." >&2
  exit 2
fi

# ----------------------------
# Config / inputs
# ----------------------------
FSL_ROOT="${FSL_ROOT:-${1:-}}"
FMRIPREP_ROOT="${FMRIPREP_ROOT:-${2:-}}"
TSV="${TSV:-${3:-smoothness_acf.tsv}}"

# Kernels to extract (0 = raw; e.g., 5 = smoothed 5mm)
KERNELS="${KERNELS:-0 5}"

# Verbose progress (0 = quiet, 1 = per-row OK)
VERBOSE="${VERBOSE:-0}"

# If roots not set, try ../derivatives/*
if [ -z "${FSL_ROOT}" ] || [ -z "${FMRIPREP_ROOT}" ]; then
  PARENT="$(cd .. 2>/dev/null && pwd || pwd)"
  [ -z "${FSL_ROOT}" ] && [ -d "${PARENT}/derivatives/fsl" ] && FSL_ROOT="${PARENT}/derivatives/fsl"
  [ -z "${FMRIPREP_ROOT}" ] && [ -d "${PARENT}/derivatives/fmriprep" ] && FMRIPREP_ROOT="${PARENT}/derivatives/fmriprep"
fi

# Validate roots
if [ ! -d "${FSL_ROOT}" ]; then
  echo "ERROR: FSL_ROOT does not exist: ${FSL_ROOT:-<unset>}" >&2
  exit 1
fi
if [ ! -d "${FMRIPREP_ROOT}" ]; then
  echo "ERROR: FMRIPREP_ROOT does not exist: ${FMRIPREP_ROOT:-<unset>}" >&2
  exit 1
fi

# Check AFNI tool
if ! command -v 3dFWHMx >/dev/null 2>&1; then
  echo "ERROR: AFNI's 3dFWHMx not found on PATH." >&2
  exit 1
fi

echo "[ $(date +'%F %T') ] Scanning L1/MNI FEATs under: ${FSL_ROOT}"
echo "[ $(date +'%F %T') ] fMRIPrep root: ${FMRIPREP_ROOT}"
echo "[ $(date +'%F %T') ] Writing to: ${TSV}"

# ----------------------------
# FEAT discovery: L1 + MNI only
# ----------------------------
FEATS_FILE="$(mktemp)"
SEEN_FILE="$(mktemp)"   # dedupe on image path
trap 'rm -f "${FEATS_FILE}" "${SEEN_FILE}"' EXIT

# Keep only L1_* that are MNI; exclude L2, gfeat, T1w, subject-level
find "${FSL_ROOT}" -type d -name "L1_*" \
  -path "*/space-*mni*_*" \
  -not -path "*/L2_*" \
  -not -path "*/gfeat/*" \
  -not -path "*/space-*t1w*/*" \
  -not -path "*/subject-level/*" \
  2>/dev/null \
| sort -u > "${FEATS_FILE}"

# ----------------------------
# Helpers
# ----------------------------

# Parse identifiers out of FEAT dir name
parse_meta_from_feat() {
  local featdir="$1" base sub ses task run acq
  base="$(basename "$featdir")"
  sub="$(echo "$base" | sed -n 's/.*sub-\([0-9][0-9]*\).*/\1/p')"
  ses="$(echo "$base" | sed -n 's/.*ses-\([0-9][0-9]*\).*/\1/p')"
  task="$(echo "$base" | sed -n 's/.*task-\([^_]*\).*/\1/p')"
  run="$(echo "$base" | sed -n 's/.*run-\([0-9][0-9]*\).*/\1/p')"

  case "$base" in
    *single-echo*) acq="single-echo" ;;  # maps to echo-2
    *multi-echo*)  acq="multiecho"   ;;
    *)             acq=""            ;;
  esac

  echo "${sub} ${ses} ${task} ${run} ${acq}"
}

# Return fMRIPrep BOLD path for given IDs and kernel (0=raw, k>0 smoothed)
fmriprep_func_file() {
  local sub="$1" ses="$2" task="$3" run="$4" acq="$5" kernel="$6"
  local space="MNI152NLin6Asym"
  local funcdir="${FMRIPREP_ROOT}/sub-${sub}/ses-${ses}/func"
  local tail="_space-${space}_desc-preproc_bold"
  if [ "$kernel" -gt 0 ]; then
    tail="${tail}_${kernel}mm"
  fi
  tail="${tail}.nii.gz"

  if [ "$acq" = "single-echo" ]; then
    # explicit echo-2
    for f in "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_"*echo-2*"${tail}"; do
      [ -f "$f" ] && { echo "$f"; return 0; }
    done
  else
    # multiecho aggregate: any that DO NOT contain "_echo-"
    shopt -s nullglob
    for f in "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_"*"${tail}"; do
      case "$f" in
        *_echo-*) : ;;    # skip echo-specific files
        *) [ -f "$f" ] && { echo "$f"; shopt -u nullglob; return 0; } ;;
      esac
    done
    shopt -u nullglob
  fi

  echo ""   # not found
  return 1
}

# Compute effective FWHM (fourth field on the 2nd numeric line from -acf)
acf_fwhm_eff() {
  local mask="$1" img="$2"
  3dFWHMx -acf NULL -detrend -mask "${mask}" -input "${img}" 2>/dev/null \
    | awk 'NR==2 {print $4}'
}

# ----------------------------
# Output header
# ----------------------------
if [ ! -s "${TSV}" ]; then
  printf "sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff\n" > "${TSV}"
fi

rows=0

# ----------------------------
# Main
# ----------------------------
while IFS= read -r feat; do
  [ -z "${feat}" ] && continue

  # Parse identifiers
  IFS=' ' read -r sub ses task run acq <<< "$(parse_meta_from_feat "$feat")"
  # Basic validity
  if [ -z "${sub}" ] || [ -z "${ses}" ] || [ -z "${task}" ] || [ -z "${run}" ] || [ -z "${acq}" ]; then
    continue
  fi

  mask="${feat}/mask.nii.gz"
  [ ! -f "${mask}" ] && continue

  for k in ${KERNELS}; do
    img="$(fmriprep_func_file "$sub" "$ses" "$task" "$run" "$acq" "$k")"
    if [ -z "$img" ] || [ ! -f "$img" ]; then
      continue
    fi

    # De-dupe: skip if we've already processed this exact image
    if grep -qxF "$img" "${SEEN_FILE}"; then
      continue
    fi

    fwhm="$(acf_fwhm_eff "$mask" "$img")"
    [ -z "$fwhm" ] && continue

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sub" "$ses" "$task" "$run" "$acq" "$k" "$fwhm" >> "${TSV}"

    echo "$img" >> "${SEEN_FILE}"
    rows=$((rows + 1))
    if [ "${VERBOSE}" = "1" ]; then
      echo "OK: sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} k=${k} -> ${fwhm}"
    fi
  done
done < "${FEATS_FILE}"

echo "[ $(date +'%F %T') ] Done. Wrote ${rows} rows to ${TSV}"
