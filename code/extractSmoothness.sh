#!/usr/bin/env bash
# extractSmoothness.sh
# Minimal, robust smoothness extractor for L1 FEATs (MNI space) from fMRIPrep outputs.
# - Outputs: sub  ses  task  run  acq  kernel_mm  fwhm_eff
# - acq values: "multiecho" (aggregate, no _echo-*) and "single-echo" (echo-2)
# - Kernels: 0 (unsmoothed) and 5 (smoothed) by default; override via KERNELS env var
# - Run with bash (not sh). Example:
#     FSL_ROOT="/.../derivatives/fsl" \
#     FMRIPREP_ROOT="/.../derivatives/fmriprep" \
#     OUT="smoothness_acf.tsv" \
#     bash extractSmoothness.sh

set -e

# -------- Config (override via environment) -----------------------------------
FSL_ROOT="${FSL_ROOT:-/home/you/path/night-owls/derivatives/fsl}"
FMRIPREP_ROOT="${FMRIPREP_ROOT:-/home/you/path/night-owls/derivatives/fmriprep}"
OUT="${OUT:-/home/you/path/night-owls/derivatives/extractions/smoothness_acf.tsv}"
# default kernels: unsmoothed (0) and smoothed (5). Add others if you used them.
KERNELS="${KERNELS:-0 5}"
# whitelist to avoid spurious tasks
TASK_WHITELIST="${TASK_WHITELIST:-mid sharedreward}"

# -------- Helpers -------------------------------------------------------------

# Build the filename tail for a given kernel (0 means unsmoothed)
make_tail() {
  # $1 = kernel mm (integer)
  local k="$1"
  local tail="_space-MNI152NLin6Asym_desc-preproc_bold"
  if [ "$k" -gt 0 ]; then
    tail="${tail}_${k}mm"
  fi
  echo "${tail}.nii.gz"
}

# Find the *aggregate multiecho* (no _echo-*) image for a run+kernel.
# Matches examples like:
#   sub-101_ses-01_task-mid_run-1_part-mag_space-MNI152NLin6Asym_desc-preproc_bold[_5mm].nii.gz
find_me_agg() {
  # $1=sub  $2=ses  $3=task  $4=run  $5=kernel
  local sub="$1" ses="$2" task="$3" run="$4" k="$5"
  local funcdir="${FMRIPREP_ROOT}/sub-${sub}/ses-${ses}/func"
  local tail
  tail="$(make_tail "$k")"

  # Enable nullglob so unmatched patterns vanish, not echo literally
  shopt -s nullglob
  # Prefer part-mag if present, but allow any part-* as a fallback
  for f in \
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag"*"$tail" \
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-"*"$tail" \
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}"*"$tail"
  do
    case "$f" in
      *_echo-*) ;;  # skip echo-specific files
      *)
        if [ -f "$f" ]; then
          echo "$f"
          shopt -u nullglob
          return 0
        fi
        ;;
    esac
  done
  shopt -u nullglob
  return 1
}

# Find the *single-echo echo-2* image for a run+kernel.
# Matches examples like:
#   sub-101_ses-01_task-mid_run-1_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold[_5mm].nii.gz
find_se_e2() {
  # $1=sub  $2=ses  $3=task  $4=run  $5=kernel
  local sub="$1" ses="$2" task="$3" run="$4" k="$5"
  local funcdir="${FMRIPREP_ROOT}/sub-${sub}/ses-${ses}/func"
  local tail
  tail="$(make_tail "$k")"

  shopt -s nullglob
  for f in \
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2_part-mag"*"$tail" \
    "${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_echo-2"*"$tail"
  do
    if [ -f "$f" ]; then
      echo "$f"
      shopt -u nullglob
      return 0
    fi
  done
  shopt -u nullglob
  return 1
}

# Extract "effective FWHM" from AFNI 3dFWHMx using accurate ACF modeling.
# We pass -acf NULL to silence side files; parse the 2nd line's 4th column.
get_fwhm_eff() {
  # $1 = mask NIfTI, $2 = BOLD NIfTI
  local mask="$1" img="$2"
  3dFWHMx -acf NULL -detrend -mask "$mask" -input "$img" 2>/dev/null | awk 'NR==2{print $4}'
}

# Write TSV header (once). If file exists but empty, add header.
ensure_header() {
  local tsv="$1"
  if [ ! -s "$tsv" ]; then
    printf "sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff\n" > "$tsv"
  fi
}

# Simple membership test for whitelist words
in_whitelist() {
  # $1=item  $2=list of words
  local needle="$1"
  shift
  for w in "$@"; do
    if [ "$needle" = "$w" ]; then
      return 0
    fi
  done
  return 1
}

# -------- Main ----------------------------------------------------------------

echo "[$(date +'%F %T')] Starting smoothness extraction"
echo "FSL_ROOT=${FSL_ROOT}"
echo "FMRIPREP_ROOT=${FMRIPREP_ROOT}"
echo "OUT=${OUT}"
echo "KERNELS=${KERNELS}"

# Ensure AFNI tool exists
if ! command -v 3dFWHMx >/dev/null 2>&1; then
  echo "ERROR: 3dFWHMx not found in PATH." >&2
  exit 2
fi

# Prepare output + a small de-dup ledger to avoid double entries
ensure_header "$OUT"
SEEN_KEYS="$(mktemp "${TMPDIR:-/tmp}/seen_keys.XXXXXX")"
trap 'rm -f "$SEEN_KEYS"' EXIT

# Find L1 FEATs in MNI space only; skip group-level or T1w/other spaces
# We filter paths like:
#   .../derivatives/fsl/sub-XXX/ses-YY/L1_*space-mni*_*.feat
# and avoid any gfeat trees.
find "$FSL_ROOT" -type d -name "*.feat" -print 2>/dev/null \
  | grep -i '/L1_' \
  | grep -i 'space-mni' \
  | grep -v '/gfeat/' \
  | while IFS= read -r FEAT; do

      # Basic fields from path and feat name
      # sub and ses from path components
      sub="$(echo "$FEAT" | sed -n 's#.*/sub-\([A-Za-z0-9]\+\)/.*#\1#p')"
      ses="$(echo "$FEAT" | sed -n 's#.*/ses-\([A-Za-z0-9]\+\)/.*#\1#p')"
      featbase="$(basename "$FEAT")"

      # task and run from the FEAT directory name
      task="$(echo "$featbase" | sed -n 's/.*task-\([A-Za-z0-9-]\+\).*/\1/p')"
      run="$(echo "$featbase" | sed -n 's/.*run-\([0-9]\+\).*/\1/p')"

      # Basic validation
      [ -n "$sub" ] || continue
      [ -n "$ses" ] || continue
      [ -n "$task" ] || continue
      [ -n "$run" ] || continue

      # Whitelist known tasks to suppress irrelevant models
      if ! in_whitelist "$task" $TASK_WHITELIST; then
        continue
      fi

      mask="${FEAT}/mask.nii.gz"
      [ -f "$mask" ] || continue

      # Per-kernel extraction for both multiecho aggregate and single-echo echo-2
      for k in $KERNELS; do
        # --- multiecho aggregate (no _echo-*) ---
        img_me="$(find_me_agg "$sub" "$ses" "$task" "$run" "$k" || true)"
        if [ -n "$img_me" ] && [ -f "$img_me" ]; then
          key="${sub}|${ses}|${task}|${run}|multiecho|${k}"
          if ! grep -qxF "$key" "$SEEN_KEYS"; then
            fwhm="$(get_fwhm_eff "$mask" "$img_me")"
            if [ -n "$fwhm" ]; then
              printf "%s\t%s\t%s\t%s\tmultiecho\t%s\t%s\n" \
                "$sub" "$ses" "$task" "$run" "$k" "$fwhm" >> "$OUT"
              echo "$key" >> "$SEEN_KEYS"
            fi
          fi
        fi

        # --- single-echo echo-2 ---
        img_se="$(find_se_e2 "$sub" "$ses" "$task" "$run" "$k" || true)"
        if [ -n "$img_se" ] && [ -f "$img_se" ]; then
          key="${sub}|${ses}|${task}|${run}|single-echo|${k}"
          if ! grep -qxF "$key" "$SEEN_KEYS"; then
            fwhm_se="$(get_fwhm_eff "$mask" "$img_se")"
            if [ -n "$fwhm_se" ]; then
              printf "%s\t%s\t%s\t%s\tsingle-echo\t%s\t%s\n" \
                "$sub" "$ses" "$task" "$run" "$k" "$fwhm_se" >> "$OUT"
              echo "$key" >> "$SEEN_KEYS"
            fi
          fi
        fi

      done
    done

echo "[$(date +'%F %T')] Done. Wrote: $OUT"
