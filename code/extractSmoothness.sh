#!/usr/bin/env bash

# ------------------------------------------------------------
# Smoothness extractor (ACF effective FWHM) for L1 FEAT runs
# - Works for multi-echo and single-echo (echo-2) acquisitions
# - Uses FEAT mask in MNI space; skips T1w and higher-level dirs
# - Writes one row per (sub, ses, task, run, acq, kernel_mm)
#   where kernel_mm ∈ {0,5} based on fMRIPrep BOLD (raw vs 5mm)
# - Output columns: sub  ses  task  run  acq  kernel_mm  fwhm_eff
# ------------------------------------------------------------

# --- project roots anchored to this script's directory ---
scriptdir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
projdir="$(dirname "$scriptdir")"            # .../<project>/{code,derivatives}

FSL_ROOT="${projdir}/derivatives/fsl"
FMRIPREP_ROOT="${projdir}/derivatives/fmriprep"
OUT="${projdir}/derivatives/extractions/smoothness_acf.tsv"

# ensure output dir exists (create each run, as before)
mkdir -p "$(dirname "$OUT")"

# kernels to evaluate (0=unsmoothed, 5=Gaussian 5mm)
KERNELS=(0 5)

# ---- tiny logger ----
log(){ printf '[%(%F %T)T] %s\n' -1 "$*"; }

# ---- preflight ----
command -v 3dFWHMx >/dev/null 2>&1 || { echo "ERR: 3dFWHMx not in PATH" >&2; exit 2; }
[ -d "$FSL_ROOT" ] || { echo "ERR: FSL_ROOT not found: $FSL_ROOT" >&2; exit 2; }
[ -d "$FMRIPREP_ROOT" ] || { echo "ERR: FMRIPREP_ROOT not found: $FMRIPREP_ROOT" >&2; exit 2; }

mkdir -p "$(dirname "$OUT")" || { echo "ERR: cannot create $(dirname "$OUT")" >&2; exit 2; }
# (Re)create TSV each run with header
printf 'sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff\n' > "$OUT" || { echo "ERR: cannot write $OUT" >&2; exit 2; }

log "Starting smoothness extraction"
echo "FSL_ROOT=$FSL_ROOT"
echo "FMRIPREP_ROOT=$FMRIPREP_ROOT"
echo "OUT=$OUT"
echo "KERNELS=${KERNELS[*]}"

# ---- helpers ----
extract_tag () {
  # usage: extract_tag <string> <tagname>
  # returns value after tag- up to next underscore
  local s="$1" tag="$2"
  echo "$s" | grep -o "${tag}-[^_/]*" | head -n1 | cut -d- -f2-
}

fwhm_from_img () {
  # usage: fwhm_from_img <mask> <img>
  local mask="$1" img="$2"
  # Use '-acf NULL' to avoid creating 1D/PNG; parse last numeric field from last numeric line
  3dFWHMx -detrend -acf NULL -mask "$mask" -input "$img" 2>/dev/null \
    | awk '/^[0-9eE.+-]/ { last=$NF } END{ if (last=="") { exit 1 } else { printf "%.5f\n", last } }'
}

# seen set for dedup (key = sub|ses|task|run|acq)
declare -A SEEN

# ---- scan L1 FEATs in MNI space only; exclude group-level and T1w ----
# We keep it tight to avoid the flood you saw previously.
# Pattern notes:
#  - L1_* to capture run-level
#  - space-*mni* (case-insensitive later)
#  - *_echo* so we touch both single-echo and multi-echo
find "$FSL_ROOT" -type d -name "L1_*" -name "*echo*.feat" \
  -not -path "*gfeat*" -not -path "*subject-level*" -print0 \
| while IFS= read -r -d '' FEAT; do
    base="$(basename "$FEAT")"
    # basic tags
    sub="$(extract_tag "$FEAT" sub)"
    ses="$(extract_tag "$FEAT" ses)"
    task="$(extract_tag "$FEAT" task)"
    run="$(extract_tag "$FEAT" run)"
    space="$(echo "$base" | grep -io 'space-[^_]*' | head -n1 | cut -d- -f2- | tr '[:upper:]' '[:lower:]')"

    # must have core tags
    if [ -z "$sub" ] || [ -z "$ses" ] || [ -z "$task" ] || [ -z "$run" ]; then
      # Skip oddly named dirs (e.g., trialwise LSS, etc.)
      continue
    fi

    # only MNI space FEATs (to match MNI fMRIPrep bold)
    if echo "$space" | grep -q 't1w'; then
      continue
    fi
    if ! echo "$space" | grep -qi 'mni'; then
      continue
    fi

    # acquisition type (from FEAT name)
    acq="multiecho"
    if echo "$base" | grep -q "single-echo"; then
      acq="singleecho"
    fi

    key="${sub}|${ses}|${task}|${run}|${acq}"
    if [[ -n "${SEEN[$key]:-}" ]]; then
      # de-dup across variant models / QC runs
      continue
    fi
    SEEN[$key]=1

    mask="$FEAT/mask.nii.gz"
    if [ ! -s "$mask" ]; then
      # if mask missing, skip to avoid grid mismatch errors
      continue
    fi

    funcdir="$FMRIPREP_ROOT/sub-${sub}/ses-${ses}/func"
    prefix="${funcdir}/sub-${sub}_ses-${ses}_task-${task}_run-${run}"

    # Build image stems for 0mm and 5mm depending on acq
    if [ "$acq" = "singleecho" ]; then
      stem_raw="${prefix}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold"
    else
      stem_raw="${prefix}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold"
    fi

    # loop kernels
    for k in "${KERNELS[@]}"; do
      if [ "$k" -eq 0 ]; then
        img="${stem_raw}.nii.gz"
      else
        img="${stem_raw}_5mm.nii.gz"
      fi

      # image must exist
      if [ ! -s "$img" ]; then
        # quietly skip missing combinations (e.g., some runs may lack 5mm file)
        continue
      fi

      # compute ACF effective FWHM
      if fwhm=$(fwhm_from_img "$mask" "$img"); then
        # progress line (compact)
        echo "OK sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} k=${k} fwhm=${fwhm}"
        # write one tidy row
        printf '%s\t%s\t%s\t%s\t%s\t%d\t%s\n' \
          "$sub" "$ses" "$task" "$run" "$acq" "$k" "$fwhm" >> "$OUT"
      else
        echo "WARN skip sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} k=${k} (3dFWHMx failed)" >&2
      fi
    done
  done

log "Done. Wrote $(wc -l < "$OUT") lines (including header) to $OUT"
