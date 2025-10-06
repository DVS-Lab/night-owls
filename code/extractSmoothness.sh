#!/usr/bin/env bash
# extractSmoothness.sh
# - Anchors paths to repo's code/ directory
# - Scans L1, MNI-space FEATs (model-1, type-act) only
# - For each run, measures FWHM_eff (gaussian_NEWmodel) from AFNI 3dFWHMx (-acf)
# - Handles both acquisitions:
#     * multiecho   -> fMRIPrep: ..._space-MNI152NLin6Asym_desc-preproc_bold[_<k>mm].nii.gz
#     * singleecho  -> prefer ..._echo-2_part-mag_space-MNI152NLin6Asym_..., else ..._part-mag_space-...
# - Outputs TSV: sub  ses  task  run  acq  kernel_mm  fwhm_eff

set -e
# Avoid -u here; some envs export unbound vars and we don't want a hard crash.

# ---------- Roots anchored to this script ----------
scriptdir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
projdir="$(dirname "$scriptdir")"    # project root containing {code,derivatives}

FSL_ROOT="${projdir}/derivatives/fsl"
FMRIPREP_ROOT="${projdir}/derivatives/fmriprep"
OUT="${projdir}/derivatives/extractions/smoothness_acf.tsv"

mkdir -p "$(dirname "$OUT")"

# ---------- Config ----------
# kernels to evaluate (0 = unsmoothed, 5 = smoothed 5mm)
KERNELS=("0" "5")

# Only these tasks
TASK_ALLOW="mid sharedreward"

# ---------- Helpers ----------
in_list() {
  # usage: in_list "needle" "hay hay2 ..."
  local x="$1"; shift
  for t in $*; do [[ "$x" == "$t" ]] && return 0; done
  return 1
}

# Pick fMRIPrep image for given acq/k (prefer the most specific candidate).
pick_fmriprep_img() {
  local sub="$1" ses="$2" task="$3" run="$4" acq="$5" k="$6"
  local dir="${FMRIPREP_ROOT}/sub-${sub}/ses-${ses}/func"
  local base="sub-${sub}_ses-${ses}_task-${task}_run-${run}_"
  local sm_sfx
  if [[ "$k" == "0" ]]; then
    sm_sfx="desc-preproc_bold.nii.gz"
  else
    sm_sfx="desc-preproc_bold_${k}mm.nii.gz"
  fi

  local cands=()
  if [[ "$acq" == "singleecho" ]]; then
    cands+=( "${dir}/${base}echo-2_part-mag_space-MNI152NLin6Asym_${sm_sfx}" )
    cands+=( "${dir}/${base}part-mag_space-MNI152NLin6Asym_${sm_sfx}" )
  else
    # multiecho (tedana-denoised lives without part-mag/echo tag in fMRIPrep outputs)
    cands+=( "${dir}/${base}space-MNI152NLin6Asym_${sm_sfx}" )
    # fallback, in case some runs were exported with part-mag label
    cands+=( "${dir}/${base}part-mag_space-MNI152NLin6Asym_${sm_sfx}" )
  fi

  for p in "${cands[@]}"; do
    [[ -r "$p" ]] && { echo "$p"; return 0; }
  done
  return 1
}

# Extract FWHM_eff (4th number of the 2nd numeric line from -acf)
measure_fwhm_eff() {
  local mask="$1" img="$2"
  # Quick grid check to avoid AFNI fatal
  if ! 3dinfo -same_grid "$mask" "$img" >/dev/null 2>&1; then
    echo "GRID_MISMATCH"
    return 0
  fi
  local eff
  # -acf prints two numeric lines (classic FWHM; then a b c FWHM_eff)
  if ! eff="$(3dFWHMx -detrend -acf -mask "$mask" -input "$img" 2>/dev/null | awk 'NR==2{print $4; exit}')" ; then
    echo "FAIL"
    return 0
  fi
  [[ -z "$eff" ]] && eff="FAIL"
  echo "$eff"
}

# ---------- Write header (fresh each run) ----------
printf "sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff\n" > "$OUT"

echo "[$(date +'%F %T')] Starting smoothness extraction"
echo "FSL_ROOT=$FSL_ROOT"
echo "FMRIPREP_ROOT=$FMRIPREP_ROOT"
echo "OUT=$OUT"

# ---------- Scan only L1, model-1, type-act, MNI-space FEATs ----------
# This avoids L2 and trial-level LSS FEATs and excludes T1w FEATs.
# Example match:
#   L1_sub-104_ses-03_task-mid_model-1_type-act_run-1_space-mni_single-echo_cnfds-fmriprep.feat
#   L1_sub-104_ses-03_task-mid_model-1_type-act_run-1_space-mni_multi-echo_cnfds-tedana.feat
declare -A seen  # key=sub|ses|task|run|acq to dedupe per run/acq

while IFS= read -r -d '' feat; do
  # Derive fields from FEAT path robustly
  bn="$(basename "$feat")"

  # Skip any FEATs that are not "space-mni"
  [[ "$bn" =~ space-mni ]] || continue

  # Parse sub/ses/task/run (order can vary in some names, so match each independently)
  [[ "$bn" =~ sub-([a-zA-Z0-9]+) ]] && sub="${BASH_REMATCH[1]}" || sub=""
  [[ "$bn" =~ ses-([0-9]+)        ]] && ses="${BASH_REMATCH[1]}" || ses=""
  [[ "$bn" =~ task-([a-zA-Z0-9]+) ]] && task="${BASH_REMATCH[1]}" || task=""
  [[ "$bn" =~ run-([0-9]+)        ]] && run="${BASH_REMATCH[1]}" || run=""

  # Acquisition
  if [[ "$bn" =~ single-echo ]]; then
    acq="singleecho"
  elif [[ "$bn" =~ multi-echo ]]; then
    acq="multiecho"
  else
    # Unlabeled (rare) -> infer multi-echo as safer default for tedana FEATs
    acq="multiecho"
  fi

  # Validate required fields
  if [[ -z "$sub" || -z "$ses" || -z "$task" || -z "$run" ]]; then
    # keep quiet; malformed FEAT directory
    continue
  fi
  # Only desired tasks
  in_list "$task" "$TASK_ALLOW" || continue

  key="${sub}|${ses}|${task}|${run}|${acq}"
  # Process each run/acq exactly once (we will loop kernels inside)
  if [[ -n "${seen[$key]:-}" ]]; then
    continue
  fi
  seen[$key]=1

  mask="${feat}/mask.nii.gz"
  if [[ ! -r "$mask" ]]; then
    echo "WARN skip sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} (missing mask)"
    continue
  fi

  # For each kernel (0 = unsmoothed; 5 = smoothed)
  for k in "${KERNELS[@]}"; do
    img="$(pick_fmriprep_img "$sub" "$ses" "$task" "$run" "$acq" "$k" || true)"
    if [[ -z "$img" ]]; then
      # Do not spam; keep a single concise line
      echo "WARN miss sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} k=${k} (no fMRIPrep file)"
      continue
    fi

    eff="$(measure_fwhm_eff "$mask" "$img")"
    case "$eff" in
      FAIL)
        echo "WARN skip sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} k=${k} (3dFWHMx failed)"
        continue
        ;;
      GRID_MISMATCH)
        echo "WARN grid mismatch sub=${sub} ses=${ses} task=${task} run=${run} acq=${acq} k=${k}"
        continue
        ;;
      *)
        # Append one clean TSV row
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sub" "$ses" "$task" "$run" "$acq" "$k" "$eff" >> "$OUT"
        ;;
    esac
  done

done < <(find "$FSL_ROOT" -type d -name "L1_*_model-1_*_space-mni_*_cnfds-*.feat" -print0)

echo "[$(date +'%F %T')] Done."
