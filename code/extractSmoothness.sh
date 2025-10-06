#!/usr/bin/env bash
# extractSmoothness.sh
# Compute ACF-based effective FWHM (gaussian_NEWmodel) from AFNI 3dFWHMx
# for both multiecho-combined and single-echo (echo-2) images, at kernels 0 and 5 mm.
# Writes a TSV: sub  ses  task  run  acq  kernel_mm  fwhm_eff  img
# Run with bash (not sh).

set -Eeuo pipefail

# ---------- tiny utils ----------
ts() { date +'%F %T'; }

in_list() {
  local needle=$1; shift
  for x in "$@"; do [[ "$x" == "$needle" ]] && return 0; done
  return 1
}

ok_row() {
  echo "OK   sub=$sub ses=$ses task=$task run=$run acq=$acq k=$k fwhm_eff=$fwhm_eff img=$img"
}

skip_row() {
  echo "SKIP sub=${sub:-NA} ses=${ses:-NA} task=${task:-NA} run=${run:-NA} acq=${acq:-NA} k=${k:-NA} :: $1"
}

# ---------- locate project roots relative to this code file ----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

FSL_ROOT="${FSL_ROOT:-$PROJECT_ROOT/derivatives/fsl}"
FMRIPREP_ROOT="${FMRIPREP_ROOT:-$PROJECT_ROOT/derivatives/fmriprep}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/derivatives/extractions}"
OUT="${OUT:-$OUT_DIR/smoothness_acf.tsv}"

# ---------- knobs ----------
# Only keep these tasks:
TASK_ALLOW=("mid" "sharedreward")
# Smoothing kernels to evaluate:
KERNELS=("0" "5")

# ---------- start ----------
echo "[$(ts)] Starting smoothness extraction"
echo "FSL_ROOT=$FSL_ROOT"
echo "FMRIPREP_ROOT=$FMRIPREP_ROOT"
echo "OUT=$OUT"
echo "KERNELS=${KERNELS[*]}"

mkdir -p "$OUT_DIR"
# header (now includes 'img' as final column)
printf 'sub\tses\ttask\trun\tacq\tkernel_mm\tfwhm_eff\timg\n' > "$OUT"

# ---------- collect L1 feat directories in MNI space (skip T1w, L2, gfeat) ----------
# We only need the FEAT to grab the run-specific mask and basic BIDS fields.
mapfile -t FEATS < <(
  find "$FSL_ROOT" -type d -name '*.feat' -ipath '*/L1_*' \
    ! -ipath '*/L2_*' ! -ipath '*/subject-level/*' ! -ipath '*/group-level/*' \
    ! -ipath '*/*.gfeat/*' ! -iname '*space-t1w*' \
    \( -iname '*space-mni*' -o -iname '*space-mni152*' \) 2>/dev/null \
  | sort -u
)

declare -A seen

# ---------- main loop ----------
for feat in "${FEATS[@]}"; do
  name="$(basename "$feat")"

  # parse sub/ses/task/run from the FEAT name (order can vary)
  sub=""; ses=""; task=""; run=""
  IFS='_' read -r -a toks <<< "${name%.feat}"
  for tok in "${toks[@]}"; do
    case "$tok" in
      sub-*)  sub="${tok#sub-}" ;;
      ses-*)  ses="${tok#ses-}" ;;
      task-*) task="${tok#task-}" ;;
      run-*)  run="${tok#run-}" ;;
      space-*)
        # safety: skip any T1w FEAT that slipped through
        [[ "${tok,,}" == *"t1w"* ]] && { skip_row "space=T1w (skip) :: $feat"; continue 2; }
        ;;
    esac
  done

  # Validate required fields
  if [[ -z "$sub" || -z "$ses" || -z "$task" || -z "$run" ]]; then
    skip_row "malformed FEAT name; cannot parse sub/ses/task/run :: $feat"
    continue
  fi

  # Only desired tasks
  in_list "$task" "${TASK_ALLOW[@]}" || { skip_row "task '$task' not in allowlist ${TASK_ALLOW[*]}"; continue; }

  mask="$feat/mask.nii.gz"
  [[ -f "$mask" ]] || { skip_row "missing mask: $mask"; continue; }

  # We'll process each run x acquisition (multiecho + singleecho) exactly once here;
  # kernels loop inside.
  for acq in multiecho singleecho; do
    key="${sub}|${ses}|${task}|${run}|${acq}"
    if [[ -n "${seen[$key]:-}" ]]; then
      skip_row "duplicate key (already processed) $key"
      continue
    fi
    seen[$key]=1

    # Build base BIDS stem
    stem="${FMRIPREP_ROOT}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}"

    # Loop kernels: derive filename suffix, check existence, compute FWHM
    for k in "${KERNELS[@]}"; do
      if [[ "$k" == "0" ]]; then
        sm_suffix=""
      else
        sm_suffix="_${k}mm"
      fi

      if [[ "$acq" == "singleecho" ]]; then
        img="${stem}_echo-2_part-mag_space-MNI152NLin6Asym_desc-preproc_bold${sm_suffix}.nii.gz"
      else
        # multiecho-combined (no echo-2 tag)
        img="${stem}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold${sm_suffix}.nii.gz"
      fi

      [[ -f "$img" ]] || { skip_row "missing image: $img"; continue; }

      # Run 3dFWHMx (quiet files), parse last numeric line last field (effective FWHM)
      if ! out="$(3dFWHMx -detrend -acf NULL -mask "$mask" -input "$img" 2>&1)"; then
        skip_row "3dFWHMx failed"
        continue
      fi

      # Pull the last numeric line, last column (effective FWHM from ACF fit)
      fwhm_eff="$(grep -E '^[[:space:]]*[0-9.+-Ee]+' <<<"$out" | tail -n 1 | awk '{print $NF}')"

      if [[ -z "${fwhm_eff:-}" ]]; then
        skip_row "could not parse fwhm_eff from 3dFWHMx output"
        continue
      fi

      ok_row
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sub" "$ses" "$task" "$run" "$acq" "$k" "$fwhm_eff" "$img" >> "$OUT"
    done
  done
done

echo "[$(ts)] Done. Wrote: $OUT"
