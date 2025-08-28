#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# --- paths (run from code/, resolve maindir one level up) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
bidsdir="$maindir/bids"
deriv_fsl="$maindir/derivatives/fsl"

# --- config ---
tasks=(mid sharedreward)
acqs=(multiecho single)
spaces=(MNI152NLin6Asym T1w)
confounds=(base tedana)   # include both
runs_per_session=2

trials_for_task() {
  case "$1" in
    mid) echo 54 ;;
    sharedreward) echo 56 ;;
    *) echo 0 ;;
  esac
}

# Expected per task = ( #unique sub-ses in BIDS that have this task ) * runs_per_session * trials_per_task
expected_for_task() {
  local task="$1"
  local ntrials; ntrials="$(trials_for_task "$task")"
  [[ "$ntrials" -gt 0 ]] || { echo 0; return; }

  # Gather unique sub-ses that have ANY bold for this task (any acq/echo), from BIDS
  local keys=()
  local f sub ses key
  for f in "$bidsdir"/sub-*/ses-*/func/*task-${task}_run-*_bold.nii.gz; do
    # extract sub-XXX and ses-YY from the path
    sub="$(basename "$(dirname "$(dirname "$(dirname "$f")")")")"   # sub-###
    ses="$(basename "$(dirname "$(dirname "$f")")")"                # ses-##
    key="${sub}_${ses}"
    keys+=("$key")
  done

  if [[ ${#keys[@]} -eq 0 ]]; then
    echo 0
  else
    local uniq_ses
    uniq_ses="$(printf "%s\n" "${keys[@]}" | sort -u | wc -l | tr -d ' ')"
    echo $(( uniq_ses * runs_per_session * ntrials ))
  fi
}

echo "SUMMARY (task acq space confounds found expected  pct)"

for task in "${tasks[@]}"; do
  expected="$(expected_for_task "$task")"

  for acq in "${acqs[@]}"; do
    for space in "${spaces[@]}"; do
      for conf in "${confounds[@]}"; do
        # Count completed LSS zstat files for this slice
        files=( "$deriv_fsl"/sub-*/LSS_task-${task}_sub-*_ses-*_run-*_acq-${acq}_space-${space}_confounds-${conf}_sm-5/zstat_trial-*.nii.gz )
        found="${#files[@]}"

        # pct relative to task-wide expected (constant across rows for the task)
        if (( expected > 0 )); then
          pct="$(awk -v f="$found" -v e="$expected" 'BEGIN{ printf "%.1f%%", (f/e)*100 }')"
        else
          pct="NA"
        fi

        printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\n" \
          "$task" "$acq" "$space" "$conf" "$found" "$expected" "$pct"
      done
    done
  done
done
