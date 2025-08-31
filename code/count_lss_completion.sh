#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob globstar

# --- anchor to this script (no cd) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
bidsdir="$maindir/bids"
deriv_fsl="$maindir/derivatives/fsl"

# --- config ---
tasks=(mid sharedreward)
acqs=(multiecho single)
spaces=(MNI152NLin6Asym T1w)
confounds=(base tedana)
subs=(101 103 104 105)
sessions=({01..12})
runs_per_session=2

# Skip specific (sub,ses) pairs
declare -A SKIP=(
  ["101:04"]=1
  ["101:05"]=1
  ["101:12"]=1
  ["103:12"]=1
)

# Return 0 (done) if output for this job exists, else 1.
# Tweak the sentinels below to match your LSS outputs if needed.
is_done() {
  local task="$1" acq="$2" space="$3" conf="$4" sub="$5" ses="$6" run="$7"

  # Canonical location (edit if your layout differs):
  local base="$deriv_fsl/L1stats_LSS/$task/$acq/$space/$conf/sub-$sub/ses-$ses/run-$run"
  [[ -d "${base}.feat" ]] && return 0
  [[ -f "$base/stats/cope1.nii.gz" || -f "$base/stats/pe1.nii.gz" || -f "$base/design.mat" ]] && return 0

  # Fallback: find any .feat matching the BIDS-style stem anywhere under derivatives/fsl
  local m
  for m in "$deriv_fsl"/**/"sub-${sub}_ses-${ses}_task-${task}_run-${run}"*.feat; do
    [[ -d "$m" ]] && return 0
  done

  return 1
}

# Precompute expected-per-task (constant across rows for that task)
declare -A EXPECTED
for task in "${tasks[@]}"; do
  pairs=0
  for sub in "${subs[@]}"; do
    for ses in "${sessions[@]}"; do
      [[ -n "${SKIP[$sub:$ses]:-}" ]] && continue
      ((pairs++))
    done
  done
  EXPECTED["$task"]=$(( pairs * runs_per_session ))
done

# --- summary header ---
printf "task\tacq\tspace\tconfounds\tfound\texpected\tpct\n"

# Collect missing jobs to list after the summary
missing=()

for task in "${tasks[@]}"; do
  expected="${EXPECTED[$task]}"
  for acq in "${acqs[@]}"; do
    for space in "${spaces[@]}"; do
      for conf in "${confounds[@]}"; do
        found=0
        for sub in "${subs[@]}"; do
          for ses in "${sessions[@]}"; do
            [[ -n "${SKIP[$sub:$ses]:-}" ]] && continue
            for run in $(seq 1 "$runs_per_session"); do
              if is_done "$task" "$acq" "$space" "$conf" "$sub" "$ses" "$run"; then
                ((found++))
              else
                # keep a detailed missing record
                missing+=( "$(printf "sub-%s\tses-%s\t%s\trun-%d\t%s\t%s\t%s" "$sub" "$ses" "$task" "$run" "$acq" "$space" "$conf")" )
              fi
            done
          done
        done

        if (( expected > 0 )); then
          pct=$(awk -v f="$found" -v e="$expected" 'BEGIN{ printf "%.1f%%", (f/e)*100 }')
        else
          pct="NA"
        fi
        printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\n" "$task" "$acq" "$space" "$conf" "$found" "$expected" "$pct"
      done
    done
  done
done

# --- detailed missing list (after summary) ---
echo -e "\nMissing (sub\tses\ttask\trun\tacq\tspace\tconfounds):"
if ((${#missing[@]}==0)); then
  echo "None"
else
  printf "%s\n" "${missing[@]}" | sort -u
fi
