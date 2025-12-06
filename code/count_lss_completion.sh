#!/usr/bin/env bash
# Fast LSS completion counter (glob-based)
# Counts trial files like:
#   sub-101/LSS_task-mid_sub-101_ses-01_run-1_acq-single_space-T1w_confounds-tedana_sm-5/zstat_trial-01.nii.gz

set -euo pipefail
shopt -s nullglob

# --- anchor & roots (no cd) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
deriv_fsl="${DERIV_FSL:-"$maindir/derivatives/fsl"}"

[[ -d "$deriv_fsl" ]] || { echo "ERROR: derivatives/fsl not found: $deriv_fsl" >&2; exit 1; }

# --- config ---
tasks=(mid sharedreward)
# Set the per-task trial counts (swap if your mapping is reversed)
declare -A TRIALS=( ["mid"]=56 ["sharedreward"]=54 )

acqs=(multiecho single)
spaces=(MNI152NLin6Asym)
confounds=(base tedana)

subs=(101 103 104 105)
sessions=({01..12})
runs_per_session=2

# Skip specific (sub,ses) pairs entirely
declare -A SKIP=(
  ["101:04"]=1
  ["101:05"]=1
  ["101:12"]=1
  ["103:12"]=1
)

# --- precompute total valid runs across sub×ses (after skips) ---
valid_runs=0
for sub in "${subs[@]}"; do
  for ses in "${sessions[@]}"; do
    [[ -n "${SKIP[$sub:$ses]:-}" ]] && continue
    (( valid_runs += runs_per_session ))
  done
done

printf "task\tacq\tspace\tconfounds\tfound_trials\texpected_trials\tpct\n"

# Track any (task,sub,ses,run) with missing trials in ANY variant
declare -A ANY_MISSING

for task in "${tasks[@]}"; do
  ntrials="${TRIALS[$task]}"
  expected_total=$(( valid_runs * ntrials ))

  for acq in "${acqs[@]}"; do
    for space in "${spaces[@]}"; do
      for conf in "${confounds[@]}"; do

        found_total=0

        # Loop concrete (sub,ses,run) and count trial files with a FAST glob
        for sub in "${subs[@]}"; do
          for ses in "${sessions[@]}"; do
            [[ -n "${SKIP[$sub:$ses]:-}" ]] && continue
            for run in $(seq 1 "$runs_per_session"); do

              files=( "$deriv_fsl"/sub-"$sub"/LSS_task-"$task"_sub-"$sub"_ses-"$ses"_run-"$run"_acq-"$acq"_space-"$space"_confounds-"$conf"_sm-0/zstat_trial-*.nii.gz )
              nfound=${#files[@]}
              (( found_total += nfound ))

              # Mark this (task,sub,ses,run) if ANY variant is incomplete
              if (( nfound < ntrials )); then
                ANY_MISSING["$task:$sub:$ses:$run"]=1
              fi

            done
          done
        done

        # Summary row for this variant
        if (( expected_total > 0 )); then
          pct=$(awk -v f="$found_total" -v e="$expected_total" 'BEGIN{printf "%.1f%%",(f/e)*100}')
        else
          pct="NA"
        fi
        printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\n" "$task" "$acq" "$space" "$conf" "$found_total" "$expected_total" "$pct"

      done
    done
  done
done

# --- list runs with missing trials (any variant incomplete) ---
echo -e "\nMissing (sub\tses\ttask\trun):"
if (( ${#ANY_MISSING[@]} == 0 )); then
  echo "None"
else
  # Print sorted unique (sub ses task run)
  for key in "${!ANY_MISSING[@]}"; do
    IFS=':' read -r t s u r <<< "$key"
    printf "sub-%s\tses-%s\t%s\trun-%s\n" "$s" "$u" "$t" "$r"
  done | sort -V
fi
