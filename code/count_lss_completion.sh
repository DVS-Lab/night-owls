#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# --- anchor to this script (no cd) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
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

# helper: extract value like sub-101 from "LSS_task-..._sub-101_ses-03_run-2_..."
val_from() { local s="$1" key="$2"; s="_${s}_"; s="${s#*_${key}-}"; printf '%s' "${s%%_*}"; }

# expected runs per (task,acq,space,confounds)
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

# Track presence across all (acq,space,conf) so we can list missing by sub/ses/task/run
declare -A PRESENT                    # key: task:sub:ses:run -> count of variants present
expected_variants=$(( ${#acqs[@]} * ${#spaces[@]} * ${#confounds[@]} ))

printf "task\tacq\tspace\tconfounds\tfound\texpected\tpct\n"

for task in "${tasks[@]}"; do
  expected="${EXPECTED[$task]}"
  for acq in "${acqs[@]}"; do
    for space in "${spaces[@]}"; do
      for conf in "${confounds[@]}"; do

        # gather unique (sub,ses,run) with at least one zstat for this variant
        declare -A SEEN=()
        files=( "$deriv_fsl"/sub-*/LSS_task-${task}_sub-*_ses-*_run-*_acq-${acq}_space-${space}_confounds-${conf}_sm-5/zstat_trial-*.nii.gz )
        for f in "${files[@]}"; do
          dir="$(dirname "$f")"
          base="$(basename "$dir")"  # e.g., LSS_task-mid_sub-101_ses-03_run-2_acq-multiecho_space-..._confounds-..._sm-5
          subv="$(val_from "$base" sub)"
          sesv="$(val_from "$base" ses)"
          runv="$(val_from "$base" run)"

          # keep only configured subjects/sessions and not skipped
          keep_sub=0; for s in "${subs[@]}";   do [[ "$s"  == "$subv" ]] && keep_sub=1 && break; done
          keep_ses=0; for s in "${sessions[@]}"; do [[ "$s" == "$sesv" ]] && keep_ses=1 && break; done
          [[ $keep_sub -eq 0 || $keep_ses -eq 0 ]] && continue
          [[ -n "${SKIP[$subv:$sesv]:-}" ]] && continue

          key="$subv:$sesv:$runv"
          SEEN["$key"]=1
        done

        # mark presence for missing-by-(sub,ses,task,run) reporting
        for key in "${!SEEN[@]}"; do
          PRESENT["$task:$key"]=$(( ${PRESENT["$task:$key"]:-0} + 1 ))
        done

        found=${#SEEN[@]}
        pct="NA"; (( expected > 0 )) && pct=$(awk -v f="$found" -v e="$expected" 'BEGIN{printf "%.1f%%",(f/e)*100}')
        printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\n" "$task" "$acq" "$space" "$conf" "$found" "$expected" "$pct"

        unset SEEN
      done
    done
  done
done

# --- detailed missing list (after summary) ---
echo -e "\nMissing (sub\tses\ttask\trun):"
missing_any=0
for task in "${tasks[@]}"; do
  for sub in "${subs[@]}"; do
    for ses in "${sessions[@]}"; do
      [[ -n "${SKIP[$sub:$ses]:-}" ]] && continue
      for run in $(seq 1 "$runs_per_session"); do
        key="$task:$sub:$ses:$run"
        have=${PRESENT[$key]:-0}
        if (( have < expected_variants )); then
          printf "sub-%s\tses-%s\t%s\trun-%d\n" "$sub" "$ses" "$task" "$run"
          missing_any=1
        fi
      done
    done
  done
done
(( missing_any == 0 )) && echo "None"
