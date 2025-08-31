#!/usr/bin/env bash
set -euo pipefail

# -----------------------
# Count LSS completion (robust matcher)
# -----------------------
# Usage: run in any dir. Optionally override root:
#   DERIV_FSL=/path/to/derivatives/fsl ./count_lss_completion.sh
# Set DEBUG=1 to print a few matched files per variant.

# --- locate derivatives/fsl (no cd) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="${scriptdir%/*}"
deriv_fsl="${DERIV_FSL:-"$maindir/derivatives/fsl"}"

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

# quick membership sets for filtering
declare -A ALLOW_SUB ALLOW_SES
for s in "${subs[@]}";     do ALLOW_SUB["$s"]=1; done
for s in "${sessions[@]}"; do ALLOW_SES["$s"]=1; done

# --- helpers ---
# extract value like sub-101 from "LSS_task-..._sub-101_ses-03_run-2_..."
val_from() { local s="$1" key="$2"; s="_${s}_"; s="${s#*_${key}-}"; printf '%s' "${s%%_*}"; }

# expected-per-task (same across variants)
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

        # Collect unique (sub,ses,run) that have at least one zstat for this variant.
        declare -A SEEN=()

        # Use find (more robust than bare globs); allow any smoothing level sm-*
        # Matches e.g.:
        #   .../sub-101/LSS_task-mid_sub-101_ses-03_run-2_acq-multiecho_space-MNI152NLin6Asym_confounds-base_sm-5/zstat_trial-0001.nii.gz
        # Adjust the *_sm-* if needed.
        while IFS= read -r -d '' f; do
          dir="$(dirname "$f")"
          base="$(basename "$dir")"
          subv="$(val_from "$base" sub)"
          sesv="$(val_from "$base" ses)"
          runv="$(val_from "$base" run)"

          # filter to configured sets and skip list
          [[ -z "${ALLOW_SUB[$subv]:-}" || -z "${ALLOW_SES[$sesv]:-}" ]] && continue
          [[ -n "${SKIP[$subv:$sesv]:-}" ]] && continue

          SEEN["$subv:$sesv:$runv"]=1
        done < <(find "$deriv_fsl" -type f -name 'zstat_trial-*.nii.gz' -path \
          "*/LSS_task-${task}_sub-*_ses-*_run-*_acq-${acq}_space-${space}_confounds-${conf}_sm-*/zstat_trial-*.nii.gz" -print0)

        # Optional debug: show a few matched files per variant
        if [[ "${DEBUG:-0}" == "1" ]]; then
          echo "DEBUG example matches for $task/$acq/$space/$conf:"
          find "$deriv_fsl" -type f -name 'zstat_trial-*.nii.gz' -path \
            "*/LSS_task-${task}_sub-*_ses-*_run-*_acq-${acq}_space-${space}_confounds-${conf}_sm-*/zstat_trial-*.nii.gz" \
            | head -n 3
        fi

        # Mark presence for global missing report
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
