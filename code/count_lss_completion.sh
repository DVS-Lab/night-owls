#!/usr/bin/env bash
# count_lss_completion.sh — run from code/

set -Eeuo pipefail
shopt -s nullglob

# Resolve paths relative to this script (code/ -> project root)
scriptdir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
projectdir="$(dirname "$scriptdir")"
fslroot="${projectdir}/derivatives/fsl"

# Expected trial counts per task (edit if yours differ)
exp_mid=56
exp_sharedreward=54

# Where to look inside each LSS dir
TARGET_GLOB='zstat_trial-*.nii.gz'

# Print detail header
printf "sub\tses\trun\ttask\tacq\tspace\tconfounds\tfound\texpected\tpct\tpath\n"

# Temp file to hold detail rows (TSV)
details="$(mktemp)"

# Walk all LSS dirs: sub-*/LSS_task-*
while IFS= read -r -d '' dir; do
  base="$(basename "$dir")"

  # Expect: LSS_task-<task>_sub-<sub>_ses-<ses>_run-<run>_acq-<acq>_space-<space>_confounds-<confounds>_sm-<s>
  if [[ "$base" =~ ^LSS_task-([^_]+)_sub-([0-9]+)_ses-([0-9]+)_run-([0-9]+)_acq-([^_]+)_space-([^_]+)_confounds-([^_]+) ]]; then
    task="${BASH_REMATCH[1]}"
    sub="${BASH_REMATCH[2]}"
    ses="${BASH_REMATCH[3]}"
    run="${BASH_REMATCH[4]}"
    acq="${BASH_REMATCH[5]}"
    space="${BASH_REMATCH[6]}"
    conf="${BASH_REMATCH[7]}"
  else
    continue
  fi

  # Count files
  files=( "$dir"/$TARGET_GLOB )
  found="${#files[@]}"

  # Expected per task
  case "$task" in
    mid)          expected="$exp_mid" ;;
    sharedreward) expected="$exp_sharedreward" ;;
    *)            expected=0 ;;
  esac

  # Fallback: infer expected from highest observed trial index
  if (( expected <= 0 )); then
    max=0
    for f in "${files[@]}"; do
      n="${f##*/zstat_trial-}"; n="${n%%.nii.gz}"
      [[ "$n" =~ ^[0-9]+$ ]] && (( n>max )) && max="$n"
    done
    expected="$max"
  fi

  # Percent (use a tiny awk just for arithmetic; works even with old awk)
  if (( expected > 0 )); then
    pct=$(awk -v a="$found" -v b="$expected" 'BEGIN{printf "%.1f%%",(a/b)*100}')
  else
    pct="NA"
  fi

  line=$(printf "sub-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s" \
    "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf" "$found" "$expected" "$pct" "$dir")

  # Echo to screen and save for summary
  printf "%s\n" "$line"
  printf "%s\n" "$line" >> "$details"

done < <(find "$fslroot" -maxdepth 2 -type d -name 'LSS_task-*' -print0 | sort -z)

# ---------- Summary without awk arrays ----------
# Sort by the grouping key: task (4), acq (5), space (6), confounds (7)
# Append a sentinel that sorts last to flush the final group.
{
  tail -n +2 "$details"
  printf "sub-\t\t\zzzz\tzzzz\tzzzz\tzzzz\t0\t0\t0\t/\n"
} | LC_ALL=C sort -t$'\t' -k4,4 -k5,5 -k6,6 -k7,7 > "${details}.sorted"

echo
echo "SUMMARY (task acq space confounds found expected  pct)"

prev_key=""
sum_found=0
sum_expected=0

# Read sorted rows and roll up when the key changes
while IFS=$'\t' read -r sub ses run task acq space conf found expected pct path; do
  key="${task}\t${acq}\t${space}\t${conf}"

  if [[ -n "$prev_key" && "$key" != "$prev_key" ]]; then
    # Print the completed group
    printf "%b\t%d\t%d\t" "$prev_key" "$sum_found" "$sum_expected"
    awk -v a="$sum_found" -v b="$sum_expected" 'BEGIN{ if (b>0) printf "%.1f%%\n",(a/b)*100; else print "0.0%"}'
    # Reset accumulators
    sum_found=0
    sum_expected=0
  fi

  # Accumulate
  [[ "$task" == "zzzz" ]] || {
    sum_found=$(( sum_found + found ))
    sum_expected=$(( sum_expected + expected ))
  }

  prev_key="$key"
done < "${details}.sorted"

rm -f "$details" "${details}.sorted"
