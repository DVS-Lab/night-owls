#!/usr/bin/env bash
# count_lss_completion.sh — run from code/

set -Eeuo pipefail
shopt -s nullglob

# Resolve paths relative to this script (code/ -> project root)
scriptdir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
projectdir="$(dirname "$scriptdir")"
fslroot="${projectdir}/derivatives/fsl"

# Expected trial counts per task (edit if needed)
declare -A EXPECTED=( ["mid"]=56 ["sharedreward"]=54 )

# Header (print to screen and to details file)
details="$(mktemp)"
header=$'sub\tses\trun\ttask\tacq\tspace\tconfounds\tfound\texpected\tpct\tpath'
printf "%s\n" "$header"
printf "%s\n" "$header" > "$details"

# Walk all LSS dirs
while IFS= read -r -d '' dir; do
  base="$(basename "$dir")"
  # LSS_task-<task>_sub-<sub>_ses-<ses>_run-<run>_acq-<acq>_space-<space>_confounds-<confounds>_sm-<s>
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

  # Count zstat files
  files=( "$dir"/zstat_trial-*.nii.gz )
  found="${#files[@]}"

  # Expected trials
  expected="${EXPECTED[$task]:-0}"
  if (( expected <= 0 )); then
    max=0
    for f in "${files[@]}"; do
      n="${f##*/zstat_trial-}"; n="${n%%.nii.gz}"
      [[ "$n" =~ ^[0-9]+$ ]] && { (( n>max )) && max="$n"; }
    done
    expected="$max"
  fi

  if (( expected > 0 )); then
    pct=$(awk -v a="$found" -v b="$expected" 'BEGIN{printf "%.1f%%",(a/b)*100}')
  else
    pct="NA"
  fi

  line=$(printf "sub-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s" \
    "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf" "$found" "$expected" "$pct" "$dir")

  printf "%s\n" "$line" | tee -a "$details" >/dev/null
done < <(find "$fslroot" -maxdepth 2 -type d -name 'LSS_task-*' -print0 | sort -z)

# Summary (portable awk: no +=, no ?:)
awk -F'\t' '
  NR==1 { next }  # skip header
  {
    key = $4 "\t" $5 "\t" $6 "\t" $7
    got[key] = (key in got ? got[key] : 0) + $8
    exp[key] = (key in exp ? exp[key] : 0) + $9
  }
  END {
    print ""
    print "SUMMARY (task\tacq\tspace\tconfounds\tfound\texpected\tpct)"
    for (k in got) {
      if (exp[k] > 0) pct = 100.0 * got[k] / exp[k]; else pct = 0
      printf "%s\t%d\t%d\t%.1f%%\n", k, got[k], exp[k], pct
    }
  }
' "$details"

rm -f "$details"
