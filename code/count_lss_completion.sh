#!/usr/bin/env bash
# count_lss_completion.sh — run from code/ (like oour other scripts)

set -Eeuo pipefail
shopt -s nullglob

# Resolve paths relative to this script (code/ -> project root)
scriptdir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
projectdir="$(dirname "$scriptdir")"
fslroot="${projectdir}/derivatives/fsl"   # /home/.../night-owls/derivatives/fsl

# Expected trial counts per task
declare -A EXPECTED=( ["mid"]=56 ["sharedreward"]=54 )

# Output file (tsv) if you want to save; otherwise stdout is fine.
# outtsv="${projectdir}/logs/lss_completion.tsv"

printf "sub\tses\trun\ttask\tacq\tspace\tconfounds\tfound\texpected\tpct\tpath\n"

# Collect detail lines into a temp file for summary
details="$(mktemp)"

# Walk all LSS dirs (sub-*/LSS_task-*)
while IFS= read -r -d '' dir; do
  base="$(basename "$dir")"

  # Parse fields from directory name:
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
    # Skip unexpected names
    continue
  fi

  # Count zstat files
  files=( "$dir"/zstat_trial-*.nii.gz )
  found="${#files[@]}"

  # Expected trials
  expected="${EXPECTED[$task]:-0}"
  if (( expected <= 0 )); then
    # Fallback: infer from highest trial index present
    max=0
    for f in "${files[@]}"; do
      n="${f##*/zstat_trial-}"; n="${n%%.nii.gz}"
      [[ "$n" =~ ^[0-9]+$ ]] && (( n>max )) && max="$n"
    done
    expected="$max"
  fi

  if (( expected > 0 )); then
    pct=$(awk -v a="$found" -v b="$expected" 'BEGIN{printf "%.1f%%",(a/b)*100}')
  else
    pct="NA"
  fi

  printf "sub-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" \
    "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf" "$found" "$expected" "$pct" "$dir" | tee -a "$details" >/dev/null

done < <(find "$fslroot" -maxdepth 2 -type d -name 'LSS_task-*' -print0 | sort -z)

# Summary: task × acq × space × confounds (weighted by expected)
awk -F'\t' 'NR==1{next} {key=$4"\t"$5"\t"$6"\t"$7; got[key]+=$8; exp[key]+=$9}
  END{
    print "";
    print "SUMMARY (task\tacq\tspace\tconfounds\tfound\texpected\tpct)";
    for(k in got){
      pct = (exp[k]>0)? 100.0*got[k]/exp[k] : 0;
      printf "%s\t%d\t%d\t%.1f%%\n", k, got[k], exp[k], pct
    }
  }' <( { printf "hdr\n"; cat "$details"; } )

rm -f "$details"
