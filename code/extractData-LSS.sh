#!/usr/bin/env bash

# ------------------------------------------------------------
# LSS extractor (MNI-only):
#   - NAcc means (zstat/cope/varcope)
#   - BRS_Cortical_3pt1 means (zstat/cope/varcope)
#   - BRS correlation (zstat vs BRS map)
#   - Last column: expected zstat path (sanity check)
# ------------------------------------------------------------

# Always run from the code directory
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

deriv_fsl="$maindir/derivatives/fsl"
maskdir="$maindir/masks"
outdir="$maindir/derivatives/extractions"
mkdir -p "$outdir"
outfile="$outdir/extractions_LSS_smoothed.tsv"

# Match your new output naming
SM_TAG="5"

# Tools check (fail fast if missing)
for cmd in fslstats fslcc fslmaths; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: $cmd not found in PATH." >&2
    exit 2
  fi
done

# MNI masks/maps
NAcc_MNI="$maskdir/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
CORT_MNI="$maskdir/BRS_Cortical_3pt1.nii.gz"
BRS_MNI="$maskdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

for f in "$NAcc_MNI" "$CORT_MNI" "$BRS_MNI"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing mask/map: $f" >&2; exit 2; }
done

# Trial counts
trials_for () {  # $1 = task
  case "$1" in
    mid) echo 56 ;;
    sharedreward) echo 54 ;;
    *) echo 0 ;;
  esac
}

# Write header
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\ttrial\tNAcc_zstat_mean\tNAcc_cope_mean\tNAcc_varcope_mean\tBRS_Cort_zstat_mean\tBRS_Cort_cope_mean\tBRS_Cort_varcope_mean\tBRS_corr\tinput_zstat" \
  > "$outfile"

# Enumerate (sub,ses) that actually have LSS data (based on trial-01 presence)
# Prune subject-level to avoid permission noise
mapfile -t sess_keys < <(
  find "$deriv_fsl" \
    -path "*/subject-level/*" -prune -o \
    -type f -name 'zstat_trial-01.nii.gz' -print 2>/dev/null \
  | sed -E 's|.*sub-([0-9]+)/.*_ses-([0-9]+)_.*|\1 \2|' \
  | sort -u
)

total_sess="${#sess_keys[@]}"
done_sess=0

# Helper: build/reuse a whole-brain mask for this combo
feat_mask_for () { # sub ses run task acq space conf
  local sub="$1" ses="$2" run="$3" task="$4" acq="$5" space="$6" conf="$7"

  local combo="$deriv_fsl/sub-${sub}/LSS_task-${task}_sub-${sub}_ses-${ses}_run-${run}_acq-${acq}_space-${space}_confounds-${conf}_sm-${SM_TAG}"
  local auto="$combo/wbmask_auto.nii.gz"

  local zref="$combo/zstat_trial-01.nii.gz"
  if [[ ! -f "$zref" ]]; then
    zref=$(ls "$combo"/zstat_trial-*.nii.gz 2>/dev/null | head -n1)
  fi

  if [[ -n "${zref:-}" && -f "$zref" ]]; then
    if [[ ! -f "$auto" ]]; then
      fslmaths "$zref" -abs -thr 0 -bin "$auto" >/dev/null
    fi
    echo "$auto"
  else
    echo ""
  fi
}

wbmask_for_combo () { # combo_dir zfile featmask
  local combo="$1" zfile="$2" featmask="$3"
  if [[ -n "$featmask" && -f "$featmask" ]]; then
    echo "$featmask"
    return
  fi
  local auto="$combo/wbmask_auto.nii.gz"
  if [[ ! -f "$auto" ]]; then
    fslmaths "$zfile" -abs -thr 0 -bin "$auto" >/dev/null
  fi
  echo "$auto"
}

# Main loop over sessions
for key in "${sess_keys[@]}"; do
  sub="${key%% *}"
  ses="${key##* }"

  for task in mid sharedreward; do
    ntrials="$(trials_for "$task")"
    [[ "$ntrials" -gt 0 ]] || continue

    for run in 1 2; do
      for acq in multiecho single; do
        for conf in base tedana; do

          space="MNI152NLin6Asym"

          combo="$deriv_fsl/sub-${sub}/LSS_task-${task}_sub-${sub}_ses-${ses}_run-${run}_acq-${acq}_space-${space}_confounds-${conf}_sm-${SM_TAG}"

          for (( t=1; t<=ntrials; t++ )); do
            trial=$(printf "%02d" "$t")

            zfile="$combo/zstat_trial-${trial}.nii.gz"
            copefile="$combo/cope_trial-${trial}.nii.gz"
            varcopefile="$combo/varcope_trial-${trial}.nii.gz"

            # Default outputs
            NAcc_z="NA"; NAcc_c="NA"; NAcc_v="NA"
            Cort_z="NA"; Cort_c="NA"; Cort_v="NA"
            BRS_corr="NA"

            # ROI means (compute each stat if present)
            [[ -f "$zfile" ]]      && NAcc_z=$(fslstats "$zfile"      -k "$NAcc_MNI" -M 2>/dev/null || echo "NA")
            [[ -f "$copefile" ]]   && NAcc_c=$(fslstats "$copefile"   -k "$NAcc_MNI" -M 2>/dev/null || echo "NA")
            [[ -f "$varcopefile" ]]&& NAcc_v=$(fslstats "$varcopefile"-k "$NAcc_MNI" -M 2>/dev/null || echo "NA")

            [[ -f "$zfile" ]]      && Cort_z=$(fslstats "$zfile"      -k "$CORT_MNI" -M 2>/dev/null || echo "NA")
            [[ -f "$copefile" ]]   && Cort_c=$(fslstats "$copefile"   -k "$CORT_MNI" -M 2>/dev/null || echo "NA")
            [[ -f "$varcopefile" ]]&& Cort_v=$(fslstats "$varcopefile"-k "$CORT_MNI" -M 2>/dev/null || echo "NA")

            # BRS correlation uses zstat (as before conceptually)
            if [[ -f "$zfile" ]]; then
              featmask="$(feat_mask_for "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf")"
              wbmask="$(wbmask_for_combo "$combo" "$zfile" "$featmask")"
              BRS_corr=$(fslcc --noabs -t -1 -m "$wbmask" -p 6 "$zfile" "$BRS_MNI" 2>/dev/null | awk '{print $NF}')
              [[ -z "$BRS_corr" ]] && BRS_corr="NA"
            fi

            # Emit row (input_zstat is last column)
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
              "$sub" "$ses" "$run" "$task" "$space" "$acq" "$conf" "$trial" \
              "$NAcc_z" "$NAcc_c" "$NAcc_v" \
              "$Cort_z" "$Cort_c" "$Cort_v" \
              "$BRS_corr" "$zfile" \
              >> "$outfile"
          done
        done
      done
    done
  done

  # Progress echo at the session level (same idea as your original)
  done_sess=$((done_sess+1))
  pct=$(( 100 * done_sess / (total_sess>0?total_sess:1) ))
  echo "$(date '+[%F %T]') ${pct}%% of sessions have been completed (${done_sess}/${total_sess})."
done

echo "$(date '+[%F %T]') Done. Wrote: $outfile"
