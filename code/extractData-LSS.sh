#!/usr/bin/env bash

# ------------------------------------------------------------
# LSS extractor: VS mean + BRS correlation per trial
# ------------------------------------------------------------

# Always run from the code directory
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

deriv_fsl="$maindir/derivatives/fsl"
deriv_fmriprep="$maindir/derivatives/fmriprep"
masks_dir="$maindir/masks"
outdir="$maindir/derivatives/extractions"
mkdir -p "$outdir"
outfile="$outdir/extractions_LSS.tsv"

# Tools check (fail fast if missing)
for cmd in fslstats fslcc fslmaths antsApplyTransforms; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: $cmd not found in PATH." >&2
    exit 2
  fi
done

# Static ROIs in MNI space
VS_MNI="$masks_dir/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"
BRS_MNI="$masks_dir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"
for f in "$VS_MNI" "$BRS_MNI"; do
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
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\ttrial\tVS_mean\tBRS_corr" > "$outfile"

# Enumerate (sub,ses) that actually have LSS data (based on trial-01 presence anywhere)
mapfile -t sess_keys < <(
  find "$deriv_fsl" -type f -name 'zstat_trial-01.nii.gz' \
  | sed -E 's|.*sub-([0-9]+)/.*_ses-([0-9]+)_.*|\1 \2|' \
  | sort -u
)

total_sess="${#sess_keys[@]}"
done_sess=0

# Helper: build/reuse a whole-brain mask for this LSS (sub/ses/run/task/acq/space/conf) combo
feat_mask_for () { # sub ses run task acq space conf
  local sub="$1" ses="$2" run="$3" task="$4" acq="$5" space="$6" conf="$7"

  # LSS combo directory matches your listing pattern
  local combo="$deriv_fsl/sub-${sub}/LSS_task-${task}_sub-${sub}_ses-${ses}_run-${run}_acq-${acq}_space-${space}_confounds-${conf}_sm-5"

  # Where we store the auto mask
  local auto="$combo/wbmask_auto.nii.gz"

  # Choose a reference zstat to define the grid (prefer trial-01; fall back to the first available)
  local zref="$combo/zstat_trial-01.nii.gz"
  if [[ ! -f "$zref" ]]; then
    zref=$(ls "$combo"/zstat_trial-*.nii.gz 2>/dev/null | head -n1)
  fi

  # If we have a reference, create the mask once; otherwise return blank so caller can handle NA
  if [[ -n "$zref" && -f "$zref" ]]; then
    if [[ ! -f "$auto" ]]; then
      # Binary mask of all nonzero voxels in the zstat grid
      fslmaths "$zref" -abs -thr 0 -bin "$auto" >/dev/null
    fi
    echo "$auto"
  else
    echo ""
  fi
}


# Helper: ensure MNI->T1w transform (prefer *MNI152NLin6Asym*)
xfm_MNI_to_T1w () { # sub ses
  local sub="$1" ses="$2" anatdir="$deriv_fmriprep/sub-${sub}/ses-${ses}/anat"
  local pref="$anatdir/sub-${sub}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"
  if [[ -f "$pref" ]]; then
    echo "$pref"; return
  fi
  local any
  any="$(find "$anatdir" -maxdepth 1 -type f -name "sub-${sub}_ses-${ses}_from-MNI*to-T1w*_xfm.h5" | sort | head -n1)"
  echo "$any"
}

# Helper: build or reuse per-combo whole-brain mask
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

  # Loop structure fixed and consistent with prior scripts
  for task in mid sharedreward; do
    ntrials="$(trials_for "$task")"
    [[ "$ntrials" -gt 0 ]] || continue
    for run in 1 2; do
      for acq in multiecho single; do
        for space in MNI152NLin6Asym T1w; do
          conf="tedana"

          # LSS combo directory (matches your listing)
          combo="$deriv_fsl/sub-${sub}/LSS_task-${task}_sub-${sub}_ses-${ses}_run-${run}_acq-${acq}_space-${space}_confounds-${conf}_sm-5"

          # We will still emit rows with NA if the combo/trial file is missing
          for (( t=1; t<=ntrials; t++ )); do
            trial=$(printf "%02d" "$t")
            zfile="$combo/zstat_trial-${trial}.nii.gz"

            VS_mean="NA"
            BRS_corr="NA"

            if [[ -f "$zfile" ]]; then
              if [[ "$space" == "MNI152NLin6Asym" ]]; then
                # No transforms: use MNI masks/maps directly
                VS_mean=$(fslstats "$zfile" -k "$VS_MNI" -M 2>/dev/null || echo "NA")

                # whole-brain mask preference: FEAT mask from matched L1, else auto
                featmask="$(feat_mask_for "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf")"
                wbmask="$(wbmask_for_combo "$combo" "$zfile" "$featmask")"

                # signed correlation, full range
                BRS_corr=$(fslcc --noabs -t -1 -m "$wbmask" -p 6 "$zfile" "$BRS_MNI" 2>/dev/null | awk '{print $NF}' )
                [[ -z "$BRS_corr" ]] && BRS_corr="NA"

              else
                # T1w: transform VS mask (NN) and BRS map (linear) from MNI -> T1w (zfile grid)
                xfm="$(xfm_MNI_to_T1w "$sub" "$ses")"
                if [[ -z "$xfm" || ! -f "$xfm" ]]; then
                  # If we cannot find a transform, leave NA but keep emitting a row
                  VS_mean="NA"; BRS_corr="NA"
                else
                  # Store transformed masks/maps for QC
                  t1qc_dir="$maindir/derivatives/masks_T1w/sub-${sub}/ses-${ses}/run-${run}_acq-${acq}"
                  mkdir -p "$t1qc_dir"
                  VS_T1="$t1qc_dir/desc-VS-Imanova_mask_space-T1w_run-${run}_acq-${acq}.nii.gz"
                  BRS_T1="$t1qc_dir/desc-BrainRewardSignature_map_space-T1w_run-${run}_acq-${acq}.nii.gz"

                  if [[ ! -f "$VS_T1" ]]; then
                    antsApplyTransforms -d 3 -i "$VS_MNI"  -r "$zfile" -o "$VS_T1"  -t "$xfm" -n NearestNeighbor >/dev/null
                  fi
                  if [[ ! -f "$BRS_T1" ]]; then
                    antsApplyTransforms -d 3 -i "$BRS_MNI" -r "$zfile" -o "$BRS_T1" -t "$xfm" >/dev/null
                  fi

                  VS_mean=$(fslstats "$zfile" -k "$VS_T1" -M 2>/dev/null || echo "NA")

                  featmask="$(feat_mask_for "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf")"
                  wbmask="$(wbmask_for_combo "$combo" "$zfile" "$featmask")"

                  BRS_corr=$(fslcc --noabs -t -1 -m "$wbmask" -p 6 "$zfile" "$BRS_T1" 2>/dev/null | awk '{print $NF}' )
                  [[ -z "$BRS_corr" ]] && BRS_corr="NA"
                fi
              fi
            fi

            # Emit row (no zstat/label columns)
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
              "$sub" "$ses" "$run" "$task" "$space" "$acq" "$conf" "$trial" "$VS_mean" "$BRS_corr" \
              >> "$outfile"
          done
        done
      done
    done
  done

  # Progress echo at the session level
  done_sess=$((done_sess+1))
  pct=$(( 100 * done_sess / (total_sess>0?total_sess:1) ))
  echo "$(date '+[%F %T]') ${pct}%% of sessions have been completed (${done_sess}/${total_sess})."
done

echo "$(date '+[%F %T]') Done. Wrote: $outfile"
