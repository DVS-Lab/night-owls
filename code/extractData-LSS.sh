#!/usr/bin/env bash

# ------------------------------------------------------------
# LSS extractor:
#   - NAcc means (zstat/cope/varcope)
#   - BRS_Cortical_3pt1 means (zstat/cope/varcope)
#   - BRS correlation (zstat vs BRS map)
# ------------------------------------------------------------

# Always run from the code directory
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

deriv_fsl="$maindir/derivatives/fsl"
deriv_fmriprep="$maindir/derivatives/fmriprep"
maskdir="$maindir/masks"
outdir="$maindir/derivatives/extractions"
mkdir -p "$outdir"
outfile="$outdir/extractions_LSS.tsv"

# Match your new output naming
SM_TAG="0"

# Tools check (fail fast if missing)
for cmd in fslstats fslcc fslmaths antsApplyTransforms; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: $cmd not found in PATH." >&2
    exit 2
  fi
done

# Static ROIs / maps in MNI space
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
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\ttrial\tNAcc_zstat_mean\tNAcc_cope_mean\tNAcc_varcope_mean\tBRS_Cort_zstat_mean\tBRS_Cort_cope_mean\tBRS_Cort_varcope_mean\tBRS_corr" \
  > "$outfile"

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

# Pick an existing stat file to use as a reference grid for transforms
pick_ref_file () { # z cope varcope
  local z="$1" c="$2" v="$3"
  if [[ -f "$z" ]]; then echo "$z"
  elif [[ -f "$c" ]]; then echo "$c"
  elif [[ -f "$v" ]]; then echo "$v"
  else echo ""
  fi
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
        for space in MNI152NLin6Asym T1w; do
          for conf in base tedana; do

            combo="$deriv_fsl/sub-${sub}/LSS_task-${task}_sub-${sub}_ses-${ses}_run-${run}_acq-${acq}_space-${space}_confounds-${conf}_sm-${SM_TAG}"

            for (( t=1; t<=ntrials; t++ )); do
              trial=$(printf "%02d" "$t")

              zfile="$combo/zstat_trial-${trial}.nii.gz"
              copefile="$combo/cope_trial-${trial}.nii.gz"
              varcopefile="$combo/varcope_trial-${trial}.nii.gz"

              NAcc_z="NA"; NAcc_c="NA"; NAcc_v="NA"
              Cort_z="NA"; Cort_c="NA"; Cort_v="NA"
              BRS_corr="NA"

              if [[ "$space" == "MNI152NLin6Asym" ]]; then

                [[ -f "$zfile" ]] && NAcc_z=$(fslstats "$zfile" -k "$NAcc_MNI" -M 2>/dev/null || echo "NA")
                [[ -f "$copefile" ]] && NAcc_c=$(fslstats "$copefile" -k "$NAcc_MNI" -M 2>/dev/null || echo "NA")
                [[ -f "$varcopefile" ]] && NAcc_v=$(fslstats "$varcopefile" -k "$NAcc_MNI" -M 2>/dev/null || echo "NA")

                [[ -f "$zfile" ]] && Cort_z=$(fslstats "$zfile" -k "$CORT_MNI" -M 2>/dev/null || echo "NA")
                [[ -f "$copefile" ]] && Cort_c=$(fslstats "$copefile" -k "$CORT_MNI" -M 2>/dev/null || echo "NA")
                [[ -f "$varcopefile" ]] && Cort_v=$(fslstats "$varcopefile" -k "$CORT_MNI" -M 2>/dev/null || echo "NA")

                if [[ -f "$zfile" ]]; then
                  featmask="$(feat_mask_for "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf")"
                  wbmask="$(wbmask_for_combo "$combo" "$zfile" "$featmask")"
                  BRS_corr=$(fslcc --noabs -t -1 -m "$wbmask" -p 6 "$zfile" "$BRS_MNI" 2>/dev/null | awk '{print $NF}')
                  [[ -z "$BRS_corr" ]] && BRS_corr="NA"
                fi

              else
                # T1w: transform masks/maps from MNI -> T1w using the stat grid as reference
                refimg="$(pick_ref_file "$zfile" "$copefile" "$varcopefile")"
                xfm="$(xfm_MNI_to_T1w "$sub" "$ses")"

                if [[ -n "$refimg" && -f "$refimg" && -n "$xfm" && -f "$xfm" ]]; then

                  t1qc_dir="$maindir/derivatives/masks_T1w/sub-${sub}/ses-${ses}/run-${run}_acq-${acq}"
                  mkdir -p "$t1qc_dir"

                  NAcc_T1="$t1qc_dir/desc-NAcc_mask_space-T1w_run-${run}_acq-${acq}.nii.gz"
                  Cort_T1="$t1qc_dir/desc-BRS_Cortical_3pt1_mask_space-T1w_run-${run}_acq-${acq}.nii.gz"
                  BRS_T1="$t1qc_dir/desc-BrainRewardSignature_map_space-T1w_run-${run}_acq-${acq}.nii.gz"

                  if [[ ! -f "$NAcc_T1" ]]; then
                    antsApplyTransforms -d 3 -i "$NAcc_MNI" -r "$refimg" -o "$NAcc_T1" -t "$xfm" -n NearestNeighbor >/dev/null
                  fi
                  if [[ ! -f "$Cort_T1" ]]; then
                    antsApplyTransforms -d 3 -i "$CORT_MNI" -r "$refimg" -o "$Cort_T1" -t "$xfm" -n NearestNeighbor >/dev/null
                  fi
                  if [[ ! -f "$BRS_T1" ]]; then
                    antsApplyTransforms -d 3 -i "$BRS_MNI" -r "$refimg" -o "$BRS_T1" -t "$xfm" >/dev/null
                  fi

                  [[ -f "$zfile" ]] && NAcc_z=$(fslstats "$zfile" -k "$NAcc_T1" -M 2>/dev/null || echo "NA")
                  [[ -f "$copefile" ]] && NAcc_c=$(fslstats "$copefile" -k "$NAcc_T1" -M 2>/dev/null || echo "NA")
                  [[ -f "$varcopefile" ]] && NAcc_v=$(fslstats "$varcopefile" -k "$NAcc_T1" -M 2>/dev/null || echo "NA")

                  [[ -f "$zfile" ]] && Cort_z=$(fslstats "$zfile" -k "$Cort_T1" -M 2>/dev/null || echo "NA")
                  [[ -f "$copefile" ]] && Cort_c=$(fslstats "$copefile" -k "$Cort_T1" -M 2>/dev/null || echo "NA")
                  [[ -f "$varcopefile" ]] && Cort_v=$(fslstats "$varcopefile" -k "$Cort_T1" -M 2>/dev/null || echo "NA")

                  if [[ -f "$zfile" ]]; then
                    featmask="$(feat_mask_for "$sub" "$ses" "$run" "$task" "$acq" "$space" "$conf")"
                    wbmask="$(wbmask_for_combo "$combo" "$zfile" "$featmask")"
                    BRS_corr=$(fslcc --noabs -t -1 -m "$wbmask" -p 6 "$zfile" "$BRS_T1" 2>/dev/null | awk '{print $NF}')
                    [[ -z "$BRS_corr" ]] && BRS_corr="NA"
                  fi
                fi
              fi

              printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$sub" "$ses" "$run" "$task" "$space" "$acq" "$conf" "$trial" \
                "$NAcc_z" "$NAcc_c" "$NAcc_v" \
                "$Cort_z" "$Cort_c" "$Cort_v" \
                "$BRS_corr" \
                >> "$outfile"
            done
          done
        done
      done
    done
  done

  done_sess=$((done_sess+1))
  pct=$(( 100 * done_sess / (total_sess>0?total_sess:1) ))
  echo "$(date '+[%F %T]') ${pct}%% of sessions have been completed (${done_sess}/${total_sess})."
done

echo "$(date '+[%F %T]') Done. Wrote: $outfile"
