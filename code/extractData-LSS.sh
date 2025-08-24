#!/usr/bin/env bash
# LSS extractor (single-trial). Runs from night-owls/code

set -u -o pipefail

# --- fixed paths (scripts always live in night-owls/code) ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"                 # -> night-owls
fsldir="$maindir/derivatives/fsl"
fmriprepdir="$maindir/derivatives/fmriprep"
masksdir="$maindir/masks"
outdir="$maindir/derivatives/extractions"
mkdir -p "$outdir"

# --- inputs (MNI-space masks/maps only; no transforms when space=MNI152NLin6Asym) ---
MNI_VS="$masksdir/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"
MNI_BRS="$masksdir/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

# --- tools check (fail fast with clear message) ---
for cmd in fslstats fslcc antsApplyTransforms ; do
  if ! command -v "$cmd" >/dev/null 2>&1 ; then
    echo "[ERROR] $cmd not found in PATH. Load your FSL/ANTs modules and re-run." >&2
    exit 1
  fi
done

# --- output file ---
tsv="$outdir/extractions_LSS.tsv"
echo -e "sub\tses\trun\ttask\tspace\tacq\tconfounds\ttrial\tzstat\tlabel\tVS_mean\tBRS_corr" > "$tsv"

# --- task→expected trial count (single-trial designs) ---
# If this mapping differs, tell me and I'll flip it.
declare -A NTRIALS=( [mid]=56 [sharedreward]=54 )

# --- helper: parse key/value from LSS directory base name ---
# Example base: LSS_task-mid_sub-101_ses-01_run-1_acq-multiecho_space-MNI152NLin6Asym_confounds-tedana_sm-5
parse_fields() {
  local base="$1"
  local k v
  for kv in ${base//_/ } ; do
    case "$kv" in
      task-*)       task="${kv#task-}" ;;
      sub-*)       subj="${kv#sub-}" ;;
      ses-*)         ses="${kv#ses-}" ;;
      run-*)         run="${kv#run-}" ;;
      acq-*)         acq="${kv#acq-}" ;;
      space-*)     space="${kv#space-}" ;;
      confounds-*) confounds="${kv#confounds-}" ;;
    esac
  done
}

# --- helper: warp MNI→T1w once per (sub,ses) and cache results ---
# Produces: $cache/space-T1w_desc-VS-Imanova_mask.nii.gz  (NN)
#           $cache/space-T1w_desc-BrainRewardSignature_map.nii.gz  (Linear)
ensure_T1w_derivs() {
  local subj="$1" ses="$2"
  local cache="$outdir/_cache/sub-${subj}_ses-${ses}"
  local doneflag="$cache/.done"
  mkdir -p "$cache"

  if [[ -f "$doneflag" ]]; then
    echo "$cache"
    return 0
  fi

  local anatdir="$fmriprepdir/sub-${subj}/ses-${ses}/anat"
  local xfm="$anatdir/sub-${subj}_ses-${ses}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"

  # Prefer preproc T1w as reference; fallbacks if needed.
  local ref=
  for cand in \
    "$anatdir/sub-${subj}_ses-${ses}_desc-preproc_T1w.nii.gz" \
    "$anatdir/sub-${subj}_ses-${ses}_desc-preproc_T1w.nii" \
    "$anatdir/sub-${subj}_ses-${ses}_T1w.nii.gz" \
    "$anatdir/sub-${subj}_ses-${ses}_T1w.nii" ; do
    [[ -f "$cand" ]] && ref="$cand" && break
  done

  if [[ ! -f "$xfm" || ! -f "$ref" ]]; then
    # Can't warp for this (sub,ses); leave cache empty and return path anyway.
    echo "$cache"
    return 0
  fi

  # VS mask (binary) — NN interpolation
  if [[ -f "$MNI_VS" && ! -f "$cache/space-T1w_desc-VS-Imanova_mask.nii.gz" ]]; then
    antsApplyTransforms -d 3 -i "$MNI_VS" -r "$ref" -t "$xfm" \
      -n NearestNeighbor -o "$cache/space-T1w_desc-VS-Imanova_mask.nii.gz" >/dev/null
  fi

  # BRS map (continuous) — Linear interpolation
  if [[ -f "$MNI_BRS" && ! -f "$cache/space-T1w_desc-BrainRewardSignature_map.nii.gz" ]]; then
    antsApplyTransforms -d 3 -i "$MNI_BRS" -r "$ref" -t "$xfm" \
      -n Linear -o "$cache/space-T1w_desc-BrainRewardSignature_map.nii.gz" >/dev/null
  fi

  touch "$doneflag"
  echo "$cache"
}

# --- collect unique (sub,ses) with any LSS file to drive progress ---
declare -A SESS_KEYS=()
while IFS= read -r -d '' f; do
  d="$(dirname "$f")"
  b="$(basename "$d")"
  unset task subj ses run acq space confounds
  parse_fields "$b"
  [[ -n "${subj:-}" && -n "${ses:-}" ]] && SESS_KEYS["${subj}-${ses}"]=1
done < <(find "$fsldir"/sub-*/LSS_task-* -maxdepth 1 -type f -name 'zstat_trial-*.nii.gz' -print0 2>/dev/null)

# Sort session keys
mapfile -t ALL_SESS <<<"$(printf "%s\n" "${!SESS_KEYS[@]}" | sort -V)"
total_sessions="${#ALL_SESS[@]}"
done_sessions=0

# --- main: process by session to report progress cleanly ---
for sk in "${ALL_SESS[@]}"; do
  IFS='-' read -r subj ses <<<"$sk"

  # All combo dirs for this (sub,ses)
  mapfile -t COMBOS < <(
    find "$fsldir/sub-${subj}/" -maxdepth 1 -type d -name 'LSS_task-*' -printf '%p\n' 2>/dev/null \
    | awk -v S="$ses" -F'/' '{
         b=$NF;
         if (b ~ "_ses-" S "_") print $0;
       }' \
    | sort -V
  )

  # Prepare T1w-space cached derivatives for this (sub,ses)
  T1W_CACHE="$(ensure_T1w_derivs "$subj" "$ses")"
  VS_T1W="$T1W_CACHE/space-T1w_desc-VS-Imanova_mask.nii.gz"
  BRS_T1W="$T1W_CACHE/space-T1w_desc-BrainRewardSignature_map.nii.gz"

  for combo in "${COMBOS[@]}"; do
    b="$(basename "$combo")"
    unset task run acq space confounds
    parse_fields "$b"

    # Skip if task not in our mapping (rare)
    [[ -z "${NTRIALS[$task]+x}" ]] && continue
    ntrials="${NTRIALS[$task]}"

    for tt in $(seq -w 01 "$ntrials"); do
      zfile="$combo/zstat_trial-${tt}.nii.gz"
      trial="${tt}"
      zstat="1"       # per your note: always zstat1 in LSS
      label=""        # no contrast label for LSS

      VS_mean="NA"
      BRS_corr="NA"

      if [[ -f "$zfile" ]]; then
        # Whole-brain FEAT mask (prefer mask alongside or within any *.feat under combo)
        wbmask="$(find "$combo" -maxdepth 2 -type f -name 'mask.nii.gz' | head -n1)"

        if [[ "$space" == "MNI152NLin6Asym" ]]; then
          # No transforms — use MNI masks directly
          if [[ -f "$MNI_VS" ]]; then
            VS_mean="$(fslstats "$zfile" -k "$MNI_VS" -M 2>/dev/null || echo NA)"
          fi
          if [[ -f "$MNI_BRS" && -n "$wbmask" ]]; then
            # Signed correlation over whole brain, full range
            # fslcc prints one line; last column is the coefficient
            BRS_corr="$(fslcc --noabs -t -1 -m "$wbmask" -p 6 "$zfile" "$MNI_BRS" 2>/dev/null | awk '{print $NF}' )"
            [[ -z "$BRS_corr" ]] && BRS_corr="NA"
          fi
        else
          # space=T1w: use pre-warped, ses-aware T1w derivatives if present
          if [[ -f "$VS_T1W" ]]; then
            VS_mean="$(fslstats "$zfile" -k "$VS_T1W" -M 2>/dev/null || echo NA)"
          fi
          if [[ -f "$BRS_T1W" && -n "$wbmask" ]]; then
            BRS_corr="$(fslcc --noabs -t -1 -m "$wbmask" -p 6 "$zfile" "$BRS_T1W" 2>/dev/null | awk '{print $NF}')"
            [[ -z "$BRS_corr" ]] && BRS_corr="NA"
          fi
        fi
      fi

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$subj" "$ses" "$run" "$task" "$space" "$acq" "$confounds" "$trial" "$zstat" "$label" "$VS_mean" "$BRS_corr" \
        >> "$tsv"
    done
  done

  # --- progress echo at session granularity ---
  done_sessions=$((done_sessions + 1))
  # integer percent with one decimal
  pct="$(awk -v d="$done_sessions" -v t="$total_sessions" 'BEGIN{if(t==0){print 100}else{printf("%.1f", (d*100.0)/t)}}')"
  echo "[`date +'%F %T'`] ${pct}% of sessions have been completed"
done

echo "[`date +'%F %T'`] Done. Wrote: $tsv"
