#!/usr/bin/env bash
set -euo pipefail

# This script lives in .../masks
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
fmriprepdir="${maindir}/derivatives/fmriprep"
maskdir="${scriptdir}"

# Pick the study MNI grid (3-D boldref). Pass as $1 or auto-detect the first one in derivatives.
REF_MNI="${1:-}"
if [[ -z "${REF_MNI}" ]]; then
  REF_MNI="$(find "${fmriprepdir}" -type f -name '*space-MNI152NLin6Asym_boldref.nii.gz' | sort | head -n1 || true)"
fi
[[ -n "${REF_MNI}" && -f "${REF_MNI}" ]] || { echo "ERROR: provide an MNI boldref:  bash $(basename "$0") /path/to/*space-MNI152NLin6Asym_boldref.nii.gz"; exit 1; }
echo "Using study MNI grid: ${REF_MNI}"

# Behavior knobs
: "${LINEAR_ALL:=0}"         # 0 = NN for VS, Linear for BRS; 1 = Linear for both (+threshold VS)
: "${THR:=0.5}"              # threshold for VS when LINEAR_ALL=1

# Inputs (as you listed)
BRS="${maskdir}/BrainRewardSignature_2mm.nii"   # continuous
VS="${maskdir}/VS-Imanova_2mm.nii"              # binary

[[ -f "$BRS" ]] || echo "WARN: missing ${BRS}"
[[ -f "$VS"  ]] || echo "WARN: missing ${VS}"

# Helper
resample () { # src, ref, out, interp
  local src="$1" ref="$2" out="$3" interp="$4"
  antsApplyTransforms -d 3 -i "$src" -r "$ref" -o "$out" -n "$interp" -t identity
}

# BrainRewardSignature (continuous) → Linear
if [[ -f "$BRS" ]]; then
  out_brs="${maskdir}/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"
  echo "Resampling BrainRewardSignature → $(basename "$out_brs") (Linear)"
  resample "$BRS" "$REF_MNI" "$out_brs" Linear
fi

# VS-Imanova (binary) → NN by default, or Linear+threshold if LINEAR_ALL=1
if [[ -f "$VS" ]]; then
  out_vs="${maskdir}/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"
  if [[ "$LINEAR_ALL" == "1" ]]; then
    echo "Resampling VS-Imanova → $(basename "$out_vs") (Linear + threshold ${THR})"
    tmp="${out_vs%.nii.gz}_tmp.nii.gz"
    resample "$VS" "$REF_MNI" "$tmp" Linear
    # binarize after interpolation
    fslmaths "$tmp" -thr "$THR" -bin "$out_vs"
    rm -f "$tmp"
  else
    echo "Resampling VS-Imanova → $(basename "$out_vs") (NearestNeighbor)"
    resample "$VS" "$REF_MNI" "$out_vs" NearestNeighbor
  fi
fi

# Quick sanity checks (only if outputs exist)
check_nonempty () {
  local f="$1"
  [[ -f "$f" ]] || return 0
  read -r _min _max < <(fslstats "$f" -R)
  if [[ "${_max:-0}" == "0" || "${_max:-0}" == "0.000000" ]]; then
    echo "ERROR: output appears empty → $f"; exit 2
  fi
}
check_nonempty "${out_brs:-/dev/null}"
check_nonempty "${out_vs:-/dev/null}"

echo "Done."
