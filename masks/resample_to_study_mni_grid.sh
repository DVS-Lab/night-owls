#!/usr/bin/env bash
set -euo pipefail

# This script lives in .../masks
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
fmriprepdir="${maindir}/derivatives/fmriprep"
maskdir="${scriptdir}"

# ---- choose the study MNI grid (3-D boldref) ----
REF_MNI="${1:-}"
if [[ -z "${REF_MNI}" ]]; then
  # auto-pick the first MNI boldref in the study (you can change this to a specific one)
  REF_MNI="$(find "${fmriprepdir}" -type f -name '*space-MNI152NLin6Asym_boldref.nii.gz' | sort | head -n1 || true)"
fi
if [[ -z "${REF_MNI}" || ! -f "${REF_MNI}" ]]; then
  echo "ERROR: Could not locate an MNI boldref. Pass one explicitly, e.g.:"
  echo "  bash $(basename "$0") /path/to/sub-XXX_ses-YY_task-ZZ_run-1_space-MNI152NLin6Asym_boldref.nii.gz"
  exit 1
fi
echo "Using study MNI grid: ${REF_MNI}"

# ---- masks to resample (add more here if needed) ----
masks=(
  "${maskdir}/BrainRewardSignature_2mm.nii"
  "${maskdir}/VS-Imanova_2mm.nii"
)

# helper to strip the trailing "_2mm" to make a clean descriptor
label_from_base () {
  local base="$1"
  base="${base%.nii}"
  echo "${base%_2mm}"
}

# tolerance checker for voxel sizes
vox_ok () {
  local img="$1" t=0.001
  local p1 p2 p3
  p1=$(fslval "$img" pixdim1); p2=$(fslval "$img" pixdim2); p3=$(fslval "$img" pixdim3)
  awk -v a="$p1" -v b="$p2" -v c="$p3" -v t="$t" '
    BEGIN {
      dx = (a-2.7); dy = (b-2.7); dz = (c-2.97);
      if (dx<0) dx=-dx; if (dy<0) dy=-dy; if (dz<0) dz=-dz;
      exit( (dx<=t && dy<=t && dz<=t) ? 0 : 1 );
    }'
}

for m in "${masks[@]}"; do
  if [[ ! -f "$m" ]]; then
    echo "WARN: mask not found → $m"
    continue
  fi

  base="$(basename "$m")"
  label="$(label_from_base "$base")"
  out="${maskdir}/space-MNI152NLin6Asym_desc-${label}_mask.nii.gz"

  echo "Resampling ${base} → $(basename "$out")"
  antsApplyTransforms -d 3 \
    -i "$m" \
    -r "$REF_MNI" \
    -o "$out" \
    -n NearestNeighbor \
    -t identity

  # sanity checks
  if [[ ! -s "$out" ]]; then
    echo "ERROR: output not created: $out"; exit 2
  fi
  read -r min max < <(fslstats "$out" -R)
  if [[ "${max:-0}" == "0" || "${max:-0}" == "0.000000" ]]; then
    echo "ERROR: output appears empty (max=0): $out"; exit 3
  fi
  if ! vox_ok "$out"; then
    p1=$(fslval "$out" pixdim1); p2=$(fslval "$out" pixdim2); p3=$(fslval "$out" pixdim3)
    echo "NOTE: output voxels are ($p1,$p2,$p3); study grid expects (~2.7,2.7,2.97)."
    echo "      This is fine if your chosen reference uses a different spacing."
  fi
done

echo "Done."
