#!/usr/bin/env bash
set -euo pipefail

# --- Paths relative to this script ---
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
fmriprepdir="${maindir}/derivatives/fmriprep"

# --- One-off entities (edit as needed) ---
sub=101
ses=03
task=mid
run=2     # run to fix; run-1 provides the reference grid

# --- Helper: does a 4D file match a 3D reference grid? (dim1-3 + pixdim1-3) ---
matches_ref_grid () {
  local ref="$1" src="$2"
  local r1 r2 r3 s1 s2 s3 rp1 rp2 rp3 sp1 sp2 sp3
  r1=$(fslval "$ref" dim1); r2=$(fslval "$ref" dim2); r3=$(fslval "$ref" dim3)
  s1=$(fslval "$src" dim1); s2=$(fslval "$src" dim2); s3=$(fslval "$src" dim3)
  rp1=$(fslval "$ref" pixdim1); rp2=$(fslval "$ref" pixdim2); rp3=$(fslval "$ref" pixdim3)
  sp1=$(fslval "$src" pixdim1); sp2=$(fslval "$src" pixdim2); sp3=$(fslval "$src" pixdim3)
  [[ "$r1" == "$s1" && "$r2" == "$s2" && "$r3" == "$s3" ]] && \
  awk -v a="$rp1" -v b="$sp1" 'BEGIN{exit ((a-b<0?b-a:a-b)<=1e-6?0:1)}' >/dev/null && \
  awk -v a="$rp2" -v b="$sp2" 'BEGIN{exit ((a-b<0?b-a:a-b)<=1e-6?0:1)}' >/dev/null && \
  awk -v a="$rp3" -v b="$sp3" 'BEGIN{exit ((a-b<0?b-a:a-b)<=1e-6?0:1)}' >/dev/null
}

# --- Core: regrid a 4D series onto a 3D boldref grid (identity transform) ---
regrid_identity () {
  local ref3d="$1" src4d="$2" label="$3"

  if [[ ! -f "$ref3d" ]]; then echo "[$label] Missing reference: $ref3d"; return 1; fi
  if [[ ! -f "$src4d" ]]; then echo "[$label] Missing source:    $src4d"; return 0;  fi

  # 1) Skip if already repaired
  if matches_ref_grid "$ref3d" "$src4d"; then
    echo "[$label] Already matches grid → $src4d"; return 0
  fi

  # 2) Backup original once
  local backup="${src4d%.nii.gz}_ORIGINAL.nii.gz"
  if [[ ! -f "$backup" ]]; then
    cp -p "$src4d" "$backup"
    echo "[$label] Backup created → $backup"
  else
    echo "[$label] Backup exists   → $backup"
  fi

  # 3) Resample with identity (write to temp, then replace)
  local tmp="${src4d%.nii.gz}_tmp.nii.gz"
  echo "[$label] Resampling (identity) → $src4d"
  antsApplyTransforms -d 3 -e 3 --float \
    -i "$backup" \
    -r "$ref3d" \
    -o "$tmp" \
    -n Linear \
    -t identity

  # 4) Sanity: nonzero output and expected voxel sizes (~2.7,2.7,2.97)
  read -r _min _max < <(fslstats "$tmp" -R)
  if [[ "${_max:-0}" == "0" || "${_max:-0}" == "0.000000" ]]; then
    echo "[$label] ERROR: Output appears empty (max=0). Keeping original."; rm -f "$tmp"; return 1
  fi
  p1=$(fslval "$tmp" pixdim1); p2=$(fslval "$tmp" pixdim2); p3=$(fslval "$tmp" pixdim3)
  awk -v a="$p1" -v b="$p2" -v c="$p3" 'BEGIN{
    dx=(a-2.7); dy=(b-2.7); dz=(c-2.97);
    if (dx<0) dx=-dx; if (dy<0) dy=-dy; if (dz<0) dz=-dz;
    if (dx>1e-3 || dy>1e-3 || dz>1e-3) exit 1; else exit 0;
  }' || echo "[$label] WARNING: pixdims are ($p1,$p2,$p3), expected (~2.7,2.7,2.97)."

  mv -f "$tmp" "$src4d"
  echo "[$label] Success → $src4d"
}

# ---------- T1w ----------
ref_t1w="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-T1w_boldref.nii.gz"
src_t1w="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-T1w_desc-preproc_bold.nii.gz"
regrid_identity "$ref_t1w" "$src_t1w" "T1w"

# ---------- MNI152NLin6Asym (no res-2) ----------
ref_mni="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-1_part-mag_space-MNI152NLin6Asym_boldref.nii.gz"
src_mni="${fmriprepdir}/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
regrid_identity "$ref_mni" "$src_mni" "MNI"
