#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <subject_id>   e.g. $0 101"
  exit 1
fi

SUBJ="$1"

for ses in $(seq -w 1 12); do
  SRC="/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/anat-only/ses-${ses}/sub-${SUBJ}"
  DST="/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/fmriprep/sub-${SUBJ}/ses-${ses}"

  # Skip entire session if the main source folder is missing
  if [ ! -d "${SRC}" ]; then
    echo "${SRC} not found, skipping ses-${ses}" >&2
    continue
  fi

  mkdir -p "${DST}"

  # 1) Copy everything except the nested ses-${ses} folder
  rsync -av --exclude="ses-${ses}/" "${SRC}/" "${DST}/"

  # 2) If the nested ses-${ses} folder exists, copy its contents too
  if [ -d "${SRC}/ses-${ses}" ]; then
    rsync -av "${SRC}/ses-${ses}/" "${DST}/"
  else
    echo "  (inner folder ${SRC}/ses-${ses} missing — skipped)" >&2
  fi
done
