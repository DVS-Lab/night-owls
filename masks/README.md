# Masks used in this project

We now keep two **canonical** mask products in this folder that match the study’s MNI grid (the same grid used by our fMRIPrep `space-MNI152NLin6Asym` boldrefs). The original “2 mm” source files remain here for provenance.

## What changed

- The source masks (`*_2mm.nii`) shipped in **MNI152NLin6Asym, res-2**.  
- We **resampled once** onto our study’s **MNI152NLin6Asym** grid (no `res-2` tag), using a 3-D MNI boldref as the reference.  
- Interpolation matched data type:
  - **Binary VS-Imanova** → `NearestNeighbor`
  - **Continuous BrainRewardSignature** → `Linear`
- Outputs are written in this directory with BIDS-like names:
  - `space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz`
  - `space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz`

> Expected voxel sizes on our study MNI grid: ~**2.7 × 2.7 × 2.97 mm**. The grid comes from the dataset’s MNI boldrefs.

## Files in this folder

**Original sources (kept unchanged)**  
- `VS-Imanova_2mm.nii` — binary mask from the Oxford-GSK-Imanova striatum atlas (see reference).  
- `BrainRewardSignature_2mm.nii` — continuous map from Speer et al. (see reference).

**Resampled study-ready products**  
- `space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz` — NN-resampled binary mask on study MNI grid.  
- `space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz` — linearly resampled continuous map on study MNI grid.

## How we generated the study-grid masks (reproducible note)

From this `masks/` directory, we used a one-off script that resamples to the study MNI grid defined by a 3-D MNI boldref:

```bash
# usage
bash resample_to_study_mni_grid.sh /full/path/to/sub-XXX_ses-YY_task-ZZ_run-1_part-mag_space-MNI152NLin6Asym_boldref.nii.gz
# or allow the script to auto-pick the first MNI boldref it finds in derivatives/
