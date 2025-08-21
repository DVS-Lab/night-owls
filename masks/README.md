# Masks Used in This Project

This directory contains two functional masks used for ROI analyses. Both were originally in `space-MNI152NLin6Asym_res-2`, but have been resampled to match the project's standard MNI152NLin6Asym grid (`pixdim: 2.7 × 2.7 × 2.97` mm).

## Included Masks

### 1. `VS-Imanova_2mm.nii`

- **Description**: Binary mask of the ventral striatum (VS) derived from the Oxford-GSK-Imanova striatal structural atlas.
- **Original Source**: [FSL Atlases: Striatum](https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/Atlases%282f%29striatumstruc.html)
- **Citation**:
  > Tziortzi, A. C., Searle, G. E., Tzimopoulou, S., Salinas, C., Beaver, J. D., Jenkinson, M., Laruelle, M., Rabiner, E. A., & Gunn, R. N. (2011). Imaging dopamine receptors in humans with [11C]-(+)-PHNO: Dissection of D3 signal and anatomy. *NeuroImage, 54*(1), 264–277. https://doi.org/10.1016/j.neuroimage.2010.06.044

### 2. `BrainRewardSignature_2mm.nii`

- **Description**: Continuous-valued multivariate brain signature of reward, downloaded from the authors’ NeuroVault repository.
- **Download Source**: [NeuroVault Image #775976](https://neurovault.org/images/775976/)
- **Citation**:
  > Speer, S. P. H., Keysers, C., Barrios, J. C., Teurlings, C. J. S., Smidts, A., Boksem, M. A. S., Wager, T. D., & Gazzola, V. (2023). A multivariate brain signature for reward. *NeuroImage, 271*, 119990. https://doi.org/10.1016/j.neuroimage.2023.119990

## Notes on Resampling

Both masks have been resampled to match the project-standard MNI grid using `antsApplyTransforms`:

- **VS-Imanova mask**: Originally binary, resampled using *nearest neighbor* interpolation to preserve labeling.
- **BrainRewardSignature**: Contains continuous values and was resampled using *linear* interpolation.

> ✅ All resampling steps are documented in the script `resample_to_project_mni.sh` located in this directory.

---

*This README was written with assistance from ChatGPT-4o (OpenAI), used in "GPT-4o Thinking" mode for clarity and technical accuracy.*
