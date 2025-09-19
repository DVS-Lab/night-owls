# night-owls  

Repository for the Night Owls Scan Club (NOSC) Project under Temple University IRB #24452.  

NOSC is a multiband, multi-echo, intensively sampled fMRI study of the reward response.  

NOSC is described in detail in Mattoni et al., 2025 (cite).  

[Study design summary](stimuli/StudyDesign.png)  

For data use please cite Mattoni et al., 2025 and (open neuro link).  

BIDS data are publicly available at: (link).  

Information on scanning sessions, behavioral data, and outputs of L1 and LSS models (described in Mattoni et al., 2025) are available on OSF at: (link).  

---

## Preprocessing  

Data were preprocessed in an HPC environment.  

All preprocessing code is in `/code`.  

Before upload to OpenNeuro:  
- Raw DICOMS were BIDS-fied using `prepdata.sh`.  
- `events.tsv` files were generated using `events_generation.R` and `convertSharedReward2BIDSevents.m`.  

### Field map generation  

- `warpkit-hpc.sh` generates fieldmap files in `/bids/sub-xx/ses-xx` from multi-echo data.  
- `addIntendedFor_fieldmap-hpc.py` edits `.json` files to include `IntendedFor` fields for the generated fieldmap files.  

### FMRIprep  

fMRIPrep was run in a 2-step process to:  
- Create a single anatomical image per subject.  
- Avoid processing multiple sessions in parallel.  

Scripts:  
- `fmriprep-hpc-anat.sh` performs anatomical-only preprocessing (`--anat-only`, `--longitudinal`).  
  - Creates one T1w image for all sessions per subject.  
- `gen_fmriprep-anat.sh` creates functional fMRIPrep commands for each session.  
  - Uses preprocessed anatomical data as an existing derivative.  
- `run_fmriprep_qsub` submits all fMRIPrep commands created in `/code/fmriprep-anat/`.  
- `fmriprepOrganize.sh` organizes fMRIPrep output in BIDS format.  
  - Removes intermediate files from `/derivatives/anat-only/`.  

### Tedana  

- `tedana-hpc.sh` estimates tedana confounds for fMRIPrepped data.  

---

## FSL Analyses  

- `gen3colfiles.sh` converts `events.tsv` files into FSL-compatible events files.  
- `MakeConfounds.py` adds fMRIPrep confounds to `/derivatives/fsl/confounds_tedana`.  
- `genTedanaMultiSes.py` adds selected tedana confounds to `/derivatives/fsl/confounds_tedana`.  

Model estimation:  
- `L1-stats-loop.sh` estimates L1 models for each run using combinations of:  
  - Space: MNI vs T1w  
  - Echo: echo-2 vs multi  
  - Confounds: base fMRIPrep vs fMRIPrep + tedana  
- `L1statsSingleTrial-${task}.sh` runs LSS models for the respective task (mid or sharedreward).  

Data extraction:  
- `extractData.sh` and `extractData-LSS.sh` return derivatives used in Mattoni et al., 2025.  
