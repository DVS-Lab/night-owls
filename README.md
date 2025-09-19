# night-owls
Repository for the Night Owls Scan Club (NOSC) Project under Temple University IRB #24452. NOSC is a mutliband, multi-echo intensively sampled fMRI study of the reward response. 

NOSC is described in detail in Mattoni et al., 2025 (cite).
[Study design summary](stimuli/StudyDesign.png)

For data use please cite Mattoni et al., 2025 and (open neuro link)

BIDS data are publicly available at:
Information of scanning sessions, behavioral data, and outputs of L1 and LSS models described in Mattoni et al., 2025 are available on OSF at: 

## Preprocessing

Data were preprocessed in an HPC environment. All preprocessing code is in /code. 

Before upload to OpenNeuro, raw DICOMS were BIDS-fied using prepdata.sh and events.tsv files were generated using events_generation.R and convertSharedReward2BIDSevents.m. 

### Field map generation 

warpkit-hpc.sh generates fieldmap files in /bids/sub-xx/ses-xx from multi-echo data. 
addIntendedFor_fieldmap-hpc.py edits json files to include IntendedFor fields for the generated fieldmap files. 

### FMRIprep

fmriprep was run in a 2-step process to create a single anatomical image per subject and avoid processing multiple sessions in parallel. 

fmriprep-hpc-anat.sh performs anatomical only fmriprep preprocessing using --anat-only and --longitudinal arguments, creating one T1w image for all sessions for each subject. 
gen_fmriprep-anat.sh creates functional fmriprep commands specific to each session, taking the preprocessed anatomical data as an existing derivative. 
run_fmriprep_qsub submits all fmriprep commands created in /code/fmriprep-anat/
fmriprepOrganize.sh puts all generated fmriprep output in BIDS format and removes intermediate files generated in /derivatives/anat-only/

### Tedana
tedana-hpc.sh estimates tedana confounds for fmriprepped data

## FSL Analyses

gen3colfiles.sh turns events.tsv files into FSL-compatible events files
MakeConfounds.py adds fmriprep confounds to /derivatives/fsl/confounds_tedana
genTedanaMultiSes.py adds selected tedana confounds to /derivatives/fsl/confounds_tedana

L1-stats-loop.sh estimates L1 models for each run using different combinations of: space (MNI vs T1w), echo (echo-2 vs. multi), confounds (base fmriprep vs fmripre + tedana)
L1statsSingleTrial-${task}.sh runs LSS models for the respective task (mid or sharedreward)

extractData.sh and extractData-LSS.sh return derivatives used in Mattoni et al., 2025. 


