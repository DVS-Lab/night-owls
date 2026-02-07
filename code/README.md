# READ ME

We largely refer readers to thr read me in the parent directory above, which describes the core processing code steps. The Mattoni2025 subdirectory also contains R code for analysis in the main night owls manuscript. All information below was mostly AI-generated and describes additional files generated for various supplemental and exploratory processing steps and analyses. It should be interpretted with substantial caution. 

```
<project-root>/
  bids/                 # BIDS-formatted dataset
  derivatives/          # fMRIPrep, MRIQC, FSL FEAT outputs, extractions, etc.
  logs/                 # job logs / command logs
  masks/                # ROI masks used for extraction (e.g., VS/NAcc)
  code/                 # (this folder) scripts in this README
```

A large fraction of the scripts are **HPC/PBS job scripts** with hard-coded paths (e.g., `/gpfs/scratch/.../night-owls`). If you copy these scripts to a new environment, your first step is to update the `maindir`/`projectdir` variables and (if needed) module loads.

---

## “What do I run, and in what order?”

A typical end-to-end workflow looks like this (not every step is always required):

1. **Get raw data into BIDS**
   - `prepdata.sh` (single subject/session) or `run_prepdata.sh` (batch)
   - (Optional fixes) `addIntendedFor_fieldmap-hpc.py`, `addUnits_func.py`, `shiftdates.py`

2. **Quality control**
   - `mriqc-hpc.sh` + `run_mriqc-hpc.sh`
   - (Optional metrics export) `gen_MRIQC-outputs.py`

3. **Preprocess with fMRIPrep**
   - **Anat-only, session-split workflow**:
     - `gen_fmriprep_anat.sh` (generates per-subject/session scripts in `code/fmriprep-anat/`)
     - `run_fmriprep-qsub.sh` (submits the generated scripts)
     - `fmriprepOrganize.sh` (rsync/organize outputs into the final `derivatives/fmriprep/` layout)
   - **Alternative “single job script” approach**:
     - `fmriprep-hpc-anat.sh` + `run_fmriprep-hpc-anat.sh`

4. ** Distortion correction with WarpKit**
   - `warpkit-hpc.sh` + `run_warpkit-hpc.sh`

5. **Multi-echo denoising with tedana (multi-echo runs)**
   - `tedana-hpc.sh` + `run_tedana-hpc.sh`
   - Then generate FSL-ready confounds:
     - `MakeConfounds.py` (baseline fMRIPrep confounds → TSV)
     - `genTedanaMultiSes.py` (adds tedana rejected components → TSV)

6. **Generate FSL EV files from BIDS events**
   - `gen3colfiles.sh` (single subject/session) or `run_gen3colfiles.sh` (batch)
   - Uses `BIDSto3col.sh` (Tom Nichols’ converter; copied into this repo)

7. **Run FSL FEAT statistics**
   - Level 1 (run-level): `L1stats*` scripts (standard, FIR, FLOBS, LSS/single-trial)
   - Level 2 (subject/session aggregations): `L2stats-hpc*.sh`
   - Level 3 (group-level tests): `L3-ttest.sh`, `L3-trend.sh` (submitted via `run_L3stats.sh`)

8. **Extraction, QC, and plotting**
   - FIR QC: `check_FIR_outputs.py`
   - ROI extraction: `extract_fir_roi_means.py`, `extract_event_metrics_raw_and_fir.py`
   - Plotting: `plot_fir_timecourses.py`
   - Similarity/reliability: `imageSimilarity.py` + `run_imageSimilarity.sh`
   - AFNI smoothness estimation: `extractSmoothness.sh`

---

## Conventions used across scripts

### Subject/session lists
- `sublist.txt` — one subject ID per line (no `sub-` prefix)
- `sublist-ses.txt` — two columns: `sub  ses` (e.g., `101 01`)

### Echo / acquisition naming
- **Single-echo** runs use `echo-2` images from fMRIPrep (`*_echo-2_*desc-preproc_bold.nii.gz`).
- **Multi-echo** runs typically use the optimally combined magnitude image (`*_part-mag_*desc-preproc_bold.nii.gz`), and `tedana` uses the individual echo series (`echo-1..4`).

### Confounds naming
You will see two parallel confound options:
- `cnfds-fmriprep` — confounds derived from fMRIPrep only
- `cnfds-tedana` — fMRIPrep confounds + tedana rejected component time series

### Output naming
FSL FEAT output directories embed analysis settings in the folder name, for example:
- `.../L1_sub-101_ses-01_task-mid_model-1_type-act_run-1_space-mni_multi-echo_cnfds-tedana.feat/`

---

## Script reference (grouped)

### 1) Data ingestion and BIDS hygiene
- `downloadXNAT.py` — downloads DICOMs from an XNAT instance (intended for local pre-BIDS staging).
- `prepdata.sh` — “prepare data” wrapper (local/Apptainer):
  1) `heudiconv` DICOM→BIDS using `heuristics_XA30.py`
  2) deface T1w with `pydeface`
  3) shift acquisition dates in `*_scans.tsv` via `shiftdates.py`
  4) (project-specific) moves localizer DICOM folders out of the way to avoid slice-index errors during conversion
- `run_prepdata.sh` — runs `prepdata.sh` in parallel over pairs in `sublist-ses.txt`.
- `heuristics_XA30.py` — HeuDiConv heuristic that maps scan names to BIDS outputs (anat/func/fmap).
- `addIntendedFor_fieldmap-hpc.py` — adds/repairs `IntendedFor` fields in fieldmap JSONs.
- `addUnits_func.py` — sets `Units` in functional JSONs (note: uses hard-coded paths; adjust before use).
- `shiftdates.py` — shifts acquisition timestamps in `*_scans.tsv` (de-identification support).
- `make_zero_duration.sh` — sets EV duration column to zero for FIR EV files (useful when your FIR model expects impulse events).

**One-off data fixes (project-specific):**
- `sub-103_ses-02_replaceRun2.sh` — rename/move run-3 files to run-2 when a run label mismatch occurred.
- `resample_sub-101_ses-03_task-mid.sh` — regrids one run to match another run’s reference grid (identity transform).
- `normEcho2.sh` — normalizes/resamples echo-2 images and reference grids (used to repair inconsistent headers/grids).

---

### 2) QC and preprocessing (MRIQC / fMRIPrep / WarpKit)
**MRIQC**
- `mriqc-hpc.sh` — PBS job script to run MRIQC in a container.
- `run_mriqc-hpc.sh` — submits `mriqc-hpc.sh` jobs across subjects from `sublist.txt`.
- `gen_MRIQC-outputs.py` — parses MRIQC `*_bold.json` outputs into a single CSV of summary metrics.

**fMRIPrep**
- `gen_fmriprep_anat.sh` — generates per-(subject, session) PBS scripts under `code/fmriprep-anat/` and matching JSON filters.
- `run_fmriprep-qsub.sh` — submits the generated scripts in `code/fmriprep-anat/`.
- `fmriprep-hpc-anat.sh` / `run_fmriprep-hpc-anat.sh` — alternative “one script submits many subjects” workflow.
- `fmriprepOrganize.sh` — organizes/rsyncs outputs from an intermediate `derivatives/anat-only/` staging layout into `derivatives/fmriprep/`.

**WarpKit**
- `warpkit-hpc.sh` — PBS job script to run WarpKit-based distortion correction for (sub, ses) pairs.
- `run_warpkit-hpc.sh` — reads `sublist-ses.txt` and submits a single WarpKit job carrying all pairs.

---

### 3) Multi-echo denoising and confounds
- `tedana-hpc.sh` — runs `tedana` across tasks/runs and writes job-specific command logs.
- `run_tedana-hpc.sh` — submits tedana jobs in chunks across `sublist.txt`.
- `genTedanaMultiSes.py` — creates `desc-TedanaPlusConfounds.tsv` by appending tedana rejected-component regressors to a confounds table.
- `MakeConfounds.py` — transforms fMRIPrep’s confounds into an FSL-friendly TSV (drops/renames columns, adds derivatives, etc.).

---

### 4) Events → EV files (FSL 3-column)
- `BIDSto3col.sh` — Tom Nichols’ script to convert a BIDS `*_events.tsv` file to FSL-style 3-column EV files.
- `gen3colfiles.sh` — calls `BIDSto3col.sh` per subject/session, for both tasks and runs; writes to `derivatives/fsl/EVFiles/...`.
- `run_gen3colfiles.sh` — example batch runner for `gen3colfiles.sh`.

**MATLAB/R event utilities**
- `convertSharedReward2BIDSevents.m` — converts task-specific behavioral logs into BIDS events (shared reward task).
- `MakeSingleTrialsEV_MID.m`, `MakeSingleTrialsEV_SR.m` — creates single-trial EV definitions for LSS analyses.
- `events_generation.R` — R-based event generation utilities (legacy/alternative to the MATLAB path).

---

### 5) FSL modeling: Level 1 (run-level)
There are multiple first-level “families.” The naming is consistent:

- **Standard GLM** (canonical HRF):
  - `L1stats-loop.sh` + `L1stats.qsub` + `run_L1stats-qsub.sh`

- **FLOBS basis set (exploratory)**:
  - `L1stats_FLOBS.sh` + `L1stats_FLOBS.qsub` + `run_L1stats-qsub_FLOBS.sh`

- **FIR model**:
  - `L1stats_FIR.sh` + `L1stats_FIR.qsub` + `run_L1stats-qsub_FIR.sh`
  - `check_FIR_outputs.py` helps verify expected FIR outputs are present.

- **Single-trial / LSS (Least-Squares Separate)**:
  - `L1statsSingleTrial-mid.sh`, `L1statsSingleTrial-sharedreward.sh`
  - FLOBS variants: `L1statsSingleTrial-mid_FLOBS.sh`, `L1statsSingleTrial-sharedreward_FLOBS.sh`
  - PBS drivers:
    - `L1stats_LSS.qsub`, `L1stats_LSS-FLOBS.qsub`
    - `run_L1stats-qsub_LSS.sh`, `run_L1stats-qsub_LSS-FLOBS.sh`
  - Utility: `count_lss_completion.sh` (counts completed single-trial FEAT directories)

- **T1w-space (smoothed) variants**:
  - `L1stats-t1w.qsub` and `L1stats-mid-t1w.sh`, `L1stats-sharedreward-t1w.sh`

---

### 6) FSL modeling: Level 2 and Level 3
- `L2stats-hpc.sh` — session/subject aggregations over L1 FEAT outputs (supports multiple tasks/spaces/echoes/confounds).
- `L2stats-hpc-subj.sh` — subject-level summaries (typically across sessions).
- `L2stats-hpc-subj-qc.sh` — QC-focused subject-level summaries.
- `L2stats-hpc-t1w.sh` — T1w-space variant.
- `run_L2stats-hpc.sh` — submits L2 scripts in chunks across `sublist.txt`.

Group-level:
- `L3-ttest.sh` — group-level t-tests (configured for a subset of options in-script).
- `L3-trend.sh` — group-level trend tests (e.g., longitudinal effects).
- `run_L3stats.sh` — submits L3 scripts across `sublist.txt`.

---

### 7) Smoothing and smoothness estimation
- `smooth-3dBlurToFWHM.sh` — applies AFNI `3dBlurToFWHM` smoothing (default 5 mm) to fMRIPrep outputs.
- `smoothing.qsub` + `run_smoothing-qsub.sh` — HPC drivers to smooth across subjects/sessions/tasks/runs/spaces.
- `extractSmoothness.sh` — estimates effective smoothness (ACF-based) using AFNI `3dFWHMx` and writes a TSV.

---

### 8) Extraction, time courses, and plotting
**FIR extraction / plotting**
- `extract_fir_roi_means.py` — extracts ROI means from FIR `cope*.nii.gz` images and writes `derivatives/extractions/fir/fir_copes_roi_means_long.csv`.
- `plot_fir_timecourses.py` — plots average FIR timecourses from the extracted CSV and saves figures to `derivatives/extractions/fir/plots/`.
- `extract_event_metrics_raw_and_fir.py` — computes summary metrics (AUC/peaks/contrasts) for raw ROI timecourses and FIR timecourses.

**“Wu-normalized” raw time courses (legacy + newer unified scripts)**
- `vs_timecourse-mid.py`, `vs_timecourse-sr.py`, `vs_timecourse_FIR.py` — earlier task-specific pipelines.
- `vs_timecourse_unified.py`, `vs_timecourse_unified_wu.py` — unified extraction scripts that read a list of FEAT directories and output long-form timecourse tables and QC plots.

**Image similarity / reliability**
- `imageSimilarity.py` — computes pairwise similarity (Spearman) across images (e.g., session-to-session) for specific contrasts.
- `run_imageSimilarity.sh` — PBS wrapper to run `imageSimilarity.py` on the cluster.

---

## Notes for new users of this code base

- **Hard-coded paths are common.** Search for `maindir=`, `projectdir=`, `bidsdir=`, and `/gpfs/` and update to your environment.
- **HPC scheduler assumptions:** many scripts use PBS directives (`#PBS ...`) and `qsub`. If your cluster uses SLURM, you will need to translate these job headers.
- **Software dependencies:** scripts assume availability of **FSL**, **AFNI** (for `3dBlurToFWHM`, `3dFWHMx`), **ANTs** (for `antsApplyTransforms`), and container runtimes (**Singularity/Apptainer**). Python scripts typically require `numpy`, `pandas`, `nibabel`, `matplotlib`, and (for similarity) `pyrelimri`.
- **Run location matters:** many scripts infer the project root as the parent of `code/`. Run them from `code/` unless the script explicitly states otherwise.

