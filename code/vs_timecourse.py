#!/usr/bin/env python3
"""
VS time-point extractor (discrete 4th-TR after event onset)
-----------------------------------------------------------

What it does
============
- Reads a list of FEAT directories from `code/feat_paths.txt` (one absolute path per line; `#` comments ok).
- For each FEAT/run, loads the VS ROI mean time series from `filtered_func_data.nii.gz`.
- Builds two normalizations of the run-wise VS time series:
    * PSC: ((ts / run_mean) - 1) * 100
    * Z:   (ts - run_mean) / run_sd
- Without any interpolation, extracts the value at the **start of the 4th TR window** after event onset:
    * With TR=1.615, this targets the interval [onset + 3*TR, onset + 4*TR) = [4.845 s, 6.460 s) after onset.
- Aggregates the **mean across trials** at that discrete time point for six conditions:
    * Anticipation: Reward, Neutral
    * Feedback: Positive-Reward, Negative-Reward, Positive-Neutral, Negative-Neutral
- Writes a per-run spreadsheet with PSC and Z values for all six conditions + trial counts.

Outputs
=======
`derivatives/extractions/summary_at_4thTR.tsv` and `.csv`
    Columns: sub, ses, run, echo,
             ANT_REWARD_PSC, ANT_NEUTRAL_PSC,
             FB_POS_REWARD_PSC, FB_NEG_REWARD_PSC, FB_POS_NEUTRAL_PSC, FB_NEG_NEUTRAL_PSC,
             ANT_REWARD_Z, ANT_NEUTRAL_Z,
             FB_POS_REWARD_Z, FB_NEG_REWARD_Z, FB_POS_NEUTRAL_Z, FB_NEG_NEUTRAL_Z,
             N_ANT_REWARD, N_ANT_NEUTRAL, N_FB_POS_REWARD, N_FB_NEG_REWARD, N_FB_POS_NEUTRAL, N_FB_NEG_NEUTRAL

Assumptions
===========
- All data are in MNI space. No transforms are performed here.
- EV files live at: <rootdir>/derivatives/fsl/EVFiles/sub-<sub>/ses-<ses>/mid/run-<run>/
  with filenames: `_anticipation_reward.txt`, `_anticipation_neutral.txt`,
                  `_feedback_positive_reward.txt`, `_feedback_negative_reward.txt`,
                  `_feedback_positive_neutral.txt`, `_feedback_negative_neutral.txt`.
- VS mask: <rootdir>/masks/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz
- TR is 1.615 s (hard-coded) with a header sanity-check.

Run from <rootdir>/code:
    python vs_tp_extract.py
"""
from __future__ import annotations

import re
import math
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import nibabel as nib

# ---------------- Fixed locations (relative to THIS script) -------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOTDIR = SCRIPT_DIR.parent
FSL_DERIV = ROOTDIR / "derivatives" / "fsl"
MASKS_DIR = ROOTDIR / "masks"
OUT_DIR   = ROOTDIR / "derivatives" / "extractions"
FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths.txt"

VS_MNI = MASKS_DIR / "space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

# Timing
TR = 1.615  # seconds, hard-coded; we sanity-check against header values
K_AFTER = 3  # 4th TR after onset => use start of [onset + 3*TR, onset + 4*TR)

# -------------------------- Utility functions --------------------------------

def parse_ids_from_feat(feat: Path) -> Tuple[str, str, str, str]:
    """Return (sub, ses, run, echo) parsed from the FEAT directory name/path.
    Expects names like: L1_sub-101_ses-01_task-mid_model-1_type-act_run-2_space-mni_multi-echo_cnfds-fmriprep+.feat
    """
    m = re.search(r"sub-(\d+).*?ses-(\d+).*?run-(\d+).*?(single-echo|multi-echo)", feat.as_posix())
    if not m:
        raise ValueError(f"Cannot parse sub/ses/run/echo from: {feat}")
    sub, ses, run, echo = m.group(1), m.group(2), m.group(3), m.group(4)
    return sub, ses, run, echo


def load_ev(ev_path: Path) -> np.ndarray:
    if not ev_path.exists():
        return np.array([], dtype=float)
    arr = []
    for line in ev_path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        parts = s.split()
        if not parts:
            continue
        try:
            arr.append(float(parts[0]))  # onset is first column
        except ValueError:
            continue
    return np.array(arr, dtype=float)


def fourth_tr_indices(onsets: np.ndarray, tr: float, n_vols: int) -> np.ndarray:
    """Indices of the volume whose start lies in [onset + 3*TR, onset + 4*TR).
    No interpolation; indexes are dropped if out-of-range.
    """
    if onsets.size == 0:
        return np.empty((0,), dtype=int)
    idx = np.floor((onsets + K_AFTER * tr) / tr).astype(int)
    return idx[(idx >= 0) & (idx < n_vols)]


def roi_meants_4d(img_4d: nib.Nifti1Image, mask_img: nib.Nifti1Image) -> np.ndarray:
    """Return ROI mean time series (float64) from a 4D fMRI image using a boolean mask.
    Assumes images are already in the same space/resolution.
    """
    data = img_4d.get_fdata(dtype=np.float32)
    mask = mask_img.get_fdata().astype(bool)
    if data.shape[:3] != mask.shape:
        raise ValueError(f"Mask shape {mask.shape} != data shape {data.shape[:3]}")
    ts = data[mask].mean(axis=0)
    return ts.astype(np.float64)


# ------------------------------- Main logic ----------------------------------

def process_one_feat(feat_dir: Path, mask_img: nib.Nifti1Image) -> Dict[str, object] | None:
    sub, ses, run, echo = parse_ids_from_feat(feat_dir)

    func_path = feat_dir / "filtered_func_data.nii.gz"
    if not func_path.exists():
        print(f"[WARN] missing filtered_func_data: {func_path}")
        return None

    img = nib.load(str(func_path))
    hdr_tr = float(img.header.get_zooms()[3]) if img.ndim == 4 and len(img.header.get_zooms()) > 3 else None
    if hdr_tr is not None and abs(hdr_tr - TR) > 1e-3:
        print(f"[WARN] Header TR={hdr_tr:.6f} differs from hard-coded TR={TR:.6f}; using {TR:.6f} s")

    ts = roi_meants_4d(img, mask_img)  # length T
    T = ts.size

    # Run-wise PSC and Z
    run_mean = float(np.nanmean(ts))
    run_sd   = float(np.nanstd(ts, ddof=1)) if T > 1 else np.nan
    if run_mean == 0 or np.isnan(run_mean):
        print(f"[WARN] bad run mean in {func_path}")
        return None
    psc = (ts / run_mean - 1.0) * 100.0
    z   = (ts - run_mean) / run_sd if run_sd > 0 else np.full_like(ts, np.nan)

    # EV paths (hard-coded locations)
    ev_base = FSL_DERIV / "EVFiles" / f"sub-{sub}" / f"ses-{ses}" / "mid" / f"run-{run}"
    ev_ANT_R = ev_base / "_anticipation_reward.txt"
    ev_ANT_N = ev_base / "_anticipation_neutral.txt"

    ev_FB_PR = ev_base / "_feedback_positive_reward.txt"
    ev_FB_NR = ev_base / "_feedback_negative_reward.txt"
    ev_FB_PN = ev_base / "_feedback_positive_neutral.txt"
    ev_FB_NN = ev_base / "_feedback_negative_neutral.txt"

    # Load onsets
    ANT_R = load_ev(ev_ANT_R); ANT_N = load_ev(ev_ANT_N)
    FB_PR = load_ev(ev_FB_PR); FB_NR = load_ev(ev_FB_NR)
    FB_PN = load_ev(ev_FB_PN); FB_NN = load_ev(ev_FB_NN)

    # Indices at the 4th TR after onset (start of window)
    idx_ANT_R = fourth_tr_indices(ANT_R, TR, T)
    idx_ANT_N = fourth_tr_indices(ANT_N, TR, T)
    idx_FB_PR = fourth_tr_indices(FB_PR, TR, T)
    idx_FB_NR = fourth_tr_indices(FB_NR, TR, T)
    idx_FB_PN = fourth_tr_indices(FB_PN, TR, T)
    idx_FB_NN = fourth_tr_indices(FB_NN, TR, T)

    # Means at those discrete indices
    def m(arr: np.ndarray) -> float | None:
        return float(np.nanmean(arr)) if arr.size else None

    out = {
        "sub": sub, "ses": ses, "run": run, "echo": echo,
        # PSC
        "ANT_REWARD_PSC": m(psc[idx_ANT_R]),
        "ANT_NEUTRAL_PSC": m(psc[idx_ANT_N]),
        "FB_POS_REWARD_PSC": m(psc[idx_FB_PR]),
        "FB_NEG_REWARD_PSC": m(psc[idx_FB_NR]),
        "FB_POS_NEUTRAL_PSC": m(psc[idx_FB_PN]),
        "FB_NEG_NEUTRAL_PSC": m(psc[idx_FB_NN]),
        # Z
        "ANT_REWARD_Z": m(z[idx_ANT_R]),
        "ANT_NEUTRAL_Z": m(z[idx_ANT_N]),
        "FB_POS_REWARD_Z": m(z[idx_FB_PR]),
        "FB_NEG_REWARD_Z": m(z[idx_FB_NR]),
        "FB_POS_NEUTRAL_Z": m(z[idx_FB_PN]),
        "FB_NEG_NEUTRAL_Z": m(z[idx_FB_NN]),
        # counts
        "N_ANT_REWARD": int(idx_ANT_R.size),
        "N_ANT_NEUTRAL": int(idx_ANT_N.size),
        "N_FB_POS_REWARD": int(idx_FB_PR.size),
        "N_FB_NEG_REWARD": int(idx_FB_NR.size),
        "N_FB_POS_NEUTRAL": int(idx_FB_PN.size),
        "N_FB_NEG_NEUTRAL": int(idx_FB_NN.size),
    }
    return out


def main():
    if not FEAT_LIST_PATH.exists():
        print(f"[WARN] FEAT list not found: {FEAT_LIST_PATH}")
        print("Create feat_paths.txt with one FEAT directory per line (absolute paths).")
        return

    feat_paths: List[Path] = []
    for line in FEAT_LIST_PATH.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        feat_paths.append(Path(s))

    if not feat_paths:
        print("No FEATs to process (feat_paths.txt empty). Exiting.")
        return

    # Load mask once
    if not VS_MNI.exists():
        print(f"[ERROR] VS mask not found: {VS_MNI}")
        return
    mask_img = nib.load(str(VS_MNI))

    rows: List[Dict[str, object]] = []
    for feat in feat_paths:
        if not feat.exists():
            print(f"[WARN] FEAT missing: {feat}")
            continue
        try:
            row = process_one_feat(feat, mask_img)
            if row is not None:
                rows.append(row)
        except Exception as e:
            print(f"[ERROR] processing {feat}: {e}")

    if not rows:
        print("Nothing processed. Check feat_paths.txt and inputs.")
        return

    # Ensure output directory
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tsv = OUT_DIR / "summary_at_4thTR.tsv"
    csv = OUT_DIR / "summary_at_4thTR.csv"

    header = [
        "sub","ses","run","echo",
        "ANT_REWARD_PSC","ANT_NEUTRAL_PSC",
        "FB_POS_REWARD_PSC","FB_NEG_REWARD_PSC","FB_POS_NEUTRAL_PSC","FB_NEG_NEUTRAL_PSC",
        "ANT_REWARD_Z","ANT_NEUTRAL_Z",
        "FB_POS_REWARD_Z","FB_NEG_REWARD_Z","FB_POS_NEUTRAL_Z","FB_NEG_NEUTRAL_Z",
        "N_ANT_REWARD","N_ANT_NEUTRAL","N_FB_POS_REWARD","N_FB_NEG_REWARD","N_FB_POS_NEUTRAL","N_FB_NEG_NEUTRAL"
    ]

    def to_row(r: Dict[str, object]) -> List[str]:
        return [
            str(r.get("sub")), str(r.get("ses")), str(r.get("run")), str(r.get("echo")),
            _fmt(r.get("ANT_REWARD_PSC")), _fmt(r.get("ANT_NEUTRAL_PSC")),
            _fmt(r.get("FB_POS_REWARD_PSC")), _fmt(r.get("FB_NEG_REWARD_PSC")), _fmt(r.get("FB_POS_NEUTRAL_PSC")), _fmt(r.get("FB_NEG_NEUTRAL_PSC")),
            _fmt(r.get("ANT_REWARD_Z")), _fmt(r.get("ANT_NEUTRAL_Z")),
            _fmt(r.get("FB_POS_REWARD_Z")), _fmt(r.get("FB_NEG_REWARD_Z")), _fmt(r.get("FB_POS_NEUTRAL_Z")), _fmt(r.get("FB_NEG_NEUTRAL_Z")),
            str(r.get("N_ANT_REWARD")), str(r.get("N_ANT_NEUTRAL")), str(r.get("N_FB_POS_REWARD")), str(r.get("N_FB_NEG_REWARD")), str(r.get("N_FB_POS_NEUTRAL")), str(r.get("N_FB_NEG_NEUTRAL")),
        ]

    with open(tsv, 'w') as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(to_row(r)) + "\n")

    with open(csv, 'w') as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(to_row(r)) + "\n")

    print(f"Wrote: {tsv}\n       {csv}")


def _fmt(x: object) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    if isinstance(x, float):
        return f"{x:.6f}"
    return str(x)


if __name__ == "__main__":
    main()
