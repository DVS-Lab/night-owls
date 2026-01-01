#!/usr/bin/env python3
"""
Extract FIR responses from FEAT cope images (ROI mean) and plot at:
- run level (sub/ses/run)
- session level (sub/ses; averaged across runs)
- subject level (sub; averaged across sessions and runs)

Inputs (hard-coded):
- A text file listing either:
    * stats/cope1.nii.gz paths, OR
    * .feat directories
  One per line, comments (#) allowed.

Notes:
- Uses design.con contrast names when present (preferred; works for MID too).
- Falls back to a sharedreward-specific cope mapping if design.con is missing.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt


# --------------------------- Config ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent

FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
MASKS_DIR  = ROOT_DIR / "masks"

# Mask (same as before)
VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"

# List of FEATs to process: can contain cope1 paths OR .feat dirs
FEAT_LIST  = SCRIPT_DIR / "FIR_cope1_paths.txt"

# Outputs
OUT_DIR    = ROOT_DIR / "derivatives" / "extractions" / "fir_copes"
PLOTS_DIR  = OUT_DIR / "plots"


# --------------------------- Fallback mapping (sharedreward) ---------------------------
# Only used if a FEAT doesn't have design.con (rare).
# Matches your screenshot: 8 blocks * 10 bins = 80 copes.
def fallback_sharedreward_contrast_names() -> Dict[int, str]:
    names: Dict[int, str] = {}
    blocks = [
        "computer", "stranger",
        "pun_c", "neu_c", "rew_c",
        "pun_s", "neu_s", "rew_s",
    ]
    idx = 1
    for cond in blocks:
        for b in range(1, 11):
            names[idx] = f"{cond} ({b})"
            idx += 1
    return names


# --------------------------- Helpers ---------------------------
def die(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def read_feat_list(path: Path) -> List[Path]:
    if not path.exists():
        die(f"FEAT_LIST not found: {path}")

    feats: List[Path] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = Path(line)

        # Allow relative paths (relative to code/)
        if not p.is_absolute():
            p = (SCRIPT_DIR / p).resolve()

        # Accept either .feat dir or stats/cope1 path
        if p.name.endswith(".feat") and p.is_dir():
            feats.append(p)
        elif p.name == "cope1.nii.gz" and p.parent.name == "stats":
            feats.append(p.parent.parent)  # .../<FEAT>.feat
        else:
            # Try to coerce if someone passed ".../.feat/stats/cope1.nii.gz" with odd formatting
            s = p.as_posix()
            m = re.search(r"(.+?\.feat)/stats/cope1\.nii\.gz$", s)
            if m:
                feats.append(Path(m.group(1)))
            else:
                print(f"WARN: skipping unrecognized entry: {line}")

    feats = sorted(set(feats))
    if not feats:
        die("No usable FEAT entries found in FEAT_LIST.")
    return feats


def parse_design_contrast_names(feat: Path) -> Optional[Dict[int, str]]:
    con = feat / "design.con"
    if not con.exists():
        return None

    names: Dict[int, str] = {}
    # Lines look like: /ContrastName1 computer (1)
    pat = re.compile(r"^/ContrastName(\d+)\s+(.*)\s*$")
    for line in con.read_text(errors="ignore").splitlines():
        m = pat.match(line.strip())
        if m:
            k = int(m.group(1))
            v = m.group(2).strip()
            names[k] = v

    return names or None


def parse_name_to_condition_bin(name: str) -> Tuple[str, int]:
    """
    Expected: "<condition> (<bin>)"
    Example: "computer (1)" or "rew_s (10)"
    """
    m = re.match(r"^(.*?)\s*\((\d+)\)\s*$", name)
    if not m:
        # If naming differs, keep condition as whole string and set bin to -1
        return name.strip(), -1
    return m.group(1).strip(), int(m.group(2))


def parse_feat_metadata(feat: Path) -> Dict[str, str]:
    """
    Pull out fields from FEAT dir name like:
    L1_sub-105_ses-12_task-sharedreward_model-1_type-act_run-2_space-mni_single-echo_cnfds-fmriprep_unsmoothed_FIR.feat
    """
    meta = {
        "feat": feat.as_posix(),
        "sub": "",
        "ses": "",
        "task": "",
        "run": "",
        "echo": "unknown",
        "confounds": "unknown",
        "space": "unknown",
    }

    s = feat.name

    def grab(pat: str) -> str:
        mm = re.search(pat, s)
        return mm.group(1) if mm else ""

    meta["sub"]  = grab(r"sub-(\d+)")
    meta["ses"]  = grab(r"ses-(\d+)")
    meta["task"] = grab(r"task-([A-Za-z0-9]+)")
    meta["run"]  = grab(r"run-(\d+)")

    if "multi-echo" in s:
        meta["echo"] = "multi-echo"
    elif "single-echo" in s:
        meta["echo"] = "single-echo"

    cf = grab(r"(?:cnfds|confounds)-([A-Za-z0-9]+)")
    if cf:
        meta["confounds"] = cf

    sp = grab(r"space-([A-Za-z0-9]+)")
    if sp:
        meta["space"] = sp

    return meta


def load_mask_bool(mask_path: Path) -> np.ndarray:
    if not mask_path.exists():
        die(f"Mask not found: {mask_path}")
    mimg = nib.load(str(mask_path))
    mdat = np.asanyarray(mimg.dataobj)
    return (mdat > 0) & np.isfinite(mdat)


def roi_mean_3d(img_path: Path, mask_bool: np.ndarray) -> float:
    img = nib.load(str(img_path))
    dat = np.asanyarray(img.dataobj)  # memmap when possible
    if dat.shape[:3] != mask_bool.shape:
        raise ValueError(f"Mask shape {mask_bool.shape} != image shape {dat.shape[:3]} for {img_path}")
    vals = dat[mask_bool]
    # coerce to float
    vals = vals.astype(np.float64, copy=False)
    if vals.size == 0:
        return float("nan")
    return float(np.nanmean(vals))


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def plot_fir(df: pd.DataFrame, out_png: Path, title: str) -> None:
    """
    df expected columns: condition, bin, value
    Plots condition-wise FIR curves.
    """
    ensure_dir(out_png.parent)

    # drop invalid bins if present
    d = df[df["bin"] > 0].copy()
    if d.empty:
        return

    plt.figure(figsize=(10, 6))
    for cond, g in d.groupby("condition", sort=True):
        g = g.sort_values("bin")
        plt.plot(g["bin"].to_numpy(), g["value"].to_numpy(), marker="o", label=cond)

    plt.axhline(0, linewidth=1)
    plt.xlabel("FIR bin (1–10)")
    plt.ylabel("ROI mean (cope units)")
    plt.title(title)
    plt.legend(fontsize=9, ncols=2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


# --------------------------- Main ---------------------------
def main() -> None:
    feats = read_feat_list(FEAT_LIST)
    print(f"Found {len(feats)} FEAT dirs in list.")

    mask_bool = load_mask_bool(VS_MNI)

    rows: List[Dict[str, object]] = []
    missing_imgs: List[str] = []

    fb_names = fallback_sharedreward_contrast_names()

    for feat in feats:
        stats_dir = feat / "stats"
        if not stats_dir.exists():
            missing_imgs.append(f"{feat} :: missing stats/")
            continue

        # Prefer design.con names (works for MID too)
        names = parse_design_contrast_names(feat)

        # If no design.con, fall back *only* if it's sharedreward-ish (or still try)
        if names is None:
            # Use fallback when task is sharedreward, otherwise we can't label bins reliably
            meta_tmp = parse_feat_metadata(feat)
            if meta_tmp["task"] == "sharedreward":
                names = fb_names
            else:
                # Still attempt to extract cope1..80 without names; label as copeN
                names = {i: f"cope{i}" for i in range(1, 81)}

        meta = parse_feat_metadata(feat)

        # Determine how many copes to expect based on the names we found
        max_cope = max(names.keys())
        # If someone has extra entries, we still only read those with actual files

        for cope_num in range(1, max_cope + 1):
            img_p = stats_dir / f"cope{cope_num}.nii.gz"
            if not img_p.exists():
                missing_imgs.append(str(img_p))
                continue

            name = names.get(cope_num, f"cope{cope_num}")
            cond, bin_ = parse_name_to_condition_bin(name)

            try:
                val = roi_mean_3d(img_p, mask_bool)
            except Exception as e:
                missing_imgs.append(f"{img_p} :: ERROR {e}")
                continue

            rows.append({
                **meta,
                "cope": cope_num,
                "cope_name": name,
                "condition": cond,
                "bin": bin_,
                "value": val,
            })

    ensure_dir(OUT_DIR)
    ensure_dir(PLOTS_DIR)

    df = pd.DataFrame(rows)
    out_csv = OUT_DIR / "fir_copes_roi_means_long.csv"
    df.to_csv(out_csv, index=False)
    print(f"Wrote: {out_csv}  (rows={len(df)})")

    if missing_imgs:
        miss_txt = OUT_DIR / "missing_or_problem_files.txt"
        miss_txt.write_text("\n".join(missing_imgs) + "\n")
        print(f"Wrote missing/problem list: {miss_txt}  (n={len(missing_imgs)})")
    else:
        print("No missing/problem files detected during extraction.")

    if df.empty:
        die("No data extracted; cannot plot.")

    # ---------------- Plots ----------------
    # Run-level
    run_cols = ["sub", "ses", "task", "run", "echo", "confounds"]
    for keys, g in df.groupby(run_cols, dropna=False):
        sub, ses, task, run, echo, conf = keys
        title = f"{task} | sub-{sub} ses-{ses} run-{run} | {echo} | {conf}"
        out_png = PLOTS_DIR / "run" / f"task-{task}" / f"sub-{sub}" / f"ses-{ses}" / f"run-{run}_{echo}_{conf}.png"
        plot_fir(g[["condition", "bin", "value"]], out_png, title)

    # Session-level (average across runs within session)
    ses_cols = ["sub", "ses", "task", "echo", "confounds"]
    df_ses = (
        df[df["bin"] > 0]
        .groupby(ses_cols + ["condition", "bin"], as_index=False)["value"]
        .mean()
    )
    for keys, g in df_ses.groupby(ses_cols, dropna=False):
        sub, ses, task, echo, conf = keys
        title = f"{task} | sub-{sub} ses-{ses} (avg runs) | {echo} | {conf}"
        out_png = PLOTS_DIR / "session" / f"task-{task}" / f"sub-{sub}" / f"ses-{ses}_{echo}_{conf}.png"
        plot_fir(g[["condition", "bin", "value"]], out_png, title)

    # Subject-level (average across sessions and runs)
    sub_cols = ["sub", "task", "echo", "confounds"]
    df_sub = (
        df[df["bin"] > 0]
        .groupby(sub_cols + ["condition", "bin"], as_index=False)["value"]
        .mean()
    )
    for keys, g in df_sub.groupby(sub_cols, dropna=False):
        sub, task, echo, conf = keys
        title = f"{task} | sub-{sub} (avg ses+runs) | {echo} | {conf}"
        out_png = PLOTS_DIR / "subject" / f"task-{task}" / f"sub-{sub}_{echo}_{conf}.png"
        plot_fir(g[["condition", "bin", "value"]], out_png, title)

    print("Done. Plots are under:", PLOTS_DIR)


if __name__ == "__main__":
    main()
