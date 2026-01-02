#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
import time
from pathlib import Path

import numpy as np
import nibabel as nib

# --------------------------- Config ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent

MASKS_DIR  = ROOT_DIR / "masks"
VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"

# Your cope1-path list (704 lines)
FEAT_LIST  = SCRIPT_DIR / "FIR_cope1_paths.txt"

OUT_DIR    = ROOT_DIR / "derivatives" / "extractions" / "fir"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV    = OUT_DIR / "fir_copes_roi_means_long.csv"

# FIR-relevant cope ranges by task
TASK_MAX_COPES = {
    "sharedreward": 80,
    "mid": 60,
}

PRINT_EVERY_FEATS = 10  # progress cadence


# ------------------------ Helpers -----------------------------
_RE_KV = {
    "sub": re.compile(r"sub-(\d+)"),
    "ses": re.compile(r"ses-(\d+)"),
    "task": re.compile(r"task-([A-Za-z0-9]+)"),
    "run": re.compile(r"run-(\d+)"),
    "space": re.compile(r"space-([A-Za-z0-9]+)"),
    "echo": re.compile(r"(single-echo|multi-echo)"),
    "confounds": re.compile(r"cnfds-([A-Za-z0-9]+)"),
}

# FEAT design.con often looks like:
#   /ContrastName1    ant_rew (1)
# sometimes it can be quoted:
#   /ContrastName1    "ant_rew (1)"
_RE_CONTRAST = re.compile(r"^/ContrastName(\d+)\s+(.+?)\s*$")
_RE_BIN = re.compile(r"^\s*(.+?)\s*\((\d+)\)\s*$")


def feat_dir_from_cope1(cope1_path: Path) -> Path:
    # .../.feat/stats/cope1.nii.gz  ->  .../.feat
    return cope1_path.parent.parent


def parse_feat_metadata(feat_dir: Path) -> dict:
    s = feat_dir.name
    out = {}
    for k, rgx in _RE_KV.items():
        m = rgx.search(s)
        out[k] = (m.group(1) if m else "")
    # strip leading zeros where relevant
    out["sub"] = out["sub"].lstrip("0") or out["sub"]
    out["ses"] = out["ses"].lstrip("0") or out["ses"]
    out["run"] = out["run"].lstrip("0") or out["run"]
    return out


def find_design_con(feat_dir: Path) -> Path:
    candidates = [feat_dir / "design.con", feat_dir / "stats" / "design.con"]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(f"No design.con found under {feat_dir}")


def read_contrast_names(design_con: Path) -> dict[int, str]:
    """
    Parse FEAT design.con ContrastName lines into {cope_index: cope_name}.
    Handles both quoted and unquoted names.
    """
    names: dict[int, str] = {}
    for raw in design_con.read_text(errors="ignore").splitlines():
        line = raw.strip()
        if not line.startswith("/ContrastName"):
            continue
        m = _RE_CONTRAST.match(line)
        if not m:
            continue
        idx = int(m.group(1))
        nm = m.group(2).strip()
        # remove surrounding quotes if present
        if (len(nm) >= 2) and ((nm[0] == nm[-1]) and nm[0] in ("'", '"')):
            nm = nm[1:-1].strip()
        names[idx] = nm

    if not names:
        raise ValueError(f"Could not parse any /ContrastName entries in {design_con}")
    return names


def load_mask_indices(mask_path: Path) -> tuple[tuple[int, int, int], np.ndarray]:
    if not mask_path.exists():
        raise FileNotFoundError(f"Missing mask: {mask_path}")
    mimg = nib.load(str(mask_path))
    mdat = np.asanyarray(mimg.dataobj)
    mbool = mdat > 0
    shape = mbool.shape
    idx = np.flatnonzero(mbool.ravel())
    if idx.size == 0:
        raise ValueError(f"Mask has 0 voxels > 0: {mask_path}")
    return shape, idx


def roi_mean_3d(img_path: Path, mask_shape: tuple[int, int, int], mask_idx: np.ndarray) -> float:
    img = nib.load(str(img_path))
    dat = np.asanyarray(img.dataobj)
    if dat.ndim != 3:
        raise ValueError(f"Expected 3D image, got shape {dat.shape} for {img_path}")
    if dat.shape != mask_shape:
        raise ValueError(f"Mask dims {mask_shape} != image dims {dat.shape} for {img_path}")

    vals = dat.ravel()[mask_idx]
    return float(np.nanmean(vals.astype(np.float32, copy=False)))


def parse_condition_and_bin(cope_name: str) -> tuple[str, int | None]:
    m = _RE_BIN.match(cope_name.strip())
    if not m:
        return cope_name.strip(), None
    return m.group(1).strip(), int(m.group(2))


# -------------------------- Main ------------------------------
def main() -> None:
    if not FEAT_LIST.exists():
        raise FileNotFoundError(f"Missing paths list: {FEAT_LIST}")

    lines = [
        ln.strip() for ln in FEAT_LIST.read_text().splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    ]
    cope1_paths = [Path(p) for p in lines]
    feat_dirs = [feat_dir_from_cope1(p) for p in cope1_paths]

    # de-dup while preserving order
    seen = set()
    feat_dirs = [f for f in feat_dirs if not (f in seen or seen.add(f))]

    print(f"Found {len(feat_dirs)} FEAT dirs in list.", flush=True)

    mask_shape, mask_idx = load_mask_indices(VS_MNI)

    header = [
        "feat", "sub", "ses", "task", "run", "echo", "confounds", "space",
        "cope", "cope_name", "condition", "bin", "value", "cope_path",
    ]

    t0 = time.time()
    n_feats = 0
    n_rows = 0
    n_skipped_task = 0
    n_warn_design = 0

    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()

        for feat_dir in feat_dirs:
            n_feats += 1
            meta = parse_feat_metadata(feat_dir)
            task = meta.get("task", "")

            if task not in TASK_MAX_COPES:
                n_skipped_task += 1
                if n_skipped_task <= 10:
                    print(f"WARN: Unknown task '{task}' in {feat_dir.name} (skipping).", flush=True)
                continue

            max_c = TASK_MAX_COPES[task]

            try:
                design_con = find_design_con(feat_dir)
                cope_names = read_contrast_names(design_con)
            except Exception as e:
                n_warn_design += 1
                print(f"WARN: {e} | skipping FEAT: {feat_dir}", flush=True)
                continue

            stats_dir = feat_dir / "stats"

            # iterate in numeric cope order, but only FIR-relevant range
            for cope_idx in sorted(i for i in cope_names.keys() if i <= max_c):
                cope_path = stats_dir / f"cope{cope_idx}.nii.gz"
                if not cope_path.exists():
                    continue

                cname = cope_names[cope_idx]
                cond, bin_ = parse_condition_and_bin(cname)

                val = roi_mean_3d(cope_path, mask_shape, mask_idx)

                w.writerow({
                    "feat": str(feat_dir),
                    "sub": meta["sub"],
                    "ses": meta["ses"],
                    "task": task,
                    "run": meta["run"],
                    "echo": meta["echo"],
                    "confounds": meta["confounds"],
                    "space": meta["space"],
                    "cope": cope_idx,
                    "cope_name": cname,
                    "condition": cond,
                    "bin": bin_,
                    "value": val,
                    "cope_path": str(cope_path),
                })
                n_rows += 1

            if (n_feats % PRINT_EVERY_FEATS) == 0:
                dt = time.time() - t0
                feats_per_min = n_feats / (dt / 60.0)
                rows_per_sec = n_rows / dt if dt > 0 else float("nan")
                print(f"[{n_feats:4d}/{len(feat_dirs)}] feats | {n_rows} rows | "
                      f"{feats_per_min:.2f} feats/min | {rows_per_sec:.1f} rows/s",
                      flush=True)

    print(f"\nWrote: {OUT_CSV}", flush=True)


if __name__ == "__main__":
    main()
