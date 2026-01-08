#!/usr/bin/env python3
"""
Extract FLOBS-based L1 metrics from FEAT outputs.

This is the L1 ("parent") analog of the LSS-FLOBS extractor you just ran.
Key behavior (per your notes):
  - We do NOT treat these as A>B contrasts.
  - For task=mid: use ONLY the first 18 contrasts in design.con (6 regressors * 3 FLOBS bases).
  - For task=sharedreward: use ONLY the first 24 contrasts in design.con (8 regressors * 3 FLOBS bases),
    even if some design.con files have more than 24 rows.
  - For each regressor-of-interest, we report metrics for:
      1) zstat for basis-1 (PE1-like) contrast
      2) scaled PE1 beta: pe[basis1_col] / PPheight[basis1_col]
      3) scaled signed RMS across bases 1..3:
           sign(PE1_scaled) * sqrt(PE1_scaled^2 + PE2_scaled^2 + PE3_scaled^2)
         where each PE is scaled by its own PPheight.

Outputs (TSV): one row per FEAT directory per regressor-of-interest (6 rows per MID FEAT, 8 per SR FEAT)
  - ROI means in NAcc mask and BRS cortical mask
  - Whole-brain correlation vs the BRS map for each metric (masked to FEAT brainmask)

Missing inputs are written as 'n/a' (e.g., missing stats/pe*, missing zstat*, missing design.mat, etc.)

Expected FEAT layout (example):
  derivatives/fsl/sub-101/ses-01/
    L1_sub-101_ses-01_task-mid_model-1_type-act_run-1_space-mni_multi-echo_cnfds-tedana_unsmoothed_FLOBS.feat/
      design.con
      design.mat
      mask.nii.gz
      stats/zstat1.nii.gz ...
      stats/pe1.nii.gz ...

Usage:
  python extractData_L1_FLOBS_EVmetrics.py
  python extractData_L1_FLOBS_EVmetrics.py --root /home/tug87422/scratch/smithlab-shared/night-owls

Notes:
  - Requires nibabel and numpy.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    import nibabel as nib
except ImportError as e:
    raise SystemExit(
        "ERROR: nibabel is required. In conda, try: conda install -c conda-forge nibabel"
    ) from e

NA = "n/a"


# ---------------------------
# small utilities
# ---------------------------

def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def safe_float(x: Optional[float]) -> str:
    if x is None or not np.isfinite(x):
        return NA
    return f"{x:.6g}"


def load_nii_float(path: Path) -> Optional[np.ndarray]:
    try:
        img = nib.load(str(path))
        return img.get_fdata(dtype=np.float32)
    except Exception:
        return None


def roi_mean(vol: Optional[np.ndarray], mask_bool: np.ndarray) -> Optional[float]:
    if vol is None:
        return None
    if vol.shape != mask_bool.shape:
        return None
    vals = vol[mask_bool]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return None
    return float(vals.mean())


def pearson_corr(a: Optional[np.ndarray], b: np.ndarray, mask_bool: np.ndarray) -> Optional[float]:
    """
    Correlation of a vs b within mask_bool, excluding non-finite.
    """
    if a is None:
        return None
    if a.shape != b.shape or a.shape != mask_bool.shape:
        return None
    x = a[mask_bool]
    y = b[mask_bool]
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if x.size < 3:
        return None
    x = x - x.mean()
    y = y - y.mean()
    den = np.sqrt((x * x).sum() * (y * y).sum())
    if den == 0:
        return None
    return float((x * y).sum() / den)


# ---------------------------
# parsing: design.mat scaling
# ---------------------------

def parse_design_mat_ppheights(design_mat: Path) -> Optional[List[float]]:
    """
    Parse /PPheights from a FEAT design.mat.

    Usually it appears on one long line, but for safety we support wrap-around
    across multiple subsequent lines until the next '/' header.

    Returns list of floats where 1-based design column j corresponds to ppheights[j-1].
    """
    try:
        lines = design_mat.read_text(errors="replace").splitlines()
    except Exception:
        return None

    # Optional: read /NumWaves to know the expected length
    num_waves = None
    for line in lines:
        if line.startswith("/NumWaves"):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    num_waves = int(parts[1])
                except ValueError:
                    num_waves = None
            break

    vals: List[float] = []
    in_pp = False
    for line in lines:
        if line.startswith("/PPheights"):
            in_pp = True
            parts = line.split()[1:]
        elif in_pp:
            if line.startswith("/"):
                break
            parts = line.split()
        else:
            continue

        for p in parts:
            try:
                vals.append(float(p))
            except ValueError:
                pass

        if num_waves is not None and len(vals) >= num_waves:
            vals = vals[:num_waves]
            break

    return vals if vals else None

    for line in txt:
        if line.startswith("/PPheights"):
            # split after the key; tolerate tabs/spaces
            parts = line.split()
            # parts[0] is /PPheights
            vals = []
            for p in parts[1:]:
                try:
                    vals.append(float(p))
                except ValueError:
                    pass
            return vals if vals else None
    return None


# ---------------------------
# parsing: design.con -> (base regressor -> basis1..3 contrast + PE column)
# ---------------------------

_CONTRAST_NAME_RE = re.compile(r"^/ContrastName(\d+)\s+(.*)$")
_BASIS_SUFFIX_RE = re.compile(r"^(.*)\s+\((\d+)\)\s*$")


def parse_design_con(design_con: Path) -> Tuple[Dict[int, str], Optional[np.ndarray]]:
    """
    Returns:
      names: {contrast_idx (1-based): name string}
      mat: np.ndarray of shape (NumContrasts, NumWaves) or None if parse fails
    """
    try:
        lines = design_con.read_text(errors="replace").splitlines()
    except Exception:
        return {}, None

    names: Dict[int, str] = {}
    for line in lines:
        m = _CONTRAST_NAME_RE.match(line)
        if m:
            idx = int(m.group(1))
            names[idx] = m.group(2).strip()

    # find /Matrix
    mat_start = None
    for i, line in enumerate(lines):
        if line.strip() == "/Matrix":
            mat_start = i + 1
            break
    if mat_start is None:
        return names, None

    rows = []
    for line in lines[mat_start:]:
        if not line.strip():
            continue
        parts = line.strip().split()
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            # stop if something non-numeric shows up
            break

    if not rows:
        return names, None
    mat = np.asarray(rows, dtype=float)
    return names, mat


@dataclass(frozen=True)
class RegressorTriplet:
    """
    One regressor-of-interest represented by a FLOBS basis triplet.

    contrast_idx_* refer to cope/zstat indices (row numbers in design.con)
    pe_col_* are the PE column indices (wave numbers in design matrix; 1-based)
    """
    base_name: str
    contrast_idx_1: int
    contrast_idx_2: int
    contrast_idx_3: int
    pe_col_1: int
    pe_col_2: int
    pe_col_3: int


def _row_to_single_col(mat: np.ndarray, row_1based: int, tol: float = 1e-8) -> Optional[int]:
    """
    If design.con row has exactly one nonzero entry, return that column index (1-based).
    Otherwise, return the column with max abs weight.
    """
    r = row_1based - 1
    if r < 0 or r >= mat.shape[0]:
        return None
    nz = np.where(np.abs(mat[r]) > tol)[0]
    if nz.size == 0:
        return None
    if nz.size == 1:
        return int(nz[0] + 1)
    # fall back to max abs
    j = int(np.argmax(np.abs(mat[r])) + 1)
    return j


def regressors_of_interest(task: str, design_con: Path) -> List[RegressorTriplet]:
    """
    Deduce the regressor triplets from design.con, following your rules:
      - mid: first 18 contrasts (6*3)
      - sharedreward: first 24 contrasts (8*3), even if there are more.

    Uses /ContrastName{n} strings with '(1)', '(2)', '(3)' suffixes to group triplets.
    """
    task = task.lower()
    if task == "mid":
        n_interest = 18
    elif task == "sharedreward":
        n_interest = 24
    else:
        raise ValueError(f"Unknown task for interest contrasts: {task}")

    names, mat = parse_design_con(design_con)
    if mat is None:
        return []

    # gather basis mapping from first N contrasts
    by_base: Dict[str, Dict[int, int]] = {}
    by_base_cols: Dict[str, Dict[int, int]] = {}
    base_order: List[str] = []

    for cidx in range(1, min(n_interest, mat.shape[0]) + 1):
        cname = names.get(cidx, f"contrast{cidx}")
        m = _BASIS_SUFFIX_RE.match(cname)
        if not m:
            # If naming drifts, skip; still might be ok if you rely on ordering elsewhere.
            continue
        base = m.group(1).strip()
        basis = int(m.group(2))
        if base not in by_base:
            base_order.append(base)
        if basis not in (1, 2, 3):
            continue

        by_base.setdefault(base, {})[basis] = cidx
        col = _row_to_single_col(mat, cidx)
        if col is not None:
            by_base_cols.setdefault(base, {})[basis] = col

    triplets: List[RegressorTriplet] = []
    for base in base_order:
        d = by_base[base]
        cols = by_base_cols.get(base, {})
        if all(b in d for b in (1, 2, 3)) and all(b in cols for b in (1, 2, 3)):
            triplets.append(
                RegressorTriplet(
                    base_name=base,
                    contrast_idx_1=d[1],
                    contrast_idx_2=d[2],
                    contrast_idx_3=d[3],
                    pe_col_1=cols[1],
                    pe_col_2=cols[2],
                    pe_col_3=cols[3],
                )
            )

    # If something went wrong with parsing names, fall back to simple ordering in blocks of 3.
    if not triplets:
        max_rows = min(n_interest, mat.shape[0])
        if max_rows % 3 != 0:
            max_rows = (max_rows // 3) * 3
        for k in range(0, max_rows, 3):
            c1, c2, c3 = k + 1, k + 2, k + 3
            col1 = _row_to_single_col(mat, c1)
            col2 = _row_to_single_col(mat, c2)
            col3 = _row_to_single_col(mat, c3)
            if None in (col1, col2, col3):
                continue
            triplets.append(
                RegressorTriplet(
                    base_name=f"EV{k//3 + 1}",
                    contrast_idx_1=c1,
                    contrast_idx_2=c2,
                    contrast_idx_3=c3,
                    pe_col_1=col1,
                    pe_col_2=col2,
                    pe_col_3=col3,
                )
            )
    return triplets


# ---------------------------
# FEAT dir parsing + discovery
# ---------------------------

# Matches your example well; still tolerant if extra tokens appear.
_FEAT_NAME_RE = re.compile(
    r"^L1_.*?_sub-(?P<sub>\d+)_ses-(?P<ses>\d+)_task-(?P<task>[A-Za-z0-9]+).*?_run-(?P<run>\d+)_space-(?P<space>[^_]+)_(?P<acq>multi-echo|single-echo|multi-echo|single|multiecho|singleecho)_(?:cnfds|confounds)-(?P<conf>[^_]+)_(?P<sm>[^_]+)_FLOBS$"
)

def parse_feat_tokens(feat_dir: Path) -> Dict[str, str]:
    """
    Parse sub/ses/task/run/space/acq/conf/sm from the FEAT directory name.

    If parsing fails, returns 'NA' for fields (but still processes).
    """
    base = feat_dir.name
    if base.endswith(".feat"):
        base = base[:-5]
    m = _FEAT_NAME_RE.match(base)
    if not m:
        # looser fallback: try to grab common tags
        def grab(tag: str) -> str:
            mm = re.search(rf"_{tag}-([^_]+)", base)
            return mm.group(1) if mm else "NA"
        acq = "NA"
        if "_multi-echo_" in base:
            acq = "multiecho"
        elif "_single-echo_" in base or "_single_" in base:
            acq = "single"
        return {
            "sub": grab("sub"),
            "ses": grab("ses"),
            "task": grab("task").lower(),
            "run": grab("run"),
            "space": grab("space"),
            "acq": acq,
            "conf": grab("cnfds") if grab("cnfds") != "NA" else grab("confounds"),
            "sm": grab("sm") if grab("sm") != "NA" else ("unsmoothed" if "unsmoothed" in base else "NA"),
        }

    gd = m.groupdict()
    acq = gd["acq"]
    acq = acq.replace("-", "")
    if acq == "multiecho":
        acq = "multiecho"
    elif acq in ("singleecho", "single"):
        acq = "single"
    return {
        "sub": gd["sub"],
        "ses": gd["ses"],
        "task": gd["task"].lower(),
        "run": gd["run"],
        "space": gd["space"],
        "acq": acq,
        "conf": gd["conf"],
        "sm": gd["sm"],
    }


def find_l1_flobs_feats(deriv_fsl: Path) -> List[Path]:
    feats = sorted(deriv_fsl.glob("sub-*/ses-*/L1_*_FLOBS.feat"))
    return feats


# ---------------------------
# core extraction
# ---------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract L1 FLOBS EV metrics (PE1 zstat + scaled PE1 + scaled signed RMS).")
    p.add_argument("--root", type=Path, default=None, help="Project root (defaults to parent of this script's directory).")
    p.add_argument("--deriv-fsl", type=Path, default=None, help="Override derivatives/fsl directory.")
    p.add_argument("--maskdir", type=Path, default=None, help="Override masks directory.")
    p.add_argument("--outdir", type=Path, default=None, help="Override derivatives/extractions directory.")
    p.add_argument("--outfile", type=Path, default=None, help="Override output TSV filename.")
    p.add_argument("--mni-only", action="store_true", help="Skip FEAT dirs whose space tag does not start with 'mni'.")
    p.add_argument("--flush-every", type=int, default=50, help="Flush TSV every N rows (default 50).")
    p.add_argument("--progress-every", type=int, default=20, help="Print progress every N FEAT dirs (default 20).")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    scriptdir = Path(__file__).resolve().parent
    root = args.root if args.root is not None else scriptdir.parent
    deriv_fsl = args.deriv_fsl if args.deriv_fsl is not None else (root / "derivatives" / "fsl")
    maskdir = args.maskdir if args.maskdir is not None else (root / "masks")
    outdir = args.outdir if args.outdir is not None else (root / "derivatives" / "extractions")
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = args.outfile if args.outfile is not None else (outdir / "extractions_L1_FLOBS_EVmetrics.tsv")

    # Masks/maps (match your bash scripts)
    nacc_mask_path = maskdir / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
    cort_mask_path = maskdir / "BRS_Cortical_3pt1.nii.gz"
    brs_map_path = maskdir / "space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

    for pth in (nacc_mask_path, cort_mask_path, brs_map_path):
        if not pth.exists():
            print(f"ERROR: missing required mask/map: {pth}", file=sys.stderr)
            return 2

    nacc_mask = load_nii_float(nacc_mask_path)
    cort_mask = load_nii_float(cort_mask_path)
    brs_map = load_nii_float(brs_map_path)

    if nacc_mask is None or cort_mask is None or brs_map is None:
        print("ERROR: failed to load one or more masks/maps.", file=sys.stderr)
        return 2

    # interpret masks as boolean (nonzero)
    nacc_bool = nacc_mask > 0
    cort_bool = cort_mask > 0

    feats = find_l1_flobs_feats(deriv_fsl)
    if not feats:
        print(f"[{now_ts()}] WARN: no L1 FLOBS FEAT dirs found under: {deriv_fsl}", file=sys.stderr)
        print(f"[{now_ts()}] Wrote header only to: {outfile}", file=sys.stderr)

    header = [
        "sub", "ses", "run", "task", "space", "acq", "confounds", "sm",
        "regressor",
        "cidx_basis1", "cidx_basis2", "cidx_basis3",
        "pecol_basis1", "pecol_basis2", "pecol_basis3",
        "NAcc_zstat_mean", "NAcc_beta_pe1_scaled_mean", "NAcc_signed_rms_scaled_mean",
        "BRS_Cort_zstat_mean", "BRS_Cort_beta_pe1_scaled_mean", "BRS_Cort_signed_rms_scaled_mean",
        "BRS_corr_zstat", "BRS_corr_pe1_scaled", "BRS_corr_signed_rms_scaled",
        "feat_dir",
    ]

    n_rows_written = 0

    with outfile.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        f.flush()

        t0 = time.time()
        for i, feat_dir in enumerate(feats, start=1):
            tokens = parse_feat_tokens(feat_dir)

            if args.mni_only:
                sp = tokens.get("space", "NA").lower()
                if not sp.startswith("mni"):
                    continue

            task = tokens.get("task", "NA").lower()
            if task not in ("mid", "sharedreward"):
                # keep it strict: masks are MNI, and our interest-contrast rules are task-specific
                continue

            design_con = feat_dir / "design.con"
            design_mat = feat_dir / "design.mat"
            stats_dir = feat_dir / "stats"
            feat_mask_path = feat_dir / "mask.nii.gz"

            # FEAT brainmask (for BRS correlations)
            feat_mask = load_nii_float(feat_mask_path)
            feat_mask_bool = (feat_mask > 0) if feat_mask is not None else None

            # scaling factors
            ppheights = parse_design_mat_ppheights(design_mat)

            # regressor triplets from design.con
            triplets = regressors_of_interest(task, design_con) if design_con.exists() else []

            if not triplets:
                # still write a minimal row indicating failure to parse
                w.writerow([
                    tokens.get("sub", "NA"), tokens.get("ses", "NA"), tokens.get("run", "NA"),
                    task, tokens.get("space", "NA"), tokens.get("acq", "NA"),
                    tokens.get("conf", "NA"), tokens.get("sm", "NA"),
                    NA, NA, NA, NA, NA, NA, NA,
                    NA, NA, NA, NA, NA, NA, NA, NA, NA,
                    str(feat_dir),
                ])
                n_rows_written += 1
                if args.flush_every > 0 and (n_rows_written % args.flush_every == 0):
                    f.flush()
                continue

            for tr in triplets:
                # defaults
                nacc_z = nacc_pe1 = nacc_rms = None
                cort_z = cort_pe1 = cort_rms = None
                brs_corr_z = brs_corr_pe1 = brs_corr_rms = None

                # zstat for basis1 contrast (PE1-like)
                zstat_path = stats_dir / f"zstat{tr.contrast_idx_1}.nii.gz"
                zstat = load_nii_float(zstat_path) if zstat_path.exists() else None

                # scaled PE volumes (basis 1..3)
                pe1 = pe2 = pe3 = None
                pe1_scaled = pe2_scaled = pe3_scaled = None

                def load_scaled_pe(col_1based: int) -> Optional[np.ndarray]:
                    pe_path = stats_dir / f"pe{col_1based}.nii.gz"
                    if not pe_path.exists():
                        return None
                    vol = load_nii_float(pe_path)
                    if vol is None:
                        return None
                    if not ppheights or col_1based > len(ppheights):
                        return None
                    scale = ppheights[col_1based - 1]
                    if scale == 0:
                        return None
                    return vol / float(scale)

                pe1_scaled = load_scaled_pe(tr.pe_col_1)
                pe2_scaled = load_scaled_pe(tr.pe_col_2)
                pe3_scaled = load_scaled_pe(tr.pe_col_3)

                # signed RMS (scaled)
                signed_rms = None
                if pe1_scaled is not None and pe2_scaled is not None and pe3_scaled is not None:
                    # magnitude
                    rms = np.sqrt(pe1_scaled * pe1_scaled + pe2_scaled * pe2_scaled + pe3_scaled * pe3_scaled)
                    # sign from pe1_scaled
                    sgn = np.sign(pe1_scaled)
                    signed_rms = rms * sgn

                # ROI means
                nacc_z = roi_mean(zstat, nacc_bool)
                cort_z = roi_mean(zstat, cort_bool)

                nacc_pe1 = roi_mean(pe1_scaled, nacc_bool)
                cort_pe1 = roi_mean(pe1_scaled, cort_bool)

                nacc_rms = roi_mean(signed_rms, nacc_bool)
                cort_rms = roi_mean(signed_rms, cort_bool)

                # BRS correlations (whole-brain; mask to FEAT brain mask if possible)
                if feat_mask_bool is not None:
                    brs_corr_z = pearson_corr(zstat, brs_map, feat_mask_bool)
                    brs_corr_pe1 = pearson_corr(pe1_scaled, brs_map, feat_mask_bool)
                    brs_corr_rms = pearson_corr(signed_rms, brs_map, feat_mask_bool)

                w.writerow([
                    tokens.get("sub", "NA"), tokens.get("ses", "NA"), tokens.get("run", "NA"),
                    task, tokens.get("space", "NA"), tokens.get("acq", "NA"),
                    tokens.get("conf", "NA"), tokens.get("sm", "NA"),
                    tr.base_name,
                    tr.contrast_idx_1, tr.contrast_idx_2, tr.contrast_idx_3,
                    tr.pe_col_1, tr.pe_col_2, tr.pe_col_3,
                    safe_float(nacc_z), safe_float(nacc_pe1), safe_float(nacc_rms),
                    safe_float(cort_z), safe_float(cort_pe1), safe_float(cort_rms),
                    safe_float(brs_corr_z), safe_float(brs_corr_pe1), safe_float(brs_corr_rms),
                    str(feat_dir),
                ])
                n_rows_written += 1

                if args.flush_every > 0 and (n_rows_written % args.flush_every == 0):
                    f.flush()

            # progress
            if args.progress_every > 0 and (i % args.progress_every == 0):
                elapsed = time.time() - t0
                print(f"[{now_ts()}] processed {i}/{len(feats)} FEAT dirs | wrote {n_rows_written} rows | elapsed {elapsed/60:.1f} min",
                      file=sys.stderr)

        f.flush()

    print(f"[{now_ts()}] Done. Wrote: {outfile}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
