#!/usr/bin/env python3
"""Extract trial-wise metrics from LSS-FLOBS FEAT outputs.

This script is meant to replace the older bash LSS extractor for the new
per-trial FEAT directory structure.

Per trial, it writes:
  - ROI means for zstat{N}
  - ROI means for FLOBS basis-1 beta (PE1), scaled by its /PPheights
  - ROI means for signed RMS across FLOBS basis betas (PE1-PE3), after scaling
    each PE by its own /PPheights
  - Whole-brain correlations against the Brain Reward Signature (BRS) map for:
      zstat{N}, scaled PE1, and scaled signed RMS

Missing trials/outputs are written as 'n/a'.

Expected directory layout (example):
  derivatives/fsl/sub-101/LSS-FLOBS/ses-01/mid/
    L1_task-mid_model-LSS-type-act_sub-101_ses-01_run-2_sm-0_trial-56_
      acq-single_space-MNI152NLin6Asym_confounds-tedana_FLOBS.feat/

Within each .feat:
  - stats/zstat{N}.nii.gz
  - stats/pe{idx}.nii.gz
  - design.mat (for /PPheights)

Masks/maps (defaults match your previous bash script):
  - masks/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz
  - masks/BRS_Cortical_3pt1.nii.gz
  - masks/space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz

Usage:
  python extractData_LSS_FLOBS.py

Helpful overrides:
  python extractData_LSS_FLOBS.py --sm 5
  python extractData_LSS_FLOBS.py --pe-indices 7 8 9
  python extractData_LSS_FLOBS.py --zstat-index 2
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

try:
    import nibabel as nib
except ImportError as e:
    raise SystemExit(
        "ERROR: nibabel is required. In conda, try `conda install -c conda-forge nibabel`."
    ) from e


NA = "n/a"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Extract zstat, scaled PE1 beta, scaled signed RMS (PE1-PE3), and BRS correlations "
            "from LSS-FLOBS trial FEAT outputs."
        )
    )

    p.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Project root (defaults to parent of this script's directory).",
    )
    p.add_argument(
        "--deriv-fsl",
        type=Path,
        default=None,
        help="Override derivatives/fsl directory.",
    )
    p.add_argument(
        "--maskdir",
        type=Path,
        default=None,
        help="Override masks directory.",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Override derivatives/extractions directory.",
    )

    p.add_argument("--sm", type=str, default="0", help="Smoothing tag used in FEAT dir names (default: 0).")
    p.add_argument(
        "--space",
        type=str,
        default="MNI152NLin6Asym",
        help="Space tag used in FEAT dir names (default: MNI152NLin6Asym).",
    )

    p.add_argument(
        "--tasks",
        nargs="+",
        default=["mid", "sharedreward"],
        help="Tasks to extract (default: mid sharedreward).",
    )
    p.add_argument(
        "--runs",
        nargs="+",
        type=int,
        default=[1, 2],
        help="Runs to extract (default: 1 2).",
    )
    p.add_argument(
        "--acqs",
        nargs="+",
        default=["multiecho", "single"],
        help="Acquisition tags to extract (default: multiecho single).",
    )
    p.add_argument(
        "--confs",
        nargs="+",
        default=["base", "tedana"],
        help="Confounds tags to extract (default: base tedana).",
    )

    p.add_argument(
        "--pe-indices",
        nargs=3,
        type=int,
        default=[1, 2, 3],
        metavar=("PE1", "PE2", "PE3"),
        help="PE indices for the FLOBS triplet (default: 1 2 3).",
    )

    p.add_argument(
        "--zstat-index",
        type=int,
        default=1,
        help="Which zstat image to read (stats/zstat{N}.nii.gz). Default: 1.",
    )

    p.add_argument(
        "--outfile",
        type=str,
        default="extractions_LSS-FLOBS_metrics.tsv",
        help="Output TSV filename (written inside outdir).",
    )

    return p.parse_args()


def trials_for(task: str) -> int:
    if task == "mid":
        return 56
    if task == "sharedreward":
        return 54
    return 0


def load_bool_mask(mask_path: Path) -> np.ndarray:
    img = nib.load(str(mask_path))
    data = img.get_fdata(dtype=np.float32)
    mask = data > 0.5
    if int(mask.sum()) == 0:
        raise SystemExit(f"ERROR: Mask has zero voxels after thresholding: {mask_path}")
    return mask


def load_float_img(img_path: Path) -> Optional[np.ndarray]:
    if not img_path.is_file():
        return None
    return nib.load(str(img_path)).get_fdata(dtype=np.float32)


def masked_mean(arr: np.ndarray, mask: np.ndarray) -> Optional[float]:
    vals = arr[mask]
    if vals.size == 0:
        return None
    # (No special NaN handling expected here, but be safe.)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return None
    return float(np.mean(vals))


def pearson_corr(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> Optional[float]:
    """Pearson correlation of two volumes under a boolean mask."""
    va = a[mask]
    vb = b[mask]
    if va.size < 2 or vb.size < 2:
        return None

    good = np.isfinite(va) & np.isfinite(vb)
    va = va[good]
    vb = vb[good]
    if va.size < 2:
        return None

    va = va.astype(np.float64, copy=False)
    vb = vb.astype(np.float64, copy=False)

    va = va - va.mean()
    vb = vb - vb.mean()

    denom = np.sqrt(np.sum(va * va) * np.sum(vb * vb))
    if denom == 0:
        return None

    return float(np.sum(va * vb) / denom)


def safe_fmt(x: Optional[float]) -> str:
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return NA
    return f"{x:.6g}"


def parse_design_mat_ppheights(design_mat: Path) -> Optional[List[float]]:
    """Parse /PPheights from a FEAT design.mat.

    Returns a list of floats in *column order* (length NumWaves) or None.
    """
    if not design_mat.is_file():
        return None

    num_waves: Optional[int] = None
    pp: List[float] = []
    in_pp = False

    with design_mat.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            m = re.match(r"^/NumWaves\s+(\d+)", s)
            if m:
                num_waves = int(m.group(1))
                continue

            if s.startswith("/PPheights"):
                in_pp = True
                parts = s.split()
                for tok in parts[1:]:
                    try:
                        pp.append(float(tok))
                    except ValueError:
                        pass
                continue

            if in_pp:
                if s.startswith("/"):
                    in_pp = False
                    break
                for tok in s.split():
                    try:
                        pp.append(float(tok))
                    except ValueError:
                        pass

    if num_waves is None:
        return pp if pp else None

    if len(pp) < num_waves:
        return None

    return pp[:num_waves]


@dataclass(frozen=True)
class SessionKey:
    sub: str
    ses: str


def discover_sessions(deriv_fsl: Path) -> List[SessionKey]:
    """Find (sub, ses) pairs that have any LSS-FLOBS outputs."""
    keys: set[Tuple[str, str]] = set()
    for ses_dir in deriv_fsl.glob("sub-*/LSS-FLOBS/ses-*"):
        try:
            sub = ses_dir.parts[-3].split("sub-")[1]
            ses = ses_dir.name.split("ses-")[1]
        except Exception:
            continue
        keys.add((sub, ses))
    return [SessionKey(sub=k[0], ses=k[1]) for k in sorted(keys)]


def feat_dir_for(
    deriv_fsl: Path,
    sub: str,
    ses: str,
    task: str,
    run: int,
    trial: int,
    acq: str,
    space: str,
    conf: str,
    sm: str,
) -> Optional[Path]:
    """Return FEAT directory for this combo, trying both unpadded and 2-digit padded trial."""

    trial_candidates = [str(trial), f"{trial:02d}"]

    for tr in trial_candidates:
        feat = (
            deriv_fsl
            / f"sub-{sub}"
            / "LSS-FLOBS"
            / f"ses-{ses}"
            / task
            / (
                f"L1_task-{task}_model-LSS-type-act_"
                f"sub-{sub}_ses-{ses}_run-{run}_sm-{sm}_trial-{tr}_"
                f"acq-{acq}_space-{space}_confounds-{conf}_FLOBS.feat"
            )
        )
        if feat.is_dir():
            return feat

    return None


def main() -> None:
    args = parse_args()

    scriptdir = Path(__file__).resolve().parent
    root = args.root if args.root is not None else scriptdir.parent

    deriv_fsl = args.deriv_fsl if args.deriv_fsl is not None else root / "derivatives" / "fsl"
    maskdir = args.maskdir if args.maskdir is not None else root / "masks"
    outdir = args.outdir if args.outdir is not None else root / "derivatives" / "extractions"
    outdir.mkdir(parents=True, exist_ok=True)

    # Paths (matching the bash script)
    nacc_path = maskdir / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
    cort_path = maskdir / "BRS_Cortical_3pt1.nii.gz"
    brs_map_path = maskdir / "space-MNI152NLin6Asym_desc-BrainRewardSignature_map.nii.gz"

    for p in (nacc_path, cort_path, brs_map_path):
        if not p.is_file():
            raise SystemExit(f"ERROR: Missing mask/map: {p}")

    nacc_mask = load_bool_mask(nacc_path)
    cort_mask = load_bool_mask(cort_path)

    brs_map = load_float_img(brs_map_path)
    if brs_map is None:
        raise SystemExit(f"ERROR: Could not load BRS map: {brs_map_path}")

    sessions = discover_sessions(deriv_fsl)
    if not sessions:
        raise SystemExit(
            f"ERROR: No sessions found under {deriv_fsl}/sub-*/LSS-FLOBS/ses-* . "
            "Check that derivatives are mounted and the path is correct."
        )

    pe1_idx, pe2_idx, pe3_idx = args.pe_indices

    outfile = outdir / args.outfile

    header = [
        "sub",
        "ses",
        "run",
        "task",
        "space",
        "acq",
        "confounds",
        "trial",
        # ROI means
        "NAcc_zstat_mean",
        "NAcc_beta_pe1_scaled_mean",
        "NAcc_signed_rms_scaled_mean",
        "BRS_Cort_zstat_mean",
        "BRS_Cort_beta_pe1_scaled_mean",
        "BRS_Cort_signed_rms_scaled_mean",
        # BRS correlations
        "BRS_corr_zstat",
        "BRS_corr_pe1_scaled",
        "BRS_corr_signed_rms_scaled",
        # PPheights used
        "PPheight_pe1",
        "PPheight_pe2",
        "PPheight_pe3",
        # sanity path
        "feat_dir",
    ]

    with outfile.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)

        total_sessions = len(sessions)
        done_sessions = 0

        for sk in sessions:
            for task in args.tasks:
                ntrials = trials_for(task)
                if ntrials <= 0:
                    continue

                for run in args.runs:
                    for acq in args.acqs:
                        for conf in args.confs:
                            for trial in range(1, ntrials + 1):
                                feat_dir = feat_dir_for(
                                    deriv_fsl=deriv_fsl,
                                    sub=sk.sub,
                                    ses=sk.ses,
                                    task=task,
                                    run=run,
                                    trial=trial,
                                    acq=acq,
                                    space=args.space,
                                    conf=conf,
                                    sm=args.sm,
                                )

                                # Defaults
                                nacc_z = nacc_b1 = nacc_rms = None
                                cort_z = cort_b1 = cort_rms = None
                                brs_corr_z = brs_corr_pe1 = brs_corr_rms = None
                                pp1 = pp2 = pp3 = None

                                if feat_dir is not None:
                                    stats_dir = feat_dir / "stats"

                                    # --- zstat (contrast-level statistic) ---
                                    zstat_path = stats_dir / f"zstat{args.zstat_index}.nii.gz"
                                    zstat = load_float_img(zstat_path)

                                    wbmask: Optional[np.ndarray] = None
                                    if zstat is not None:
                                        if zstat.shape != brs_map.shape:
                                            raise SystemExit(
                                                f"ERROR: Shape mismatch between zstat and BRS map.\n"
                                                f"  zstat: {zstat_path} shape={zstat.shape}\n"
                                                f"  brs:   {brs_map_path} shape={brs_map.shape}\n"
                                                "This usually means the files are not in the same space/grid."
                                            )

                                        # Whole-brain mask like the bash script: abs(z) > 0
                                        wbmask = (np.abs(zstat) > 0) & np.isfinite(zstat)

                                        nacc_z = masked_mean(zstat, nacc_mask)
                                        cort_z = masked_mean(zstat, cort_mask)

                                        brs_corr_z = pearson_corr(zstat, brs_map, wbmask)

                                    # --- design scaling (PPheights) ---
                                    ppheights = parse_design_mat_ppheights(feat_dir / "design.mat")
                                    if ppheights and max(pe1_idx, pe2_idx, pe3_idx) <= len(ppheights):
                                        pp1 = ppheights[pe1_idx - 1]
                                        pp2 = ppheights[pe2_idx - 1]
                                        pp3 = ppheights[pe3_idx - 1]

                                    # --- FLOBS metrics (scaled PE1 and scaled signed RMS) ---
                                    pe1_path = stats_dir / f"pe{pe1_idx}.nii.gz"
                                    pe2_path = stats_dir / f"pe{pe2_idx}.nii.gz"
                                    pe3_path = stats_dir / f"pe{pe3_idx}.nii.gz"

                                    if (
                                        pp1 is not None
                                        and pp2 is not None
                                        and pp3 is not None
                                        and pp1 != 0
                                        and pp2 != 0
                                        and pp3 != 0
                                        and pe1_path.is_file()
                                        and pe2_path.is_file()
                                        and pe3_path.is_file()
                                    ):
                                        pe1 = load_float_img(pe1_path)
                                        pe2 = load_float_img(pe2_path)
                                        pe3 = load_float_img(pe3_path)

                                        if pe1 is None or pe2 is None or pe3 is None:
                                            # Leave as NA
                                            pass
                                        else:
                                            if pe1.shape != brs_map.shape:
                                                raise SystemExit(
                                                    f"ERROR: Shape mismatch between PE and BRS map.\n"
                                                    f"  pe1: {pe1_path} shape={pe1.shape}\n"
                                                    f"  brs: {brs_map_path} shape={brs_map.shape}"
                                                )

                                            pe1_s = pe1 / float(pp1)
                                            pe2_s = pe2 / float(pp2)
                                            pe3_s = pe3 / float(pp3)

                                            signed_rms = np.sign(pe1_s) * np.sqrt(pe1_s * pe1_s + pe2_s * pe2_s + pe3_s * pe3_s)

                                            nacc_b1 = masked_mean(pe1_s, nacc_mask)
                                            cort_b1 = masked_mean(pe1_s, cort_mask)

                                            nacc_rms = masked_mean(signed_rms, nacc_mask)
                                            cort_rms = masked_mean(signed_rms, cort_mask)

                                            # For correlations, prefer the same wbmask used for zstat
                                            # (matches the older bash behavior). If zstat was missing,
                                            # fall back to a PE-based mask.
                                            if wbmask is None:
                                                wbmask = (np.abs(pe1_s) > 0) & np.isfinite(pe1_s)

                                            brs_corr_pe1 = pearson_corr(pe1_s, brs_map, wbmask)
                                            brs_corr_rms = pearson_corr(signed_rms, brs_map, wbmask)

                                w.writerow(
                                    [
                                        sk.sub,
                                        sk.ses,
                                        str(run),
                                        task,
                                        args.space,
                                        acq,
                                        conf,
                                        str(trial),
                                        safe_fmt(nacc_z),
                                        safe_fmt(nacc_b1),
                                        safe_fmt(nacc_rms),
                                        safe_fmt(cort_z),
                                        safe_fmt(cort_b1),
                                        safe_fmt(cort_rms),
                                        safe_fmt(brs_corr_z),
                                        safe_fmt(brs_corr_pe1),
                                        safe_fmt(brs_corr_rms),
                                        safe_fmt(pp1),
                                        safe_fmt(pp2),
                                        safe_fmt(pp3),
                                        str(feat_dir) if feat_dir is not None else NA,
                                    ]
                                )

            done_sessions += 1
            pct = int(round(100 * done_sessions / max(total_sessions, 1)))
            print(f"[{done_sessions}/{total_sessions}] {pct}% sessions processed (last: sub-{sk.sub} ses-{sk.ses})")

    print(f"Done. Wrote: {outfile}")


if __name__ == "__main__":
    main()
