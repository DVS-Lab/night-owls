#!/usr/bin/env python3
"""
VS time-course QC for MID: event-locked ROI timecourses with:
- Cue-offset anchored anticipation (t=0 at cue offset; cue onset at -0.75s)
- Wu-analog stimulus window shading (0–2s after cue offset)
- Clarifying timing markers: cue onset/offset, median target onset, median feedback onset, +6s
- Optional nuisance regression using FEAT confoundevs.txt
- Dual curves on each plot: raw (solid) and confound-regressed (dashed), matched colors

CRITICAL: Keep existing conventions from the original script
  FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-unsmoothed.txt"
  VS_MNI         = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
"""

from __future__ import annotations

import argparse
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

try:
    from nibabel.processing import resample_from_to
except Exception:
    resample_from_to = None


# ------------------------------- Paths / Conventions -------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
MASKS_DIR  = (SCRIPT_DIR / "../masks").resolve()

# Keep your existing conventions from the original script
FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-unsmoothed.txt"
VS_MNI         = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"

# EV layout (kept from prior versions that were already working)
EV_ROOT = (SCRIPT_DIR / "../EVFiles").resolve()

# Outputs
OUT_DIR = (SCRIPT_DIR / "../derivatives/extractions/timecourses-mid-unsmoothed_markers_dualconf").resolve()


# ------------------------------- Analysis settings -------------------------------

# Cue duration (s) used for cue-offset anchoring (t=0 at cue offset)
CUE_DUR_S = 0.75

# Time axis for event-locked curves (seconds)
T_MIN_S = -4.0
T_MAX_S = 16.0
DT_S = 0.10

# Baseline subtraction window (seconds, relative to t=0)
BASELINE_WIN = (-4.0, -2.0)

# Shade Wu "anticipation stimulus epoch" (seconds, relative to cue offset)
WU_ANT_WINDOW = (0.0, 2.0)

# TR sanity check: compare header TR to design.fsf TR (warn if mismatch)
TR_MISMATCH_TOL = 1e-3  # seconds

# Confound regression from FEAT confoundevs.txt
APPLY_CONFOUND_REGRESSION = True


# ------------------------------- Data structures -------------------------------

@dataclass
class RunResult:
    sub: str
    ses: str
    run: str
    echo: str  # "single-echo" or "multi-echo"
    feat_dir: Path

    tr_header: float
    tr_design: Optional[float]
    tr_used: float

    time_axis: np.ndarray

    # per-event timing diagnostics (seconds relative to cue offset)
    isi_s: np.ndarray           # target onset ≈ (ant_dur - cue_dur)
    fb_delay_s: np.ndarray      # feedback onset relative to cue offset (next feedback after cue offset)

    # Event-locked curves: dict condition -> (mean, sem, n)
    ant_psc_raw: Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    ant_psc_conf: Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_psc_raw: Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_psc_conf: Dict[str, Tuple[np.ndarray, np.ndarray, int]]

    # bookkeeping
    confounds_used: bool
    n_confounds: int


# ------------------------------- Helpers -------------------------------

def _read_text(path: Path) -> str:
    return path.read_text(errors="ignore")


def normalize_feat_path(p: str) -> Path:
    """Accept lines that may point inside a .feat directory; return the .feat dir itself."""
    s = p.strip().strip('"').strip("'")
    if not s:
        return Path()
    pp = Path(s)
    if not pp.is_absolute():
        # resolve relative to the script dir (matches typical "run from code/" workflows)
        pp = (SCRIPT_DIR / pp).resolve()

    # If the path includes ".feat/..." trim to the .feat directory
    parts = pp.parts
    if any(x.endswith(".feat") for x in parts):
        # find first part ending with .feat
        idx = next(i for i, x in enumerate(parts) if x.endswith(".feat"))
        pp = Path(*parts[: idx + 1]).resolve()

    return pp


def load_feat_list(path: Path) -> List[Path]:
    if not path.exists():
        raise FileNotFoundError(f"Missing FEAT list: {path}")
    feats: List[Path] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fp = normalize_feat_path(line)
        if fp and fp.exists() and fp.is_dir() and str(fp).endswith(".feat"):
            feats.append(fp)
    return feats


def parse_ids_from_feat(feat: Path) -> Tuple[str, str, str, str]:
    """
    Try to parse sub/ses/run and echo type from FEAT path.
    Echo mapping:
      - contains "tedana" or "multi" -> multi-echo
      - contains "fmriprep" or "single" -> single-echo
      - fallback: "single-echo"
    """
    p = str(feat)
    msub = re.search(r"sub-([A-Za-z0-9]+)", p)
    mses = re.search(r"ses-([A-Za-z0-9]+)", p)
    mrun = re.search(r"run-([0-9]+)", p)
    if not (msub and mses and mrun):
        raise ValueError(f"Could not parse sub/ses/run from: {feat}")
    sub, ses, run = msub.group(1), mses.group(1), mrun.group(1)

    low = p.lower()
    if "tedana" in low or "multi" in low:
        echo = "multi-echo"
    elif "fmriprep" in low or "single" in low:
        echo = "single-echo"
    else:
        echo = "single-echo"
    return sub, ses, run, echo


def load_tr_from_design_fsf(feat: Path) -> Optional[float]:
    fsf = feat / "design.fsf"
    if not fsf.exists():
        return None
    txt = _read_text(fsf)
    m = re.search(r"set\s+fmri\(tr\)\s+([0-9]*\.?[0-9]+)", txt)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def load_ev_3col(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not path.exists():
        return np.array([], float), np.array([], float), np.array([], float)
    try:
        arr = np.loadtxt(path, ndmin=2)
    except Exception:
        return np.array([], float), np.array([], float), np.array([], float)
    if arr.size == 0:
        return np.array([], float), np.array([], float), np.array([], float)
    if arr.shape[1] < 3:
        raise ValueError(f"EV file must have 3 columns: {path}")
    return arr[:, 0].astype(float), arr[:, 1].astype(float), arr[:, 2].astype(float)


def get_ev_dir(sub: str, ses: str, run: str) -> Path:
    return EV_ROOT / f"sub-{sub}" / f"ses-{ses}" / f"run-{run}"


def _mask_to_bool(mask_img: nib.Nifti1Image, data_img: nib.Nifti1Image) -> np.ndarray:
    """
    Return a boolean mask aligned to data_img space.
    - If shapes/affines match, use directly.
    - Else, attempt nearest-neighbor resampling (if nibabel.processing is available).
    """
    mdata = mask_img.get_fdata()
    if mdata.ndim == 4:
        mdata = mdata[..., 0]
    mbool = mdata > 0.5

    if mbool.shape == data_img.shape[:3] and np.allclose(mask_img.affine, data_img.affine, atol=1e-3):
        return mbool

    if resample_from_to is None:
        raise RuntimeError(
            "Mask/data grid mismatch and nibabel.processing.resample_from_to is unavailable.\n"
            f"Mask shape={mbool.shape}, data shape={data_img.shape[:3]}"
        )

    # Resample mask into data grid using nearest-neighbor
    res = resample_from_to(mask_img, (data_img.shape[:3], data_img.affine), order=0)
    rdat = res.get_fdata()
    if rdat.ndim == 4:
        rdat = rdat[..., 0]
    return (rdat > 0.5)


def extract_roi_timeseries(feat: Path, vs_mask_path: Path) -> Tuple[np.ndarray, float]:
    func_path = feat / "filtered_func_data.nii.gz"
    if not func_path.exists():
        raise FileNotFoundError(f"Missing: {func_path}")
    if not vs_mask_path.exists():
        raise FileNotFoundError(f"Missing VS mask: {vs_mask_path}")

    img = nib.load(str(func_path))
    data = img.get_fdata()
    tr = float(img.header.get_zooms()[3])

    mimg = nib.load(str(vs_mask_path))
    mask = _mask_to_bool(mimg, img)

    ts = data[mask, :].mean(axis=0)
    return ts.astype(float), tr


def regress_out_confounds(ts: np.ndarray, conf_path: Path) -> Tuple[np.ndarray, bool, int]:
    """
    OLS regression: ts ~ [1, confounds], return residualized ts with mean added back.
    - Drops columns with ~0 variance.
    - If confounds missing or length mismatch, returns original ts.
    """
    if not conf_path.exists():
        return ts, False, 0
    try:
        X = np.loadtxt(conf_path, ndmin=2)
    except Exception:
        return ts, False, 0

    if X.size == 0:
        return ts, False, 0
    if X.shape[0] != ts.shape[0]:
        return ts, False, 0

    # Drop (near) constant columns
    keep = np.nanstd(X, axis=0) > 1e-8
    Xk = X[:, keep] if keep.any() else np.zeros((X.shape[0], 0), float)
    n_conf = int(Xk.shape[1])

    # Build design with intercept
    Xd = np.column_stack([np.ones((X.shape[0], 1)), Xk])  # (T, 1+n_conf)
    y = ts.reshape(-1, 1)

    # Solve least squares
    try:
        beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
        yhat = Xd @ beta
        resid = (y - yhat).ravel()
    except Exception:
        return ts, False, 0

    # Add back mean to keep PSC well-behaved
    resid = resid + np.mean(ts)
    return resid.astype(float), True, n_conf


def to_psc(ts: np.ndarray) -> np.ndarray:
    m = float(np.mean(ts))
    if not np.isfinite(m) or abs(m) < 1e-12:
        return ts * np.nan
    return (ts - m) / m * 100.0


def event_locked_mean_sem(
    ts: np.ndarray,
    tr: float,
    onsets_s: np.ndarray,
    time_axis: np.ndarray,
    baseline_win: Tuple[float, float]
) -> Tuple[np.ndarray, np.ndarray, int]:
    if onsets_s.size == 0:
        return np.full_like(time_axis, np.nan), np.full_like(time_axis, np.nan), 0

    n_tp = ts.shape[0]
    tpoints = np.arange(n_tp) * tr

    b0, b1 = baseline_win
    bmask = (time_axis >= b0) & (time_axis <= b1)
    if not np.any(bmask):
        raise ValueError("Baseline window does not overlap the time axis.")

    curves = []
    for o in onsets_s:
        sample_t = o + time_axis
        samp = np.interp(sample_t, tpoints, ts, left=np.nan, right=np.nan)
        b = np.nanmean(samp[bmask])
        curves.append(samp - b)

    mat = np.stack(curves, axis=0)
    mean = np.nanmean(mat, axis=0)
    if mat.shape[0] > 1:
        sem = np.nanstd(mat, axis=0, ddof=1) / math.sqrt(mat.shape[0])
    else:
        sem = np.zeros_like(mean)
    return mean, sem, int(mat.shape[0])


def weighted_mean_and_sem(curves: List[Tuple[np.ndarray, int]]) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Weighted aggregation across runs.
      curves: list of (mean_curve, n_events)
    SEM here is across runs (not trials), weighted by n_events.
    """
    curves = [(c, n) for c, n in curves if n > 0 and c is not None and np.any(np.isfinite(c))]
    if not curves:
        # Return NaNs to make "no data" obvious.
        return None, None, 0  # type: ignore

    mats = np.stack([c for c, _ in curves], axis=0)
    w = np.array([n for _, n in curves], dtype=float)
    w = w / np.sum(w)
    total_n = int(np.sum([n for _, n in curves]))

    wmean = np.sum(mats * w[:, None], axis=0)
    if mats.shape[0] == 1:
        sem = np.zeros_like(wmean)
    else:
        wvar = np.sum(w[:, None] * (mats - wmean[None, :]) ** 2, axis=0)
        sem = np.sqrt(wvar) / math.sqrt(mats.shape[0])

    return wmean, sem, total_n


def _nanmedian(a: np.ndarray) -> float:
    if a is None or a.size == 0:
        return float("nan")
    return float(np.nanmedian(a))


def plot_dual_two_conditions(
    time_axis: np.ndarray,
    rawA: Tuple[np.ndarray, np.ndarray, int], confA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
    rawB: Tuple[np.ndarray, np.ndarray, int], confB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
    title: str, ylabel: str, out_png: Path,
    vlines: List[Tuple[float, str]] | None = None,
    spans: List[Tuple[float, float, str]] | None = None,
    show_sem_for: str = "raw"
) -> None:
    """
    Two conditions, each with raw (solid) and confound-regressed (dashed).
    Colors are matched between raw and confound for each condition.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Shaded spans
    if spans:
        for x0, x1, lab in spans:
            ax.axvspan(x0, x1, alpha=0.10, label=lab)

    # Condition A
    mA, sA, nA = rawA
    lineA = ax.plot(time_axis, mA, label=f"{labelA} (raw; n={nA})")[0]
    cA = lineA.get_color()
    if show_sem_for == "raw":
        ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.20, color=cA)
    mAc, sAc, nAc = confA
    if nAc > 0 and np.any(np.isfinite(mAc)):
        ax.plot(time_axis, mAc, ls="--", color=cA, label=f"{labelA} (conf; n={nAc})")

    # Condition B
    mB, sB, nB = rawB
    lineB = ax.plot(time_axis, mB, label=f"{labelB} (raw; n={nB})")[0]
    cB = lineB.get_color()
    if show_sem_for == "raw":
        ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.20, color=cB)
    mBc, sBc, nBc = confB
    if nBc > 0 and np.any(np.isfinite(mBc)):
        ax.plot(time_axis, mBc, ls="--", color=cB, label=f"{labelB} (conf; n={nBc})")

    # Vertical markers
    if vlines:
        for x, _lab in vlines:
            if np.isfinite(x):
                ax.axvline(x, ls=":", lw=1)

        # Labels at top
        y_top = ax.get_ylim()[1]
        for x, lab in vlines:
            if np.isfinite(x) and lab:
                ax.text(x, y_top, lab, rotation=90, va="top", ha="right", fontsize=8)

    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ------------------------------- Core processing -------------------------------

def process_feat(feat: Path, time_axis: np.ndarray) -> RunResult:
    sub, ses, run, echo = parse_ids_from_feat(feat)

    # ROI timeseries
    ts_raw, tr_header = extract_roi_timeseries(feat, VS_MNI)

    # TR sanity check vs design.fsf
    tr_design = load_tr_from_design_fsf(feat)
    tr_used = float(tr_header)
    if tr_design is not None and abs(tr_design - tr_header) > TR_MISMATCH_TOL:
        print(f"[WARN] TR mismatch: {feat}\n  header={tr_header:.6f}, design={tr_design:.6f}")

    # Confounds
    confound_path = feat / "confoundevs.txt"
    ts_conf = ts_raw
    conf_used = False
    n_conf = 0
    if APPLY_CONFOUND_REGRESSION:
        ts_conf, conf_used, n_conf = regress_out_confounds(ts_raw, confound_path)

    # Convert to PSC
    psc_raw = to_psc(ts_raw)
    psc_conf = to_psc(ts_conf)

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)

    # Anticipation events: cue onset + duration. Anchor at cue offset (add 0.75s)
    antR_on, antR_dur, _ = load_ev_3col(ev_dir / "_anticipation_reward.txt")
    antN_on, antN_dur, _ = load_ev_3col(ev_dir / "_anticipation_neutral.txt")
    antR = antR_on + CUE_DUR_S
    antN = antN_on + CUE_DUR_S

    # Target onset relative to cue offset ≈ (ant_duration - cue_dur)
    isi_s = np.concatenate([antR_dur, antN_dur]) - CUE_DUR_S
    isi_s = isi_s[np.isfinite(isi_s)]

    # Feedback onsets (pooled by valence)
    fb_PR_on = load_ev_3col(ev_dir / "_feedback_positive_reward.txt")[0]
    fb_NR_on = load_ev_3col(ev_dir / "_feedback_negative_reward.txt")[0]
    fb_PN_on = load_ev_3col(ev_dir / "_feedback_positive_neutral.txt")[0]
    fb_NN_on = load_ev_3col(ev_dir / "_feedback_negative_neutral.txt")[0]

    fb_pos = np.sort(np.concatenate([fb_PR_on, fb_PN_on]))
    fb_neg = np.sort(np.concatenate([fb_NR_on, fb_NN_on]))

    # Feedback delay diagnostic: next feedback after cue offset
    cueoffs = np.sort(np.concatenate([antR, antN]))
    fb_all = np.sort(np.concatenate([fb_pos, fb_neg]))
    if cueoffs.size == 0 or fb_all.size == 0:
        fb_delay_s = np.array([], float)
    else:
        idx_next = np.searchsorted(fb_all, cueoffs, side="right")
        valid = idx_next < fb_all.size
        fb_delay_s = fb_all[idx_next[valid]] - cueoffs[valid]

    # Event-locked curves
    ant_psc_raw = {
        "Reward": event_locked_mean_sem(psc_raw, tr_used, antR, time_axis, BASELINE_WIN),
        "Neutral": event_locked_mean_sem(psc_raw, tr_used, antN, time_axis, BASELINE_WIN),
    }
    ant_psc_conf = {
        "Reward": event_locked_mean_sem(psc_conf, tr_used, antR, time_axis, BASELINE_WIN),
        "Neutral": event_locked_mean_sem(psc_conf, tr_used, antN, time_axis, BASELINE_WIN),
    }
    fb_psc_raw = {
        "Positive": event_locked_mean_sem(psc_raw, tr_used, fb_pos, time_axis, BASELINE_WIN),
        "Negative": event_locked_mean_sem(psc_raw, tr_used, fb_neg, time_axis, BASELINE_WIN),
    }
    fb_psc_conf = {
        "Positive": event_locked_mean_sem(psc_conf, tr_used, fb_pos, time_axis, BASELINE_WIN),
        "Negative": event_locked_mean_sem(psc_conf, tr_used, fb_neg, time_axis, BASELINE_WIN),
    }

    return RunResult(
        sub=sub, ses=ses, run=run, echo=echo, feat_dir=feat,
        tr_header=float(tr_header), tr_design=tr_design, tr_used=tr_used,
        time_axis=time_axis,
        isi_s=isi_s, fb_delay_s=fb_delay_s,
        ant_psc_raw=ant_psc_raw, ant_psc_conf=ant_psc_conf,
        fb_psc_raw=fb_psc_raw, fb_psc_conf=fb_psc_conf,
        confounds_used=conf_used, n_confounds=n_conf
    )


def save_run_plots(r: RunResult) -> None:
    out = OUT_DIR / "runs" / f"sub-{r.sub}" / f"ses-{r.ses}" / f"run-{r.run}" / r.echo

    med_isi = _nanmedian(r.isi_s)
    med_fb = _nanmedian(r.fb_delay_s)
    ant_vlines = [
        (-CUE_DUR_S, "Cue onset"),
        (0.0, "Cue offset"),
        (med_isi, "Target (median)"),
        (med_fb, "Feedback (median)"),
        (6.0, "+6s"),
    ]
    ant_spans = [(WU_ANT_WINDOW[0], WU_ANT_WINDOW[1], "Wu ANT window (0–2s)")]

    fb_vlines = [
        (0.0, "Feedback onset"),
        (6.0, "+6s"),
    ]

    plot_dual_two_conditions(
        r.time_axis,
        r.ant_psc_raw["Reward"], r.ant_psc_conf["Reward"], "Reward",
        r.ant_psc_raw["Neutral"], r.ant_psc_conf["Neutral"], "Neutral",
        title=f"VS — sub {r.sub} [{r.echo}] (ANT, PSC) | ses {r.ses} run {r.run}",
        ylabel="% signal change (PSC)",
        out_png=out / "anticipation_psc.png",
        vlines=ant_vlines,
        spans=ant_spans,
    )

    plot_dual_two_conditions(
        r.time_axis,
        r.fb_psc_raw["Positive"], r.fb_psc_conf["Positive"], "Feedback +",
        r.fb_psc_raw["Negative"], r.fb_psc_conf["Negative"], "Feedback −",
        title=f"VS — sub {r.sub} [{r.echo}] (FB, PSC) | ses {r.ses} run {r.run}",
        ylabel="% signal change (PSC)",
        out_png=out / "feedback_psc.png",
        vlines=fb_vlines,
        spans=None,
    )


def aggregate_subject(results: List[RunResult], sub: str, echo: str, time_axis: np.ndarray) -> None:
    runs = [r for r in results if r.sub == sub and r.echo == echo]
    if not runs:
        return

    # pooled timing diagnostics across runs
    isi_all = np.concatenate([r.isi_s for r in runs if r.isi_s is not None and r.isi_s.size > 0]) if runs else np.array([], float)
    fb_all = np.concatenate([r.fb_delay_s for r in runs if r.fb_delay_s is not None and r.fb_delay_s.size > 0]) if runs else np.array([], float)

    med_isi = _nanmedian(isi_all)
    med_fb = _nanmedian(fb_all)

    ant_vlines = [
        (-CUE_DUR_S, "Cue onset"),
        (0.0, "Cue offset"),
        (med_isi, "Target (median)"),
        (med_fb, "Feedback (median)"),
        (6.0, "+6s"),
    ]
    ant_spans = [(WU_ANT_WINDOW[0], WU_ANT_WINDOW[1], "Wu ANT window (0–2s)")]

    fb_vlines = [
        (0.0, "Feedback onset"),
        (6.0, "+6s"),
    ]

    # Weighted aggregation (raw + conf) for ANT
    ant_raw_R = weighted_mean_and_sem([(r.ant_psc_raw["Reward"][0], r.ant_psc_raw["Reward"][2]) for r in runs])
    ant_raw_N = weighted_mean_and_sem([(r.ant_psc_raw["Neutral"][0], r.ant_psc_raw["Neutral"][2]) for r in runs])
    ant_con_R = weighted_mean_and_sem([(r.ant_psc_conf["Reward"][0], r.ant_psc_conf["Reward"][2]) for r in runs])
    ant_con_N = weighted_mean_and_sem([(r.ant_psc_conf["Neutral"][0], r.ant_psc_conf["Neutral"][2]) for r in runs])

    fb_raw_P = weighted_mean_and_sem([(r.fb_psc_raw["Positive"][0], r.fb_psc_raw["Positive"][2]) for r in runs])
    fb_raw_N = weighted_mean_and_sem([(r.fb_psc_raw["Negative"][0], r.fb_psc_raw["Negative"][2]) for r in runs])
    fb_con_P = weighted_mean_and_sem([(r.fb_psc_conf["Positive"][0], r.fb_psc_conf["Positive"][2]) for r in runs])
    fb_con_N = weighted_mean_and_sem([(r.fb_psc_conf["Negative"][0], r.fb_psc_conf["Negative"][2]) for r in runs])

    out = OUT_DIR / "subjects" / f"sub-{sub}" / echo
    plot_dual_two_conditions(
        time_axis,
        ant_raw_R, ant_con_R, "Reward",
        ant_raw_N, ant_con_N, "Neutral",
        title=f"VS — subject {sub} [{echo}] (ANT, PSC)",
        ylabel="% signal change (PSC)",
        out_png=out / "anticipation_psc.png",
        vlines=ant_vlines,
        spans=ant_spans,
    )
    plot_dual_two_conditions(
        time_axis,
        fb_raw_P, fb_con_P, "Feedback +",
        fb_raw_N, fb_con_N, "Feedback −",
        title=f"VS — subject {sub} [{echo}] (FB, PSC)",
        ylabel="% signal change (PSC)",
        out_png=out / "feedback_psc.png",
        vlines=fb_vlines,
        spans=None,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only-sub", default=None, help="Only process a subject id (e.g., 105)")
    ap.add_argument("--no-run-plots", action="store_true", help="Skip per-run plots")
    args = ap.parse_args()

    feats = load_feat_list(FEAT_LIST_PATH)
    if args.only_sub:
        feats = [f for f in feats if f"sub-{args.only_sub}" in str(f)]
    if not feats:
        raise SystemExit("No FEAT directories found. Check FEAT_LIST_PATH and its contents.")

    time_axis = np.arange(T_MIN_S, T_MAX_S + 1e-9, DT_S)

    results: List[RunResult] = []
    for feat in feats:
        r = process_feat(feat, time_axis)
        results.append(r)
        if not args.no_run_plots:
            save_run_plots(r)

    # Subject-level aggregation
    subjects = sorted({r.sub for r in results})
    echos = sorted({r.echo for r in results})
    for sub in subjects:
        for echo in echos:
            aggregate_subject(results, sub, echo, time_axis)

    # Print TR distribution and confound usage summary
    trs = np.array([r.tr_header for r in results], float)
    if trs.size:
        vals, counts = np.unique(np.round(trs, 6), return_counts=True)
        print("TR(header) distribution:", ", ".join([f"{v:.6f} ({c})" for v, c in zip(vals, counts)]))

    conf_used = sum(1 for r in results if r.confounds_used)
    print(f"Confounds used in {conf_used}/{len(results)} runs (APPLY_CONFOUND_REGRESSION={APPLY_CONFOUND_REGRESSION}).")
    print(f"Output: {OUT_DIR}")


if __name__ == "__main__":
    main()
