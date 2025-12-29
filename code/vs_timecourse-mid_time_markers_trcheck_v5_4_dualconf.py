#!/usr/bin/env python3
"""
MID VS time-course extraction (event-locked selective averaging) with timing markers
and side-by-side comparison of curves WITH vs WITHOUT confound regression.

Core choices for Wu comparability
- ANT is anchored at *cue offset / delay onset* by shifting cue-onset EVs by +0.75 s.
- Wu's "anticipation period" is a 0–2 s stimulus-time window after cue offset.
  This script shades that window (0–2 s) for interpretability, but note: the BOLD
  peak for that window is expected later (roughly 4–8 s post-window onset).

What this script produces
- Subject-level plots (pooled across runs/sessions), separately for multi-echo vs single-echo:
    * ANT (PSC) and FB (PSC), each with:
        - solid lines: raw ROI time series
        - dashed lines: confound-regressed ROI time series (if confounds available)
      Shaded bands are SEM across runs (weighted), shown for the RAW curves.
- A run-level TSV summary including TR sanity checks and confound usage.

Notes on confound regression
- We regress confounds out of the ROI time series using OLS with an intercept.
- We add the original mean back to the residual so PSC remains well-defined.
- This is not a full FEAT-equivalent residualization; it is a lightweight nuisance cleanup.
"""

from __future__ import annotations

import glob
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt


# ------------------------------- Configuration -------------------------------

# Where FEATs live (relative to where you run the script)
FSL_DERIV = Path("derivatives/fsl")

# Where EV files live (relative)
EV_ROOT = Path("EVFiles")

# Output root
OUT_DIR = Path("derivatives/extractions/timecourses-mid-unsmoothed_markers_dualconf")
SUMMARY_TSV = OUT_DIR / "summary_mid_markers_dualconf.tsv"

# Optional: provide a plain-text list of FEAT directories (one per line) to override glob discovery.
FEAT_LIST_PATH = Path("feat_list.txt")  # if exists, will be used

# Within each FEAT directory
FEAT_FUNC = Path("filtered_func_data.nii.gz")
VS_MASK_REL = Path("mask_vstriatum.nii.gz")
CONFOUNDEVS_NAME = "confoundevs.txt"  # expected under <feat>/confoundevs.txt OR <feat>/confoundevs/confoundevs.txt

# Timing
TR_MISMATCH_WARN_EPS = 1e-3

# Cue duration (s): used to shift cue-onset EVs to cue-offset / delay onset.
CUE_DUR_S = 0.75

# Time grid for interpolation (seconds)
T_MIN_S = -4.0
T_MAX_S = 16.0
DT_S = 0.10

# Baseline subtraction window (in event-locked time)
BASELINE_WIN = (-4.0, -2.0)

# Prefer TR from design.fsf if it disagrees with NIfTI header TR?
PREFER_TR_FROM_DESIGN = False

# Compute confound-regressed variant (dashed lines) if confounds exist and match length
COMPUTE_CONFOUND_VARIANT = True


# ------------------------------- Data structures -----------------------------

Curve = Tuple[np.ndarray, np.ndarray, int]  # (mean, sem, n_events)

@dataclass
class RunResult:
    sub: str
    ses: str
    run: str
    echo: str  # e.g., "fmriprep" or "tedana" or "single-echo" / "multi-echo" labels in your tree
    time_axis: np.ndarray

    tr_header: float
    tr_design: Optional[float]
    tr_used: float

    # Diagnostics for markers (seconds)
    isi_s: np.ndarray        # target-onset relative to cue-offset (from ANT EV durations)
    fb_delay_s: np.ndarray   # feedback-onset relative to cue-offset (paired as next feedback after cue-offset)

    # Confound info
    confounds_used: int      # 1 if dashed curves are true nuisance-regressed curves, else 0
    n_confounds: int         # number of retained confound columns

    # Curves (RAW)
    ant_psc_raw: Dict[str, Curve]  # Reward/Neutral
    fb_psc_raw: Dict[str, Curve]   # Positive/Negative

    # Curves (CONF-REG) - present only if confounds_used==1; otherwise empty dicts
    ant_psc_conf: Dict[str, Curve]
    fb_psc_conf: Dict[str, Curve]


# ------------------------------- Helpers -------------------------------------

def nanmedian_or_nan(a: np.ndarray) -> float:
    if a is None or a.size == 0:
        return float("nan")
    return float(np.nanmedian(a))


def _fmt(x: Optional[float]) -> str:
    if x is None:
        return ""
    if not np.isfinite(x):
        return ""
    return f"{x:.6f}"


def to_psc(ts: np.ndarray) -> np.ndarray:
    m = float(np.mean(ts))
    if not np.isfinite(m) or m == 0.0:
        return np.full_like(ts, np.nan, dtype=float)
    return (ts - m) / m * 100.0


def load_ev_3col(ev_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load an FSL 3-col EV file: onset(s), duration(s), weight. Returns empty arrays if missing."""
    if not ev_path.exists():
        return (np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float))
    try:
        arr = np.loadtxt(ev_path, ndmin=2)
    except Exception:
        return (np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float))
    if arr.size == 0 or arr.shape[0] == 0:
        return (np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float))
    if arr.shape[1] < 3:
        raise ValueError(f"EV file does not have 3 columns: {ev_path}")
    return (arr[:, 0].astype(float), arr[:, 1].astype(float), arr[:, 2].astype(float))


def load_ev_onsets(ev_path: Path) -> np.ndarray:
    on, _, _ = load_ev_3col(ev_path)
    return on


def get_ev_dir(sub: str, ses: str, run: str) -> Path:
    return EV_ROOT / f"sub-{sub}" / f"ses-{ses}" / f"run-{run}"


def load_tr_from_design_fsf(feat: Path) -> Optional[float]:
    fsf = feat / "design.fsf"
    if not fsf.exists():
        return None
    txt = fsf.read_text(errors="ignore")
    m = re.search(r"set\s+fmri\(tr\)\s+([0-9]*\.?[0-9]+)", txt)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def load_vs_timeseries(feat: Path) -> Tuple[np.ndarray, float]:
    func_path = feat / FEAT_FUNC
    mask_path = feat / VS_MASK_REL
    if not func_path.exists():
        raise FileNotFoundError(f"Missing: {func_path}")
    if not mask_path.exists():
        raise FileNotFoundError(f"Missing: {mask_path}")

    img = nib.load(str(func_path))
    data = img.get_fdata()
    tr = float(img.header.get_zooms()[3])

    mimg = nib.load(str(mask_path))
    mdat = mimg.get_fdata()
    mask = mdat > 0.5
    if mask.ndim == 4:
        mask = mask[..., 0]

    ts = data[mask, :].mean(axis=0).astype(float)
    return ts, tr


def find_confoundevs_file(feat: Path) -> Optional[Path]:
    """Try a few common locations for FEAT confound EVs."""
    c1 = feat / CONFOUNDEVS_NAME
    c2 = feat / "confoundevs" / CONFOUNDEVS_NAME
    c3 = feat / "confoundevs.txt"  # redundant, but explicit
    for c in (c1, c2, c3):
        if c.exists():
            return c
    return None


def load_confound_matrix(feat: Path, n_tp: int) -> Tuple[Optional[np.ndarray], int]:
    """
    Load confounds as an (n_tp, k) matrix; drop constant/near-constant columns.
    Returns (X, k_retained). If unavailable or mismatched, returns (None, 0).
    """
    cpath = find_confoundevs_file(feat)
    if cpath is None:
        return None, 0

    try:
        X = np.loadtxt(cpath, ndmin=2).astype(float)
    except Exception:
        return None, 0

    if X.size == 0 or X.shape[0] == 0:
        return None, 0

    # Some confoundevs files can be transposed; enforce rows=time
    if X.shape[0] != n_tp and X.shape[1] == n_tp:
        X = X.T

    if X.shape[0] != n_tp:
        # Length mismatch: do not use
        return None, 0

    # Drop constant/near-constant columns
    keep = []
    for j in range(X.shape[1]):
        col = X[:, j]
        if np.nanstd(col) > 1e-8:
            keep.append(j)
    if not keep:
        return None, 0

    X = X[:, keep]
    return X, X.shape[1]


def regress_out(ts: np.ndarray, X: np.ndarray) -> np.ndarray:
    """
    OLS residualization of ts against confounds X with intercept; add mean(ts) back.
    """
    y = ts.astype(float)
    mu = float(np.mean(y))
    # Add intercept
    X1 = np.column_stack([np.ones((X.shape[0],), dtype=float), X.astype(float)])
    # Solve least squares
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    yhat = X1 @ beta
    resid = y - yhat
    return resid + mu


def extract_event_locked(
    ts: np.ndarray,
    tr: float,
    onsets_s: np.ndarray,
    time_axis: np.ndarray,
    baseline_win: Tuple[float, float],
) -> Curve:
    """
    Event-locked selective averaging using linear interpolation on a fixed time grid.
    Baseline subtract per-event.
    """
    if onsets_s.size == 0:
        mean = np.full_like(time_axis, np.nan, dtype=float)
        sem = np.full_like(time_axis, np.nan, dtype=float)
        return mean, sem, 0

    n_tp = ts.shape[0]
    tpoints = np.arange(n_tp, dtype=float) * tr

    b0, b1 = baseline_win
    bmask = (time_axis >= b0) & (time_axis <= b1)
    if not np.any(bmask):
        raise ValueError("Baseline window does not overlap time axis")

    curves = []
    for o in onsets_s:
        sample_t = o + time_axis
        samp = np.interp(sample_t, tpoints, ts, left=np.nan, right=np.nan)
        b = np.nanmean(samp[bmask])
        curves.append(samp - b)

    mat = np.stack(curves, axis=0)
    mean = np.nanmean(mat, axis=0)
    sem = np.nanstd(mat, axis=0, ddof=1) / math.sqrt(mat.shape[0]) if mat.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem, int(mat.shape[0])


def weighted_mean_and_sem(curves: List[Curve]) -> Curve:
    """
    Weighted aggregation across runs.
    - mean: weighted by n_events
    - sem: SEM across runs (weighted SD / sqrt(n_runs))
    - n: sum of n_events
    """
    if not curves:
        mean = np.full((int(round((T_MAX_S - T_MIN_S) / DT_S)) + 1,), np.nan, dtype=float)
        sem = mean.copy()
        return mean, sem, 0

    mats = np.stack([c[0] for c in curves], axis=0)
    weights = np.array([c[2] for c in curves], dtype=float)
    total_n = int(np.sum(weights))

    if np.sum(weights) <= 0:
        w = np.ones_like(weights) / len(weights)
    else:
        w = weights / np.sum(weights)

    wmean = np.nansum(mats * w[:, None], axis=0)

    if len(curves) == 1:
        sem = np.zeros_like(wmean)
    else:
        wvar = np.nansum(w[:, None] * (mats - wmean[None, :]) ** 2, axis=0)
        wstd = np.sqrt(wvar)
        sem = wstd / math.sqrt(len(curves))

    return wmean, sem, total_n


def markers_for_ant(isi_s: np.ndarray, fb_delay_s: np.ndarray) -> Tuple[List[Tuple[float, str]], List[Tuple[float, float, str]]]:
    """
    Build vertical markers and spans for ANT plots (t=0 is cue offset).
    """
    med_isi = nanmedian_or_nan(isi_s)
    med_fb = nanmedian_or_nan(fb_delay_s)
    vlines: List[Tuple[float, str]] = [
        (-CUE_DUR_S, "Cue onset"),
        (0.0, "Cue offset"),
        (med_isi, "Target (median)"),
        (med_fb, "Feedback (median)"),
        (6.0, "+6s"),
    ]
    spans = [(0.0, 2.0, "Wu ANT window (0–2s)")]
    return vlines, spans


def plot_two_conditions_dual(
    time_axis: np.ndarray,
    A_raw: Curve, labelA: str,
    B_raw: Curve, labelB: str,
    title: str,
    ylabel: str,
    out_png: Path,
    vlines: Optional[List[Tuple[float, str]]] = None,
    spans: Optional[List[Tuple[float, float, str]]] = None,
    A_conf: Optional[Curve] = None,
    B_conf: Optional[Curve] = None,
    conf_label: str = "conf-reg",
) -> None:
    """
    Plot two conditions with:
    - RAW (solid + SEM shading)
    - CONF (dashed; no shading), if provided.
    Colors are matched between raw and conf curves by re-using the raw line's color.
    """
    mA, sA, nA = A_raw
    mB, sB, nB = B_raw

    fig, ax = plt.subplots(figsize=(8, 4.8))

    # Shaded windows
    if spans:
        for x0, x1, lab in spans:
            ax.axvspan(x0, x1, alpha=0.10, label=lab if lab else None)

    # RAW curves
    lineA = ax.plot(time_axis, mA, label=f"{labelA} raw (n={nA})")[0]
    ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)

    lineB = ax.plot(time_axis, mB, label=f"{labelB} raw (n={nB})")[0]
    ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)

    # CONF curves (dashed, same colors)
    if A_conf is not None:
        mAc, _, nAc = A_conf
        if nAc > 0:
            ax.plot(time_axis, mAc, ls="--", color=lineA.get_color(),
                    label=f"{labelA} {conf_label} (n={nAc})")
    if B_conf is not None:
        mBc, _, nBc = B_conf
        if nBc > 0:
            ax.plot(time_axis, mBc, ls="--", color=lineB.get_color(),
                    label=f"{labelB} {conf_label} (n={nBc})")

    # Vertical lines
    if vlines:
        for x, _lab in vlines:
            if np.isfinite(x):
                ax.axvline(x, ls=":", lw=1)
        # annotate labels at top
        y_top = ax.get_ylim()[1]
        for x, lab in vlines:
            if lab and np.isfinite(x):
                ax.text(x, y_top, lab, rotation=90, va="top", ha="right", fontsize=8)

    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def parse_ids(feat_path: Path) -> Tuple[str, str, str, str]:
    """
    Parse sub/ses/run/echo from a FEAT path. This expects patterns somewhere in the path like:
      sub-XXX / ses-YYY / ... run-<N> ... (fmriprep|tedana|single-echo|multi-echo).feat

    If your echo labels differ, adjust the regex below.
    """
    p = str(feat_path)
    msub = re.search(r"sub-([A-Za-z0-9]+)", p)
    mses = re.search(r"ses-([A-Za-z0-9]+)", p)
    mrun = re.search(r"run-([0-9]+)", p)
    # echo label: common cases
    mecho = re.search(r"(fmriprep|tedana|single-echo|multi-echo)\.feat", p)
    if not (msub and mses and mrun and mecho):
        raise ValueError(f"Could not parse ids from: {feat_path}")
    return msub.group(1), mses.group(1), mrun.group(1), mecho.group(1)


def process_one_feat(feat: Path, time_axis: np.ndarray) -> RunResult:
    sub, ses, run, echo = parse_ids(feat)

    ts_raw, tr_header = load_vs_timeseries(feat)
    n_tp = ts_raw.shape[0]

    tr_design = load_tr_from_design_fsf(feat)
    tr_used = float(tr_header)
    if tr_design is not None and abs(tr_design - tr_header) > TR_MISMATCH_WARN_EPS:
        print(f"[WARN] TR mismatch for {feat}: header={tr_header:.6f}s, design.fsf={tr_design:.6f}s")
    if tr_design is not None and PREFER_TR_FROM_DESIGN:
        tr_used = float(tr_design)

    # Confound regression (optional)
    confounds_used = 0
    n_confounds = 0
    ts_conf = None
    if COMPUTE_CONFOUND_VARIANT:
        X, k = load_confound_matrix(feat, n_tp)
        if X is not None and k > 0:
            ts_conf = regress_out(ts_raw, X)
            confounds_used = 1
            n_confounds = k

    # PSC time series
    ts_psc_raw = to_psc(ts_raw)
    ts_psc_conf = to_psc(ts_conf) if (ts_conf is not None) else None

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)

    # Anticipation EVs: these are cue-onset with duration; shift to cue-offset anchor
    ant_R_on, ant_R_dur, _ = load_ev_3col(ev_dir / "_anticipation_reward.txt")
    ant_N_on, ant_N_dur, _ = load_ev_3col(ev_dir / "_anticipation_neutral.txt")
    ant_R = ant_R_on + CUE_DUR_S
    ant_N = ant_N_on + CUE_DUR_S

    # ISI diagnostic (target onset relative to cue-offset) ~ dur - cue_dur
    isi_s = np.concatenate([ant_R_dur, ant_N_dur]) - CUE_DUR_S
    isi_s = isi_s[np.isfinite(isi_s)]

    # Feedback EVs (for plotting pooled by valence)
    fb_pos = np.sort(np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_positive_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_positive_neutral.txt"),
    ]))
    fb_neg = np.sort(np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_negative_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_negative_neutral.txt"),
    ]))
    fb_all = np.sort(np.concatenate([fb_pos, fb_neg]))

    # Feedback timing diagnostic: next feedback after each cue-offset
    cueoffs = np.sort(np.concatenate([ant_R, ant_N]))
    if cueoffs.size == 0 or fb_all.size == 0:
        fb_delay_s = np.array([], dtype=float)
    else:
        idx_next = np.searchsorted(fb_all, cueoffs, side="right")
        valid = idx_next < fb_all.size
        fb_delay_s = fb_all[idx_next[valid]] - cueoffs[valid]

    # Event-locked curves (RAW)
    ant_psc_raw = {
        "Reward": extract_event_locked(ts_psc_raw, tr_used, ant_R, time_axis, BASELINE_WIN),
        "Neutral": extract_event_locked(ts_psc_raw, tr_used, ant_N, time_axis, BASELINE_WIN),
    }
    fb_psc_raw = {
        "Positive": extract_event_locked(ts_psc_raw, tr_used, fb_pos, time_axis, BASELINE_WIN),
        "Negative": extract_event_locked(ts_psc_raw, tr_used, fb_neg, time_axis, BASELINE_WIN),
    }

    # Event-locked curves (CONF), only if confounds were actually used
    ant_psc_conf: Dict[str, Curve] = {}
    fb_psc_conf: Dict[str, Curve] = {}
    if confounds_used == 1 and ts_psc_conf is not None:
        ant_psc_conf = {
            "Reward": extract_event_locked(ts_psc_conf, tr_used, ant_R, time_axis, BASELINE_WIN),
            "Neutral": extract_event_locked(ts_psc_conf, tr_used, ant_N, time_axis, BASELINE_WIN),
        }
        fb_psc_conf = {
            "Positive": extract_event_locked(ts_psc_conf, tr_used, fb_pos, time_axis, BASELINE_WIN),
            "Negative": extract_event_locked(ts_psc_conf, tr_used, fb_neg, time_axis, BASELINE_WIN),
        }

    return RunResult(
        sub=sub, ses=ses, run=run, echo=echo,
        time_axis=time_axis,
        tr_header=float(tr_header), tr_design=tr_design, tr_used=float(tr_used),
        isi_s=isi_s, fb_delay_s=fb_delay_s,
        confounds_used=confounds_used, n_confounds=n_confounds,
        ant_psc_raw=ant_psc_raw, fb_psc_raw=fb_psc_raw,
        ant_psc_conf=ant_psc_conf, fb_psc_conf=fb_psc_conf,
    )


def find_feat_dirs() -> List[Path]:
    if FEAT_LIST_PATH.exists():
        feats: List[Path] = []
        for line in FEAT_LIST_PATH.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            feats.append(Path(line))
        return feats
    # Default discovery
    return sorted(FSL_DERIV.glob("sub-*/ses-*/L1_*.feat"))


def aggregate_subject(results: List[RunResult], subject: str, time_axis: np.ndarray) -> None:
    subj_runs = [r for r in results if r.sub == subject]
    if not subj_runs:
        return

    # group by echo label
    by_echo: Dict[str, List[RunResult]] = {}
    for r in subj_runs:
        by_echo.setdefault(r.echo, []).append(r)

    for echo, runs in by_echo.items():
        # markers pooled across runs
        isi_all = np.concatenate([r.isi_s for r in runs if r.isi_s.size > 0]) if runs else np.array([], dtype=float)
        fb_all = np.concatenate([r.fb_delay_s for r in runs if r.fb_delay_s.size > 0]) if runs else np.array([], dtype=float)
        vlines, spans = markers_for_ant(isi_all, fb_all)

        # RAW aggregates (always)
        antR_raw = weighted_mean_and_sem([r.ant_psc_raw["Reward"] for r in runs])
        antN_raw = weighted_mean_and_sem([r.ant_psc_raw["Neutral"] for r in runs])
        fbP_raw = weighted_mean_and_sem([r.fb_psc_raw["Positive"] for r in runs])
        fbN_raw = weighted_mean_and_sem([r.fb_psc_raw["Negative"] for r in runs])

        # CONF aggregates (only runs where confounds_used==1)
        runs_conf = [r for r in runs if r.confounds_used == 1]
        antR_conf = weighted_mean_and_sem([r.ant_psc_conf["Reward"] for r in runs_conf]) if runs_conf else None
        antN_conf = weighted_mean_and_sem([r.ant_psc_conf["Neutral"] for r in runs_conf]) if runs_conf else None
        fbP_conf = weighted_mean_and_sem([r.fb_psc_conf["Positive"] for r in runs_conf]) if runs_conf else None
        fbN_conf = weighted_mean_and_sem([r.fb_psc_conf["Negative"] for r in runs_conf]) if runs_conf else None

        out_dir = OUT_DIR / "subjects" / f"sub-{subject}" / echo

        plot_two_conditions_dual(
            time_axis,
            antR_raw, "Reward",
            antN_raw, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=out_dir / "anticipation_psc.png",
            vlines=vlines,
            spans=spans,
            A_conf=antR_conf,
            B_conf=antN_conf,
            conf_label="conf-reg",
        )

        plot_two_conditions_dual(
            time_axis,
            fbP_raw, "Feedback +",
            fbN_raw, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=out_dir / "feedback_psc.png",
            vlines=[(0.0, "Feedback onset"), (6.0, "+6s")],
            spans=None,
            A_conf=fbP_conf,
            B_conf=fbN_conf,
            conf_label="conf-reg",
        )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    time_axis = np.arange(T_MIN_S, T_MAX_S + 1e-9, DT_S)

    feats = find_feat_dirs()
    if not feats:
        raise SystemExit("No FEAT directories found. Check FSL_DERIV or feat_list.txt.")

    results: List[RunResult] = []
    for feat in feats:
        try:
            res = process_one_feat(feat, time_axis)
        except Exception as e:
            print(f"[ERROR] Skipping {feat} due to: {e}")
            continue
        results.append(res)

    # Subject plots
    for sub in sorted(set(r.sub for r in results)):
        aggregate_subject(results, sub, time_axis)

    # Summary TSV (run-level)
    base_cols = [
        "sub", "ses", "run", "echo",
        "tr_header", "tr_design", "tr_used",
        "confounds_used", "n_confounds",
        "med_isi", "med_fbdelay", "n_isi", "n_fbdelay",
    ]

    with open(SUMMARY_TSV, "w") as f:
        f.write("\t".join(base_cols) + "\n")
        for r in results:
            med_isi = nanmedian_or_nan(r.isi_s)
            med_fb = nanmedian_or_nan(r.fb_delay_s)
            row = [
                r.sub, r.ses, r.run, r.echo,
                f"{r.tr_header:.6f}",
                "" if r.tr_design is None else f"{r.tr_design:.6f}",
                f"{r.tr_used:.6f}",
                str(r.confounds_used),
                str(r.n_confounds),
                "" if np.isnan(med_isi) else f"{med_isi:.6f}",
                "" if np.isnan(med_fb) else f"{med_fb:.6f}",
                str(int(r.isi_s.size)),
                str(int(r.fb_delay_s.size)),
            ]
            f.write("\t".join(row) + "\n")

    # TR quick summary
    trs = np.array([r.tr_header for r in results], dtype=float)
    if trs.size > 0:
        vals, counts = np.unique(np.round(trs, 6), return_counts=True)
        pairs = ", ".join([f"{v:.6f} ({c})" for v, c in zip(vals, counts)])
        print(f"TR (header) values across processed runs: {pairs}")

    print(f"Done.\n  Plots: {OUT_DIR}\n  Summary: {SUMMARY_TSV}")


if __name__ == "__main__":
    main()
