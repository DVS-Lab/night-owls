#!/usr/bin/env python3
"""
VS time‑course extraction, plotting, and discrete 4th‑TR summaries
------------------------------------------------------------------
Assumptions
- Run from <rootdir>/code (this file lives in that directory).
- All data are already in MNI152NLin6Asym space.
- Input FEAT directories are listed one per line in ./feat_paths.txt (absolute
  paths), exactly as provided by the user. The script will ignore commented
  lines (starting with '#') and blank lines.
- EV files live at:
    <rootdir>/derivatives/fsl/EVFiles/sub-XXX/ses-YY/mid/run-ZZ/
  with the following filenames:
    _anticipation_reward.txt
    _anticipation_neutral.txt
    _feedback_positive_reward.txt
    _feedback_negative_reward.txt
    _feedback_positive_neutral.txt
    _feedback_negative_neutral.txt
- Ventral Striatum (VS) mask lives at:
    <rootdir>/masks/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz

What it does
- Extracts the VS mean time series from each FEAT's filtered_func_data.nii.gz.
- Computes two normalized versions per run:
    • PSC   : ((ts / ts.mean()) - 1) * 100
    • Z     : (ts - ts.mean()) / ts.std(ddof=1)
- Builds event‑locked peri‑stimulus windows with linear interpolation
  (fractional onsets are handled) for two analysis families:
    • ANTICIPATION: Reward vs Neutral (cue‑locked)
    • FEEDBACK (valence): Positive vs Negative (feedback‑locked, pooling across
      incentive conditions)
- Produces per‑run plots and subject‑level aggregated plots split by echo
  (single‑echo vs multi‑echo) for PSC and Z (four figures per subject per family).
- Additionally, writes **discrete 4th‑TR‑after‑onset summaries** (no interpolation)
  for **six conditions** (ANT_R, ANT_N, FB_POS_R, FB_NEG_R, FB_POS_N, FB_NEG_N)
  with PSC and Z, per run, aggregated into a single CSV/TSV.

Outputs
  <rootdir>/derivatives/extractions/timecourses/
    ├── runs/  (per FEAT plots + per‑run TP TSV)
    └── subjects/  (per subject, echo‑split plots)
  <rootdir>/derivatives/extractions/summary_at_4thTR.tsv  (and .csv)

Dependencies: nibabel, numpy, matplotlib
"""
from __future__ import annotations

import re
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# --------------------------- Fixed project structure ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
MASKS_DIR  = ROOT_DIR / "masks"
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-unsmoothed"
SUMMARY_DIR= ROOT_DIR / "derivatives" / "extractions"
EV_BASE    = FSL_DERIV / "EVFiles"  # EVs are arranged by sub/ses/task/run below

VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"

# Input list of FEAT directories (one per line, absolute). You provided this.
FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-unsmoothed.txt"

# --------------------------- Analysis parameters ------------------------------
TASK         = "mid"
TR_HARDCODE  = 1.615  # sec (we also check header and warn on mismatch)
K_AFTER      = 3      # 4th TR after onset -> index at floor((onset + 3*TR)/TR)
TMIN         = -4.0   # sec relative to event onset (pre‑stim baseline window)
TMAX         = 16.0   # sec after event onset
# Vertical reference lines to annotate typical MID timings (sec from event)
VERT_LINES   = [0.0, 4.0, 6.0]

# Plot style
plt.rcParams.update({
    "figure.dpi": 120,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 12,
})

# ------------------------------- Data classes ---------------------------------
@dataclass
class RunResult:
    feat_path: Path
    sub: str
    ses: str
    run: str
    echo: str  # 'single-echo' or 'multi-echo'
    time_axis: np.ndarray
    # Each dict maps condition label -> (mean, sem, n_trials)
    ant_psc: Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    ant_z:   Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_psc:  Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_z:    Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    # Discrete 4th‑TR summaries (no interpolation)
    tp_means_psc: Dict[str, float | None]
    tp_means_z:   Dict[str, float | None]
    tp_counts:    Dict[str, int]

    # Optional confound-regressed versions (same structure as above).
    # Populated only when confoundevs.txt exists AND matches the fMRI time series length.
    ant_psc_cr: Dict[str, Tuple[np.ndarray, np.ndarray, int]] | None = None
    ant_z_cr:   Dict[str, Tuple[np.ndarray, np.ndarray, int]] | None = None
    fb_psc_cr:  Dict[str, Tuple[np.ndarray, np.ndarray, int]] | None = None
    fb_z_cr:    Dict[str, Tuple[np.ndarray, np.ndarray, int]] | None = None
    confounds_applied: bool = False


# ------------------------------ Helper functions ------------------------------

def parse_sub_ses_run_from_feat(feat: Path) -> Tuple[str, str, str]:
    """Parse sub, ses, run from FEAT path based on your naming scheme."""
    m_sub = re.search(r"sub-(\d+)", feat.as_posix())
    m_ses = re.search(r"ses-(\d+)", feat.as_posix())
    m_run = re.search(r"run-(\d+)", feat.as_posix())
    if not (m_sub and m_ses and m_run):
        raise ValueError(f"Could not parse sub/ses/run from: {feat}")
    return m_sub.group(1), m_ses.group(1), m_run.group(1)


def echo_from_feat(feat: Path) -> str:
    if "single-echo" in feat.name or "single-echo" in feat.as_posix():
        return "single-echo"
    if "multi-echo" in feat.name or "multi-echo" in feat.as_posix():
        return "multi-echo"
    return "unknown-echo"


def get_ev_dir(sub: str, ses: str, run: str) -> Path:
    return EV_BASE / f"sub-{sub}" / f"ses-{ses}" / TASK / f"run-{run}"


def load_ev(ev_path: Path) -> np.ndarray:
    """Load a 3‑column EV (onset, duration, amplitude) -> onsets (seconds)."""
    if not ev_path.exists():
        return np.array([], dtype=float)
    try:
        arr = np.loadtxt(ev_path, ndmin=2)
    except Exception:
        return np.array([], dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float)
    return np.asarray(arr[:, 0], dtype=float)


def warn_tr(nifti_img: nib.Nifti1Image, feat: Path) -> float:
    hdr_tr = float(nifti_img.header.get_zooms()[3])
    if abs(hdr_tr - TR_HARDCODE) > 1e-3:
        print(f"[WARN] TR mismatch in {feat}: header={hdr_tr:.6f}s vs hardcoded={TR_HARDCODE:.6f}s. Using hardcoded.")
    return TR_HARDCODE


def load_vs_timeseries(feat: Path) -> Tuple[np.ndarray, float]:
    img_path = feat / "filtered_func_data.nii.gz"
    if not img_path.exists():
        raise FileNotFoundError(f"Missing filtered_func_data.nii.gz in {feat}")
    if not VS_MNI.exists():
        raise FileNotFoundError(f"Missing VS mask at {VS_MNI}")
    img = nib.load(str(img_path))
    mask_img = nib.load(str(VS_MNI))

    data = img.get_fdata()  # X×Y×Z×T
    mask = mask_img.get_fdata() > 0
    if data.shape[:3] != mask.shape:
        raise ValueError(f"Mask dims {mask.shape} do not match data dims {data.shape[:3]} for {feat}")

    vox_ts = data[mask, :]  # V×T
    vs_ts = np.nanmean(vox_ts, axis=0)  # T
    tr = warn_tr(img, feat)
    return vs_ts.astype(float), tr



def load_confound_evs(feat: Path, n_tp: int) -> np.ndarray | None:
    """Load FEAT confounds (confoundevs.txt) if available and well-formed.

    Returns
    -------
    conf : (n_tp, n_conf) array or None
        None if file missing/empty or if row-count doesn't match n_tp.
    """
    f = feat / "confoundevs.txt"
    if not f.exists():
        return None
    try:
        conf = np.loadtxt(f)
    except Exception:
        return None

    if conf.size == 0:
        return None

    conf = np.atleast_2d(conf)
    # If loadtxt returns shape (n_conf,) for single-row files, fix to (1, n_conf)
    if conf.shape[0] == 1 and n_tp != 1 and conf.shape[1] == n_tp:
        conf = conf.T

    if conf.shape[0] != n_tp:
        print(f"[WARN] confoundevs rows ({conf.shape[0]}) != n_tp ({n_tp}) in: {f}")
        return None

    # Replace NaNs/Infs defensively
    conf = np.nan_to_num(conf, nan=0.0, posinf=0.0, neginf=0.0)
    return conf


def regress_out_confounds(ts: np.ndarray, conf: np.ndarray) -> np.ndarray:
    """Regress confounds (plus intercept) from a 1D time series and preserve its mean."""
    y = np.asarray(ts, dtype=float)
    if y.ndim != 1:
        raise ValueError("ts must be 1D")
    if conf.ndim != 2 or conf.shape[0] != y.size:
        raise ValueError("confounds must be (n_tp, n_conf)")

    # Design matrix: intercept + confounds
    X = np.column_stack([np.ones(y.size), conf])

    # Least-squares fit
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ beta
    resid = y - y_hat

    # Add back original mean for interpretability (keeps PSC scale comparable)
    mu = np.nanmean(y)
    return resid + mu


def to_psc(ts: np.ndarray) -> np.ndarray:
    mean = np.nanmean(ts)
    return ((ts / mean) - 1.0) * 100.0


def to_z(ts: np.ndarray) -> np.ndarray:
    mu = np.nanmean(ts)
    sd = np.nanstd(ts, ddof=1)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(ts)
    return (ts - mu) / sd


def build_time_axis(tmin: float, tmax: float, tr: float) -> np.ndarray:
    n = int(np.floor((tmax - tmin) / tr)) + 1
    return (tmin + np.arange(n) * tr).astype(float)


def sample_windows(ts: np.ndarray, onsets_s: np.ndarray, tr: float,
                   tmin: float, tmax: float) -> np.ndarray:
    """Return array of shape (n_trials, n_time) sampled by linear interpolation."""
    if onsets_s.size == 0:
        return np.empty((0, 0))
    t_series = np.arange(ts.size) * tr
    t_axis = build_time_axis(tmin, tmax, tr)
    windows = []
    for onset in onsets_s:
        t_abs = onset + t_axis
        vals = np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)
        windows.append(vals)
    return np.vstack(windows)


def mean_and_sem(windows: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    if windows.size == 0:
        return np.array([]), np.array([]), 0
    valid = ~np.all(~np.isfinite(windows), axis=0)
    w = np.where(valid, windows, np.nan)
    mean = np.nanmean(w, axis=0)
    n = np.sum(np.isfinite(w), axis=0)
    sem = np.nanstd(w, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))
    n_trials = int(np.sum(~np.all(np.isnan(w), axis=1)))
    return mean, sem, n_trials


def fourth_tr_indices(onsets: np.ndarray, tr: float, n_vols: int) -> np.ndarray:
    """Indices for the volume whose start lies in [onset + 3*TR, onset + 4*TR)."""
    if onsets.size == 0:
        return np.empty((0,), dtype=int)
    idx = np.floor((onsets + K_AFTER * tr) / tr).astype(int)
    return idx[(idx >= 0) & (idx < n_vols)]


def plot_two_conditions(time_axis: np.ndarray,
                        condA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
                        condB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
                        title: str, ylabel: str, out_png: Path) -> None:
    mA, sA, nA = condA
    mB, sB, nB = condB
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(time_axis, mA, label=f"{labelA} (n={nA})")
    ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)
    ax.plot(time_axis, mB, label=f"{labelB} (n={nB})")
    ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)
    for v in VERT_LINES:
        ax.axvline(v, ls=":", lw=1)
    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_two_conditions_dualconf(
    time_axis: np.ndarray,
    rawA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
    rawB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
    crA: Tuple[np.ndarray, np.ndarray, int] | None,
    crB: Tuple[np.ndarray, np.ndarray, int] | None,
    title: str,
    ylabel: str,
    out_png: Path,
    vlines: List[float] | None = None,
):
    """Plot raw (solid + band) and confound-regressed (dashed) overlays for two conditions.

    If crA/crB is None, this falls back to the raw-only visualization.
    """
    out_png.parent.mkdir(parents=True, exist_ok=True)

    mA, seA, nA = rawA
    mB, seB, nB = rawB

    # If either condition is empty, make a minimal plot (prevents hard crashes in edge cases)
    if mA.size == 0 or mB.size == 0:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.set_title(title)
        ax.set_xlabel("Time from event (s)")
        ax.set_ylabel(ylabel)
        ax.text(0.5, 0.5, "No trials for one or both conditions", ha="center", va="center", transform=ax.transAxes)
        if vlines:
            for x in vlines:
                ax.axvline(x, linestyle=":", alpha=0.6)
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    # Raw (solid + band)
    (lineA,) = ax.plot(time_axis, mA, label=f"{labelA} (raw; n={nA})")
    ax.fill_between(time_axis, mA - seA, mA + seA, alpha=0.2)
    (lineB,) = ax.plot(time_axis, mB, label=f"{labelB} (raw; n={nB})")
    ax.fill_between(time_axis, mB - seB, mB + seB, alpha=0.2)

    # Confound-regressed (dashed; no band to avoid clutter)
    if crA is not None and crB is not None:
        mA2, _, nA2 = crA
        mB2, _, nB2 = crB
        if mA2.size == mA.size:
            ax.plot(time_axis, mA2, linestyle="--", color=lineA.get_color(), label=f"{labelA} (confound-reg; n={nA2})")
        if mB2.size == mB.size:
            ax.plot(time_axis, mB2, linestyle="--", color=lineB.get_color(), label=f"{labelB} (confound-reg; n={nB2})")

    if vlines:
        for x in vlines:
            ax.axvline(x, linestyle=":", alpha=0.6)

    ax.set_title(title)
    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

# ------------------------------- Core pipeline --------------------------------

def process_one_feat(feat: Path) -> RunResult | None:
    try:
        sub, ses, run = parse_sub_ses_run_from_feat(feat)
        echo = echo_from_feat(feat)
        ts_raw, tr = load_vs_timeseries(feat)
    except Exception as e:
        print(f"[WARN] Skipping FEAT due to error: {feat}\n  -> {e}")
        return None

    ts_psc = to_psc(ts_raw)
    ts_z   = to_z(ts_raw)
    T = ts_raw.size

    # Optional confound regression (uses FEAT confoundevs.txt when present)
    confounds = load_confound_evs(feat, T)
    confounds_applied = False
    ts_psc_cr = ts_z_cr = None
    if confounds is not None:
        try:
            ts_clean = regress_out_confounds(ts_raw, confounds)
            ts_psc_cr = to_psc(ts_clean)
            ts_z_cr   = to_z(ts_clean)
            confounds_applied = True
        except Exception as e:
            print(f"[WARN] Confound regression failed in {feat}: {e}")
            confounds_applied = False
            ts_psc_cr = ts_z_cr = None

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)
    ant_R = load_ev(ev_dir / "_anticipation_reward.txt")
    ant_N = load_ev(ev_dir / "_anticipation_neutral.txt")

    # Feedback (pooled by valence for plotting)
    fb_pos = np.sort(np.concatenate([
        load_ev(ev_dir / "_feedback_positive_reward.txt"),
        load_ev(ev_dir / "_feedback_positive_neutral.txt"),
    ]))
    fb_neg = np.sort(np.concatenate([
        load_ev(ev_dir / "_feedback_negative_reward.txt"),
        load_ev(ev_dir / "_feedback_negative_neutral.txt"),
    ]))

    # Feedback (separate conditions for discrete TP summaries)
    fb_PR = load_ev(ev_dir / "_feedback_positive_reward.txt")
    fb_NR = load_ev(ev_dir / "_feedback_negative_reward.txt")
    fb_PN = load_ev(ev_dir / "_feedback_positive_neutral.txt")
    fb_NN = load_ev(ev_dir / "_feedback_negative_neutral.txt")

    # Build time axis for peri‑event plots
    t_axis = build_time_axis(TMIN, TMAX, tr)

    # ANT windows
    R_psc_w = sample_windows(ts_psc, ant_R, tr, TMIN, TMAX)
    N_psc_w = sample_windows(ts_psc, ant_N, tr, TMIN, TMAX)
    R_z_w   = sample_windows(ts_z,   ant_R, tr, TMIN, TMAX)
    N_z_w   = sample_windows(ts_z,   ant_N, tr, TMIN, TMAX)

    ant_psc = {
        "Reward": mean_and_sem(R_psc_w),
        "Neutral": mean_and_sem(N_psc_w),
    }
    ant_z = {
        "Reward": mean_and_sem(R_z_w),
        "Neutral": mean_and_sem(N_z_w),
    }

    # Optional confound-regressed ANT windows
    ant_psc_cr = ant_z_cr = None
    if confounds_applied and ts_psc_cr is not None and ts_z_cr is not None:
        R_psc_w_cr = sample_windows(ts_psc_cr, ant_R, tr, TMIN, TMAX)
        N_psc_w_cr = sample_windows(ts_psc_cr, ant_N, tr, TMIN, TMAX)
        R_z_w_cr   = sample_windows(ts_z_cr,   ant_R, tr, TMIN, TMAX)
        N_z_w_cr   = sample_windows(ts_z_cr,   ant_N, tr, TMIN, TMAX)

        ant_psc_cr = {
            "Reward": mean_and_sem(R_psc_w_cr),
            "Neutral": mean_and_sem(N_psc_w_cr),
        }
        ant_z_cr = {
            "Reward": mean_and_sem(R_z_w_cr),
            "Neutral": mean_and_sem(N_z_w_cr),
        }

    # FB windows (valence pooled)
    P_psc_w = sample_windows(ts_psc, fb_pos, tr, TMIN, TMAX)
    G_psc_w = sample_windows(ts_psc, fb_neg, tr, TMIN, TMAX)
    P_z_w   = sample_windows(ts_z,   fb_pos, tr, TMIN, TMAX)
    G_z_w   = sample_windows(ts_z,   fb_neg, tr, TMIN, TMAX)

    fb_psc = {
        "Positive": mean_and_sem(P_psc_w),
        "Negative": mean_and_sem(G_psc_w),
    }
    fb_z = {
        "Positive": mean_and_sem(P_z_w),
        "Negative": mean_and_sem(G_z_w),
    }

    # Optional confound-regressed FB windows (valence pooled)
    fb_psc_cr = fb_z_cr = None
    if confounds_applied and ts_psc_cr is not None and ts_z_cr is not None:
        P_psc_w_cr = sample_windows(ts_psc_cr, fb_pos, tr, TMIN, TMAX)
        G_psc_w_cr = sample_windows(ts_psc_cr, fb_neg, tr, TMIN, TMAX)
        P_z_w_cr   = sample_windows(ts_z_cr,   fb_pos, tr, TMIN, TMAX)
        G_z_w_cr   = sample_windows(ts_z_cr,   fb_neg, tr, TMIN, TMAX)

        fb_psc_cr = {
            "Positive": mean_and_sem(P_psc_w_cr),
            "Negative": mean_and_sem(G_psc_w_cr),
        }
        fb_z_cr = {
            "Positive": mean_and_sem(P_z_w_cr),
            "Negative": mean_and_sem(G_z_w_cr),
        }

    # ---------------- Discrete 4th‑TR summaries (no interpolation) -------------
    def m(arr: np.ndarray) -> float | None:
        return float(np.nanmean(arr)) if arr.size else None

    def vals_at_tp(ts: np.ndarray, onsets: np.ndarray) -> np.ndarray:
        idx = fourth_tr_indices(onsets, tr, T)
        return ts[idx]

    # Anticipation (Reward/Neutral)
    ant_R_psc_tp = vals_at_tp(ts_psc, ant_R)
    ant_N_psc_tp = vals_at_tp(ts_psc, ant_N)
    ant_R_z_tp   = vals_at_tp(ts_z,   ant_R)
    ant_N_z_tp   = vals_at_tp(ts_z,   ant_N)

    # Feedback (six conditions)
    fb_PR_psc_tp = vals_at_tp(ts_psc, fb_PR)
    fb_NR_psc_tp = vals_at_tp(ts_psc, fb_NR)
    fb_PN_psc_tp = vals_at_tp(ts_psc, fb_PN)
    fb_NN_psc_tp = vals_at_tp(ts_psc, fb_NN)

    fb_PR_z_tp   = vals_at_tp(ts_z,   fb_PR)
    fb_NR_z_tp   = vals_at_tp(ts_z,   fb_NR)
    fb_PN_z_tp   = vals_at_tp(ts_z,   fb_PN)
    fb_NN_z_tp   = vals_at_tp(ts_z,   fb_NN)

    tp_means_psc = {
        "ANT_REWARD":  m(ant_R_psc_tp),
        "ANT_NEUTRAL": m(ant_N_psc_tp),
        "FB_POS_REWARD":  m(fb_PR_psc_tp),
        "FB_NEG_REWARD":  m(fb_NR_psc_tp),
        "FB_POS_NEUTRAL": m(fb_PN_psc_tp),
        "FB_NEG_NEUTRAL": m(fb_NN_psc_tp),
    }
    tp_means_z = {
        "ANT_REWARD":  m(ant_R_z_tp),
        "ANT_NEUTRAL": m(ant_N_z_tp),
        "FB_POS_REWARD":  m(fb_PR_z_tp),
        "FB_NEG_REWARD":  m(fb_NR_z_tp),
        "FB_POS_NEUTRAL": m(fb_PN_z_tp),
        "FB_NEG_NEUTRAL": m(fb_NN_z_tp),
    }
    tp_counts = {
        "ANT_REWARD":  int(ant_R_psc_tp.size),
        "ANT_NEUTRAL": int(ant_N_psc_tp.size),
        "FB_POS_REWARD":  int(fb_PR_psc_tp.size),
        "FB_NEG_REWARD":  int(fb_NR_psc_tp.size),
        "FB_POS_NEUTRAL": int(fb_PN_psc_tp.size),
        "FB_NEG_NEUTRAL": int(fb_NN_psc_tp.size),
    }

    return RunResult(
        feat_path=feat,
        sub=sub, ses=ses, run=run, echo=echo,
        time_axis=t_axis,
        ant_psc=ant_psc, ant_z=ant_z,
        fb_psc=fb_psc, fb_z=fb_z,
        tp_means_psc=tp_means_psc,
        tp_means_z=tp_means_z,
        tp_counts=tp_counts,
        ant_psc_cr=ant_psc_cr,
        ant_z_cr=ant_z_cr,
        fb_psc_cr=fb_psc_cr,
        fb_z_cr=fb_z_cr,
        confounds_applied=confounds_applied,
    )



def save_run_plots(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

    # Anticipation: PSC and Z
    if res.confounds_applied and res.ant_psc_cr is not None and res.ant_z_cr is not None:
        plot_two_conditions_dualconf(
            res.time_axis,
            res.ant_psc["Reward"], "Reward",
            res.ant_psc["Neutral"], "Neutral",
            res.ant_psc_cr["Reward"],
            res.ant_psc_cr["Neutral"],
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=run_out / "anticipation_psc.png",
        )
        plot_two_conditions_dualconf(
            res.time_axis,
            res.ant_z["Reward"], "Reward",
            res.ant_z["Neutral"], "Neutral",
            res.ant_z_cr["Reward"],
            res.ant_z_cr["Neutral"],
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, Z)",
            ylabel="Z-scored signal",
            out_png=run_out / "anticipation_z.png",
        )
    else:
        plot_two_conditions(
            res.time_axis,
            res.ant_psc["Reward"], "Reward",
            res.ant_psc["Neutral"], "Neutral",
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=run_out / "anticipation_psc.png",
        )
        plot_two_conditions(
            res.time_axis,
            res.ant_z["Reward"], "Reward",
            res.ant_z["Neutral"], "Neutral",
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, Z)",
            ylabel="Z-scored signal",
            out_png=run_out / "anticipation_z.png",
        )

    # Feedback (valence): PSC and Z
    if res.confounds_applied and res.fb_psc_cr is not None and res.fb_z_cr is not None:
        plot_two_conditions_dualconf(
            res.time_axis,
            res.fb_psc["Positive"], "Feedback +",
            res.fb_psc["Negative"], "Feedback −",
            res.fb_psc_cr["Positive"],
            res.fb_psc_cr["Negative"],
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=run_out / "feedback_psc.png",
        )
        plot_two_conditions_dualconf(
            res.time_axis,
            res.fb_z["Positive"], "Feedback +",
            res.fb_z["Negative"], "Feedback −",
            res.fb_z_cr["Positive"],
            res.fb_z_cr["Negative"],
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, Z)",
            ylabel="Z-scored signal",
            out_png=run_out / "feedback_z.png",
        )
    else:
        plot_two_conditions(
            res.time_axis,
            res.fb_psc["Positive"], "Feedback +",
            res.fb_psc["Negative"], "Feedback −",
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=run_out / "feedback_psc.png",
        )
        plot_two_conditions(
            res.time_axis,
            res.fb_z["Positive"], "Feedback +",
            res.fb_z["Negative"], "Feedback −",
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, Z)",
            ylabel="Z-scored signal",
            out_png=run_out / "feedback_z.png",
        )

    # Save run-level summary: points at +6s (FSL-aligned) and +4th TR (classic MID)
    # This is intentionally raw-only (confound-regressed curves are for visualization/QC).
    summ_path = run_out / "summary_points.tsv"
    summ_path.parent.mkdir(parents=True, exist_ok=True)

    with summ_path.open("w") as f:
        f.write("sub    ses run echo    ")
        f.write("ant_R_psc  ant_N_psc   ant_R_z ant_N_z ")
        f.write("fb_pos_psc fb_neg_psc  fb_pos_z    fb_neg_z    ")
        f.write("tp_total   tp_good\n")
        f.write(f"{res.sub}\t{res.ses}\t{res.run}\t{res.echo}\t")
        f.write(f"{res.points_psc['ant_R']:.6f}\t{res.points_psc['ant_N']:.6f}\t")
        f.write(f"{res.points_z['ant_R']:.6f}\t{res.points_z['ant_N']:.6f}\t")
        f.write(f"{res.points_psc['fb_pos']:.6f}\t{res.points_psc['fb_neg']:.6f}\t")
        f.write(f"{res.points_z['fb_pos']:.6f}\t{res.points_z['fb_neg']:.6f}\t")
        f.write(f"{res.tp_counts['total']}\t{res.tp_counts['good']}\n")



def aggregate_subject(results: List[RunResult], sub: str) -> None:
    """Aggregate across runs within each echo type, and emit subject-level plots + TSV."""
    by_echo: Dict[str, List[RunResult]] = {}
    for res in results:
        by_echo.setdefault(res.echo, []).append(res)

    for echo, runs in by_echo.items():
        outdir = OUT_TC_DIR / "subjects" / f"sub-{sub}" / echo
        outdir.mkdir(parents=True, exist_ok=True)

        # Helper: stack run-level mean curves (each run contributes one curve)
        def stack(items: List[Tuple[np.ndarray, np.ndarray, int]]) -> Tuple[np.ndarray, np.ndarray, int]:
            return weighted_mean_and_sem(items)

        # ----- ANT (raw)
        ant_psc_R = stack([r.ant_psc["Reward"] for r in runs])
        ant_psc_N = stack([r.ant_psc["Neutral"] for r in runs])
        ant_z_R   = stack([r.ant_z["Reward"] for r in runs])
        ant_z_N   = stack([r.ant_z["Neutral"] for r in runs])

        # Optional ANT (confound-regressed): only include runs where it was actually applied
        runs_cr = [r for r in runs if r.confounds_applied and r.ant_psc_cr is not None and r.ant_z_cr is not None]
        ant_psc_R_cr = ant_psc_N_cr = ant_z_R_cr = ant_z_N_cr = None
        if len(runs_cr) > 0:
            ant_psc_R_cr = stack([r.ant_psc_cr["Reward"] for r in runs_cr])  # type: ignore[index]
            ant_psc_N_cr = stack([r.ant_psc_cr["Neutral"] for r in runs_cr])  # type: ignore[index]
            ant_z_R_cr   = stack([r.ant_z_cr["Reward"] for r in runs_cr])  # type: ignore[index]
            ant_z_N_cr   = stack([r.ant_z_cr["Neutral"] for r in runs_cr])  # type: ignore[index]

        # ----- FB (raw)
        fb_psc_P = stack([r.fb_psc["Positive"] for r in runs])
        fb_psc_G = stack([r.fb_psc["Negative"] for r in runs])
        fb_z_P   = stack([r.fb_z["Positive"] for r in runs])
        fb_z_G   = stack([r.fb_z["Negative"] for r in runs])

        # Optional FB (confound-regressed)
        runs_cr_fb = [r for r in runs if r.confounds_applied and r.fb_psc_cr is not None and r.fb_z_cr is not None]
        fb_psc_P_cr = fb_psc_G_cr = fb_z_P_cr = fb_z_G_cr = None
        if len(runs_cr_fb) > 0:
            fb_psc_P_cr = stack([r.fb_psc_cr["Positive"] for r in runs_cr_fb])  # type: ignore[index]
            fb_psc_G_cr = stack([r.fb_psc_cr["Negative"] for r in runs_cr_fb])  # type: ignore[index]
            fb_z_P_cr   = stack([r.fb_z_cr["Positive"] for r in runs_cr_fb])  # type: ignore[index]
            fb_z_G_cr   = stack([r.fb_z_cr["Negative"] for r in runs_cr_fb])  # type: ignore[index]

        # Use the first run's time axis (all runs for an echo should match)
        t_axis = runs[0].time_axis

        # ---- Plots
        if ant_psc_R_cr is not None and ant_psc_N_cr is not None:
            plot_two_conditions_dualconf(
                t_axis,
                ant_psc_R, "Reward",
                ant_psc_N, "Neutral",
                ant_psc_R_cr, ant_psc_N_cr,
                title=f"VS — subject {sub} [{echo}] (ANT, PSC)",
                ylabel="% signal change (PSC)",
                out_png=outdir / "anticipation_psc.png",
            )
            plot_two_conditions_dualconf(
                t_axis,
                ant_z_R, "Reward",
                ant_z_N, "Neutral",
                ant_z_R_cr, ant_z_N_cr,
                title=f"VS — subject {sub} [{echo}] (ANT, Z)",
                ylabel="Z-scored signal",
                out_png=outdir / "anticipation_z.png",
            )
        else:
            plot_two_conditions(
                t_axis,
                ant_psc_R, "Reward",
                ant_psc_N, "Neutral",
                title=f"VS — subject {sub} [{echo}] (ANT, PSC)",
                ylabel="% signal change (PSC)",
                out_png=outdir / "anticipation_psc.png",
            )
            plot_two_conditions(
                t_axis,
                ant_z_R, "Reward",
                ant_z_N, "Neutral",
                title=f"VS — subject {sub} [{echo}] (ANT, Z)",
                ylabel="Z-scored signal",
                out_png=outdir / "anticipation_z.png",
            )

        if fb_psc_P_cr is not None and fb_psc_G_cr is not None:
            plot_two_conditions_dualconf(
                t_axis,
                fb_psc_P, "Feedback +",
                fb_psc_G, "Feedback −",
                fb_psc_P_cr, fb_psc_G_cr,
                title=f"VS — subject {sub} [{echo}] (FB, PSC)",
                ylabel="% signal change (PSC)",
                out_png=outdir / "feedback_psc.png",
            )
            plot_two_conditions_dualconf(
                t_axis,
                fb_z_P, "Feedback +",
                fb_z_G, "Feedback −",
                fb_z_P_cr, fb_z_G_cr,
                title=f"VS — subject {sub} [{echo}] (FB, Z)",
                ylabel="Z-scored signal",
                out_png=outdir / "feedback_z.png",
            )
        else:
            plot_two_conditions(
                t_axis,
                fb_psc_P, "Feedback +",
                fb_psc_G, "Feedback −",
                title=f"VS — subject {sub} [{echo}] (FB, PSC)",
                ylabel="% signal change (PSC)",
                out_png=outdir / "feedback_psc.png",
            )
            plot_two_conditions(
                t_axis,
                fb_z_P, "Feedback +",
                fb_z_G, "Feedback −",
                title=f"VS — subject {sub} [{echo}] (FB, Z)",
                ylabel="Z-scored signal",
                out_png=outdir / "feedback_z.png",
            )

        # ---- Subject-level summary TSV (raw-only; same columns as before)
        summ = outdir / "summary_points.tsv"
        with summ.open("w") as f:
            f.write("sub\techo\tant_R_psc\tant_N_psc\tant_R_z\tant_N_z\tfb_pos_psc\tfb_neg_psc\tfb_pos_z\tfb_neg_z\n")
            # mean of run-level point estimates (raw)
            ant_R = float(np.nanmean([r.points_psc["ant_R"] for r in runs]))
            ant_N = float(np.nanmean([r.points_psc["ant_N"] for r in runs]))
            ant_Rz = float(np.nanmean([r.points_z["ant_R"] for r in runs]))
            ant_Nz = float(np.nanmean([r.points_z["ant_N"] for r in runs]))
            fb_P = float(np.nanmean([r.points_psc["fb_pos"] for r in runs]))
            fb_G = float(np.nanmean([r.points_psc["fb_neg"] for r in runs]))
            fb_Pz = float(np.nanmean([r.points_z["fb_pos"] for r in runs]))
            fb_Gz = float(np.nanmean([r.points_z["fb_neg"] for r in runs]))
            f.write(f"{sub}\t{echo}\t{ant_R:.6f}\t{ant_N:.6f}\t{ant_Rz:.6f}\t{ant_Nz:.6f}\t{fb_P:.6f}\t{fb_G:.6f}\t{fb_Pz:.6f}\t{fb_Gz:.6f}\n")


def main():
    if not FEAT_LIST_PATH.exists():
        print(f"[WARN] FEAT list not found: {FEAT_LIST_PATH}. Nothing to do.")
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

    # Process all FEATs
    results: List[RunResult] = []
    rows_for_tp: List[List[str]] = []

    for feat in feat_paths:
        res = process_one_feat(feat)
        if res is None:
            continue
        save_run_plots(res)
        results.append(res)
        rows_for_tp.append([
            res.sub,res.ses,res.run,res.echo,
            _fmt(res.tp_means_psc.get("ANT_REWARD")),
            _fmt(res.tp_means_psc.get("ANT_NEUTRAL")),
            _fmt(res.tp_means_psc.get("FB_POS_REWARD")),
            _fmt(res.tp_means_psc.get("FB_NEG_REWARD")),
            _fmt(res.tp_means_psc.get("FB_POS_NEUTRAL")),
            _fmt(res.tp_means_psc.get("FB_NEG_NEUTRAL")),
            _fmt(res.tp_means_z.get("ANT_REWARD")),
            _fmt(res.tp_means_z.get("ANT_NEUTRAL")),
            _fmt(res.tp_means_z.get("FB_POS_REWARD")),
            _fmt(res.tp_means_z.get("FB_NEG_REWARD")),
            _fmt(res.tp_means_z.get("FB_POS_NEUTRAL")),
            _fmt(res.tp_means_z.get("FB_NEG_NEUTRAL")),
            str(res.tp_counts.get("ANT_REWARD",0)),
            str(res.tp_counts.get("ANT_NEUTRAL",0)),
            str(res.tp_counts.get("FB_POS_REWARD",0)),
            str(res.tp_counts.get("FB_NEG_REWARD",0)),
            str(res.tp_counts.get("FB_POS_NEUTRAL",0)),
            str(res.tp_counts.get("FB_NEG_NEUTRAL",0)),
        ])

    # Aggregate by subject (split by echo)
    subjects = sorted({r.sub for r in results})
    for sub in subjects:
        aggregate_subject(results, sub)

    # Write the spreadsheet for the discrete 4th‑TR summaries
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    header = [
        "sub","ses","run","echo",
        "ANT_REWARD_PSC","ANT_NEUTRAL_PSC",
        "FB_POS_REWARD_PSC","FB_NEG_REWARD_PSC","FB_POS_NEUTRAL_PSC","FB_NEG_NEUTRAL_PSC",
        "ANT_REWARD_Z","ANT_NEUTRAL_Z",
        "FB_POS_REWARD_Z","FB_NEG_REWARD_Z","FB_POS_NEUTRAL_Z","FB_NEG_NEUTRAL_Z",
        "N_ANT_REWARD","N_ANT_NEUTRAL","N_FB_POS_REWARD","N_FB_NEG_REWARD","N_FB_POS_NEUTRAL","N_FB_NEG_NEUTRAL"
    ]
    tsv_path = SUMMARY_DIR / "summary_at_4thTR_mid-unsmoothed.tsv"

    with open(tsv_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in rows_for_tp:
            f.write("\t".join(row) + "\n")


    print(f"Done. Outputs in: {OUT_TC_DIR}\n  - Discrete summaries: {tsv_path}")


if __name__ == "__main__":
    main()
