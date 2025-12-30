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
from typing import Dict, List, Tuple, Optional

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# --------------------------- Fixed project structure ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
MASKS_DIR  = ROOT_DIR / "masks"
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-unsmoothed-new"
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
    # Optional nuisance-regressed curves (same format as above).
    # When confounds_used==0 these dicts are empty and plotting falls back to raw curves only.
    confounds_used: int = 0
    n_confounds: int = 0
    ant_psc_conf: Dict[str, Tuple[np.ndarray, np.ndarray, int]] = None  # type: ignore[assignment]
    ant_z_conf:   Dict[str, Tuple[np.ndarray, np.ndarray, int]] = None  # type: ignore[assignment]
    fb_psc_conf:  Dict[str, Tuple[np.ndarray, np.ndarray, int]] = None  # type: ignore[assignment]
    fb_z_conf:    Dict[str, Tuple[np.ndarray, np.ndarray, int]] = None  # type: ignore[assignment]

    # Timing markers (seconds; relative to the event used for extraction/plotting).
    cue_dur_s: float = 0.75  # cue display duration (used to mark cue-onset at -cue_dur_s)
    isi_median_s: float = float('nan')   # median target-onset relative to event (from ANT EV durations)
    fb_delay_median_s: float = float('nan')  # median feedback-onset relative to event (paired as next FB after ANT onset)

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


# --------------------------- Optional confound regression ----------------------

def find_confoundevs_file(feat: Path) -> Optional[Path]:
    """Return confoundevs.txt if present inside a FEAT directory."""
    cand = feat / "confoundevs.txt"
    return cand if cand.exists() else None


def load_confound_matrix(conf_path: Path, T: int) -> Optional[np.ndarray]:
    """Load confound matrix (T x K). Returns None if missing/mismatched."""
    try:
        X = np.loadtxt(conf_path)
    except Exception as e:
        print(f"[WARN] Could not read confounds: {conf_path} ({e})")
        return None

    if X.ndim == 1:
        X = X.reshape(-1, 1)

    if X.shape[0] != T:
        print(f"[WARN] Confound rows != timeseries length for {conf_path}: {X.shape[0]} vs {T}. Skipping conf-reg.")
        return None

    # Drop constant columns (std==0) to avoid rank issues.
    keep = np.nanstd(X, axis=0) > 0
    if keep.ndim == 0:
        keep = np.array([bool(keep)])
    X = X[:, keep]
    if X.shape[1] == 0:
        return None
    return X


def regress_out(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Regress X (plus intercept) out of y and return cleaned y preserving mean."""
    if X is None or X.size == 0:
        return y

    # Add intercept
    X_ = np.column_stack([X, np.ones(X.shape[0])])
    try:
        beta, *_ = np.linalg.lstsq(X_, y, rcond=None)
        y_hat = X_ @ beta
        resid = y - y_hat
        return resid + np.nanmean(y)
    except Exception as e:
        print(f"[WARN] Confound regression failed: {e}. Returning raw y.")
        return y


# --------------------------- EV timing helpers (markers) -----------------------

def load_ev_3col(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load an EV file with onset, duration, amplitude."""
    if not path.exists():
        return np.array([]), np.array([]), np.array([])
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 3:
        raise ValueError(f"EV file does not have 3 columns: {path}")
    return arr[:, 0], arr[:, 1], arr[:, 2]


def median_isi_from_ant_evs(antR_path: Path, antN_path: Path) -> float:
    """Median target-onset relative to ANT event, based on EV durations."""
    _, dR, _ = load_ev_3col(antR_path)
    _, dN, _ = load_ev_3col(antN_path)
    d = np.concatenate([dR, dN]) if (dR.size or dN.size) else np.array([])
    return float(np.nanmedian(d)) if d.size else float("nan")


def median_fb_delay_from_evs(ant_onsets: np.ndarray, fb_onsets: np.ndarray) -> float:
    """Median feedback-onset delay relative to ANT event (paired as next FB after ANT onset)."""
    if ant_onsets.size == 0 or fb_onsets.size == 0:
        return float("nan")
    fb_onsets = np.sort(fb_onsets)
    ant_onsets = np.sort(ant_onsets)

    delays = []
    for t in ant_onsets:
        j = np.searchsorted(fb_onsets, t, side="right")
        if j < fb_onsets.size:
            dt = fb_onsets[j] - t
            if dt >= 0:
                delays.append(dt)
    if not delays:
        return float("nan")
    return float(np.nanmedian(np.array(delays)))


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


def plot_two_conditions(
    time_axis: np.ndarray,
    condA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
    condB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
    title: str,
    ylabel: str,
    out_png: Path,
    vlines: Optional[List[Tuple[float, str]]] = None,
    overlayA: Optional[Tuple[np.ndarray, np.ndarray, int]] = None,
    overlayB: Optional[Tuple[np.ndarray, np.ndarray, int]] = None,
    overlay_label: str = "conf-reg",
) -> None:
    """Plot two conditions, optionally overlaying nuisance-regressed curves as dashed lines."""
    mA, sA, nA = condA
    mB, sB, nB = condB

    fig, ax = plt.subplots(figsize=(8, 4.8))

    # Handle empty inputs gracefully (write a diagnostic plot instead of crashing).
    if mA is None or mB is None or len(mA) == 0 or len(mB) == 0:
        ax.text(0.5, 0.5, "No data to plot (empty condition curves).", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        ax.set_xlabel("Time from event (s)")
        ax.set_ylabel(ylabel)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)
        return

    ax.plot(time_axis, mA, label=f"{labelA} (n={nA})")
    ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)
    ax.plot(time_axis, mB, label=f"{labelB} (n={nB})")
    ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)

    # Optional overlays (dashed, no shading to keep legible)
    if overlayA is not None and overlayB is not None:
        oA, _, onA = overlayA
        oB, _, onB = overlayB
        if oA is not None and oB is not None and len(oA) == len(time_axis) and len(oB) == len(time_axis):
            ax.plot(time_axis, oA, ls="--", lw=1.5, label=f"{labelA} ({overlay_label}; n={onA})")
            ax.plot(time_axis, oB, ls="--", lw=1.5, label=f"{labelB} ({overlay_label}; n={onB})")

    # Vertical markers
    if vlines is None:
        # fall back to the original defaults
        for v in VERT_LINES:
            ax.axvline(v, ls=":", lw=1)
    else:
        for x, lab in vlines:
            ax.axvline(x, ls=":", lw=1)
            if lab:
                # put label near the top of the axis
                y_top = ax.get_ylim()[1]
                ax.text(x, y_top, f" {lab}", rotation=90, va="bottom", ha="left", fontsize=8)

    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
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

    # Optional nuisance regression using FEAT confound EVs (if available).
    confounds_used = 0
    n_confounds = 0
    ts_psc_conf = None
    ts_z_conf = None

    conf_path = find_confoundevs_file(feat)
    if conf_path is not None:
        X = load_confound_matrix(conf_path, T)
        if X is not None:
            n_confounds = int(X.shape[1])
            ts_clean = regress_out(ts_raw, X)
            ts_psc_conf = to_psc(ts_clean)
            ts_z_conf   = to_z(ts_clean)
            confounds_used = 1

    T = ts_raw.size

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)
    ant_R = load_ev(ev_dir / "_anticipation_reward.txt")
    ant_N = load_ev(ev_dir / "_anticipation_neutral.txt")
    # For timing markers, we also read the 3-column ANT EVs (durations = cue-offset -> target-onset).
    ant_R_path = ev_dir / "_anticipation_reward.txt"
    ant_N_path = ev_dir / "_anticipation_neutral.txt"

    # Feedback (pooled by valence for plotting)
    fb_pos = np.sort(np.concatenate([
        load_ev(ev_dir / "_feedback_positive_reward.txt"),
        load_ev(ev_dir / "_feedback_positive_neutral.txt"),
    ]))
    fb_neg = np.sort(np.concatenate([
        load_ev(ev_dir / "_feedback_negative_reward.txt"),
        load_ev(ev_dir / "_feedback_negative_neutral.txt"),
    ]))


    # Timing markers (all relative to the ANT event onset used in extraction).
    all_ant_onsets = np.sort(np.concatenate([ant_R, ant_N])) if (ant_R.size or ant_N.size) else np.array([])
    all_fb_onsets = np.sort(np.concatenate([fb_pos, fb_neg])) if (fb_pos.size or fb_neg.size) else np.array([])

    isi_median_s = median_isi_from_ant_evs(ant_R_path, ant_N_path)
    fb_delay_median_s = median_fb_delay_from_evs(all_ant_onsets, all_fb_onsets)

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


    # Confound-regressed curves (optional). We compute the same windows if ts_psc_conf/ts_z_conf exist.
    ant_psc_conf = {}
    ant_z_conf   = {}
    fb_psc_conf  = {}
    fb_z_conf    = {}
    if confounds_used and ts_psc_conf is not None and ts_z_conf is not None:
        R_psc_w_c = sample_windows(ts_psc_conf, ant_R, tr, TMIN, TMAX)
        N_psc_w_c = sample_windows(ts_psc_conf, ant_N, tr, TMIN, TMAX)
        R_z_w_c   = sample_windows(ts_z_conf,   ant_R, tr, TMIN, TMAX)
        N_z_w_c   = sample_windows(ts_z_conf,   ant_N, tr, TMIN, TMAX)

        ant_psc_conf = {
            "Reward": mean_and_sem(R_psc_w_c),
            "Neutral": mean_and_sem(N_psc_w_c),
        }
        ant_z_conf = {
            "Reward": mean_and_sem(R_z_w_c),
            "Neutral": mean_and_sem(N_z_w_c),
        }

        P_psc_w_c = sample_windows(ts_psc_conf, fb_pos, tr, TMIN, TMAX)
        G_psc_w_c = sample_windows(ts_psc_conf, fb_neg, tr, TMIN, TMAX)
        P_z_w_c   = sample_windows(ts_z_conf,   fb_pos, tr, TMIN, TMAX)
        G_z_w_c   = sample_windows(ts_z_conf,   fb_neg, tr, TMIN, TMAX)

        fb_psc_conf = {
            "Positive": mean_and_sem(P_psc_w_c),
            "Negative": mean_and_sem(G_psc_w_c),
        }
        fb_z_conf = {
            "Positive": mean_and_sem(P_z_w_c),
            "Negative": mean_and_sem(G_z_w_c),
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
        confounds_used=confounds_used,
        n_confounds=n_confounds,
        ant_psc_conf=ant_psc_conf,
        ant_z_conf=ant_z_conf,
        fb_psc_conf=fb_psc_conf,
        fb_z_conf=fb_z_conf,
        cue_dur_s=0.75,
        isi_median_s=isi_median_s,
        fb_delay_median_s=fb_delay_median_s,
    )


def save_run_plots(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

    # Build plot-specific timing markers.
    # ANT plots are anchored at the ANT event onset used in extraction (t=0).
    vlines_ant: List[Tuple[float, str]] = [
        (-res.cue_dur_s, "cue onset"),
        (0.0, "cue offset / ANT onset"),
    ]
    if not np.isnan(res.isi_median_s):
        vlines_ant.append((res.isi_median_s, "median target onset"))
    if not np.isnan(res.fb_delay_median_s):
        vlines_ant.append((res.fb_delay_median_s, "median feedback onset"))
    # Wu-style comparability: 6 s after cue onset.
    vlines_ant.append((6.0 - res.cue_dur_s, "Wu: +6s from cue onset"))
    vlines_ant.append((6.0, "+6s from ANT onset"))

    # FB plots are anchored at feedback onset (t=0).
    vlines_fb: List[Tuple[float, str]] = [
        (0.0, "feedback onset"),
        (6.0, "+6s"),
    ]

    # Optional overlay curves (confound-regressed) — plotted as dashed lines.
    ant_psc_oR = res.ant_psc_conf.get("Reward") if (res.confounds_used and res.ant_psc_conf) else None
    ant_psc_oN = res.ant_psc_conf.get("Neutral") if (res.confounds_used and res.ant_psc_conf) else None
    ant_z_oR   = res.ant_z_conf.get("Reward")   if (res.confounds_used and res.ant_z_conf)   else None
    ant_z_oN   = res.ant_z_conf.get("Neutral")  if (res.confounds_used and res.ant_z_conf)   else None

    fb_psc_oP  = res.fb_psc_conf.get("Positive") if (res.confounds_used and res.fb_psc_conf) else None
    fb_psc_oN  = res.fb_psc_conf.get("Negative") if (res.confounds_used and res.fb_psc_conf) else None
    fb_z_oP    = res.fb_z_conf.get("Positive")   if (res.confounds_used and res.fb_z_conf)   else None
    fb_z_oN    = res.fb_z_conf.get("Negative")   if (res.confounds_used and res.fb_z_conf)   else None



    # Anticipation: PSC and Z
    plot_two_conditions(
        res.time_axis,
        res.ant_psc["Reward"], "Reward",
        res.ant_psc["Neutral"], "Neutral",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, PSC)",
        ylabel="% signal change (PSC)",
        out_png=run_out / "anticipation_psc.png",
        vlines=vlines_ant,
        overlayA=ant_psc_oR,
        overlayB=ant_psc_oN,
        overlay_label="conf-reg",
    )
    plot_two_conditions(
        res.time_axis,
        res.ant_z["Reward"], "Reward",
        res.ant_z["Neutral"], "Neutral",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, Z)",
        ylabel="Z (SD units)",
        out_png=run_out / "anticipation_z.png",
        vlines=vlines_ant,
        overlayA=ant_z_oR,
        overlayB=ant_z_oN,
        overlay_label="conf-reg",
    )

    # Feedback (valence): PSC and Z
    plot_two_conditions(
        res.time_axis,
        res.fb_psc["Positive"], "Feedback +",
        res.fb_psc["Negative"], "Feedback −",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, PSC)",
        ylabel="% signal change (PSC)",
        out_png=run_out / "feedback_psc.png",
        vlines=vlines_fb,
        overlayA=fb_psc_oP,
        overlayB=fb_psc_oN,
        overlay_label="conf-reg",
    )
    plot_two_conditions(
        res.time_axis,
        res.fb_z["Positive"], "Feedback +",
        res.fb_z["Negative"], "Feedback −",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, Z)",
        ylabel="Z (SD units)",
        out_png=run_out / "feedback_z.png",
    )

    # Save a summary file
    summary = (
        f"Run: sub-{res.sub} ses-{res.ses} run-{res.run} echo={res.echo}\n"
        f"ANT  Reward n={res.ant_psc['Reward'][2]}, Neutral n={res.ant_psc['Neutral'][2]}\n"
        f"FB   Positive n={res.fb_psc['Positive'][2]}, Negative n={res.fb_psc['Negative'][2]}\n"
    )
    txt_path = run_out / "summary.txt"
    txt_path.parent.mkdir(parents=True, exist_ok=True)
    txt_path.write_text(summary)

    # Save per‑run discrete 4th‑TR summaries (TSV)
    tp_path = run_out / "tp_4thTR.tsv"
    header = [
        "sub","ses","run","echo",
        "ANT_REWARD_PSC","ANT_NEUTRAL_PSC",
        "FB_POS_REWARD_PSC","FB_NEG_REWARD_PSC","FB_POS_NEUTRAL_PSC","FB_NEG_NEUTRAL_PSC",
        "ANT_REWARD_Z","ANT_NEUTRAL_Z",
        "FB_POS_REWARD_Z","FB_NEG_REWARD_Z","FB_POS_NEUTRAL_Z","FB_NEG_NEUTRAL_Z",
        "N_ANT_REWARD","N_ANT_NEUTRAL","N_FB_POS_REWARD","N_FB_NEG_REWARD","N_FB_POS_NEUTRAL","N_FB_NEG_NEUTRAL"
    ]
    row = [
        res.sub, res.ses, res.run, res.echo,
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
    ]
    tp_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tp_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")


# -------------------------- Subject‑level aggregation --------------------------

def aggregate_subject(results: List[RunResult], subject: str) -> None:
    pool: Dict[str, Dict[str, Dict[str, List[Tuple[np.ndarray, int]]]]] = {}
    time_axis = None

    # Collect timing markers across runs (for ANT plots).
    marker_pool: Dict[str, Dict[str, List[float]]] = {}  # echo -> {'isi':[], 'fb':[], 'cue':[]}


    for res in results:
        if res.sub != subject:
            continue
        if time_axis is None:
            time_axis = res.time_axis
        echo = res.echo
        marker_pool.setdefault(echo, {'cue': [], 'isi': [], 'fb': []})
        marker_pool[echo]['cue'].append(res.cue_dur_s)
        if not np.isnan(res.isi_median_s):
            marker_pool[echo]['isi'].append(res.isi_median_s)
        if not np.isnan(res.fb_delay_median_s):
            marker_pool[echo]['fb'].append(res.fb_delay_median_s)

        pool.setdefault(echo, {
            "ANT_PSC": {"Reward": [], "Neutral": []},
            "ANT_Z":   {"Reward": [], "Neutral": []},
            "FB_PSC":  {"Positive": [], "Negative": []},
            "FB_Z":    {"Positive": [], "Negative": []},
            "ANT_PSC_CONF": {"Reward": [], "Neutral": []},
            "ANT_Z_CONF":   {"Reward": [], "Neutral": []},
            "FB_PSC_CONF":  {"Positive": [], "Negative": []},
            "FB_Z_CONF":    {"Positive": [], "Negative": []},
        })
        pool[echo]["ANT_PSC"]["Reward"].append((res.ant_psc["Reward"][0], res.ant_psc["Reward"][2]))
        pool[echo]["ANT_PSC"]["Neutral"].append((res.ant_psc["Neutral"][0], res.ant_psc["Neutral"][2]))
        pool[echo]["ANT_Z"]["Reward"].append((res.ant_z["Reward"][0], res.ant_z["Reward"][2]))
        pool[echo]["ANT_Z"]["Neutral"].append((res.ant_z["Neutral"][0], res.ant_z["Neutral"][2]))
        pool[echo]["FB_PSC"]["Positive"].append((res.fb_psc["Positive"][0], res.fb_psc["Positive"][2]))
        pool[echo]["FB_PSC"]["Negative"].append((res.fb_psc["Negative"][0], res.fb_psc["Negative"][2]))
        pool[echo]["FB_Z"]["Positive"].append((res.fb_z["Positive"][0], res.fb_z["Positive"][2]))
        pool[echo]["FB_Z"]["Negative"].append((res.fb_z["Negative"][0], res.fb_z["Negative"][2]))
        # Optional confound-regressed curves
        if res.confounds_used and res.ant_psc_conf and res.fb_psc_conf:
            pool[echo]["ANT_PSC_CONF"]["Reward"].append((res.ant_psc_conf["Reward"][0], res.ant_psc_conf["Reward"][2]))
            pool[echo]["ANT_PSC_CONF"]["Neutral"].append((res.ant_psc_conf["Neutral"][0], res.ant_psc_conf["Neutral"][2]))
            pool[echo]["ANT_Z_CONF"]["Reward"].append((res.ant_z_conf["Reward"][0], res.ant_z_conf["Reward"][2]))
            pool[echo]["ANT_Z_CONF"]["Neutral"].append((res.ant_z_conf["Neutral"][0], res.ant_z_conf["Neutral"][2]))
            pool[echo]["FB_PSC_CONF"]["Positive"].append((res.fb_psc_conf["Positive"][0], res.fb_psc_conf["Positive"][2]))
            pool[echo]["FB_PSC_CONF"]["Negative"].append((res.fb_psc_conf["Negative"][0], res.fb_psc_conf["Negative"][2]))
            pool[echo]["FB_Z_CONF"]["Positive"].append((res.fb_z_conf["Positive"][0], res.fb_z_conf["Positive"][2]))
            pool[echo]["FB_Z_CONF"]["Negative"].append((res.fb_z_conf["Negative"][0], res.fb_z_conf["Negative"][2]))


    if time_axis is None:
        print(f"[INFO] No runs aggregated for subject {subject}.")
        return

    subj_out = OUT_TC_DIR / "subjects" / f"sub-{subject}"

    def weighted_mean_and_sem(items: List[Tuple[np.ndarray, int]]) -> Tuple[np.ndarray, np.ndarray, int]:
        if len(items) == 0:
            return np.array([]), np.array([]), 0
        means = np.vstack([m for (m, n) in items])
        sem = np.nanstd(means, axis=0, ddof=1) / np.sqrt(max(len(items), 1))
        n_total = int(np.sum([max(n, 0) for (m, n) in items]))
        mean = np.nanmean(means, axis=0)
        return mean, sem, n_total

    for echo, families in pool.items():
        # Subject-level timing markers (ANT plots).
        cue_dur = float(np.nanmedian(np.array(marker_pool.get(echo, {}).get('cue', [0.75]))))
        isi_med = float(np.nanmedian(np.array(marker_pool.get(echo, {}).get('isi', [])))) if marker_pool.get(echo, {}).get('isi') else float('nan')
        fb_med  = float(np.nanmedian(np.array(marker_pool.get(echo, {}).get('fb', []))))  if marker_pool.get(echo, {}).get('fb')  else float('nan')

        vlines_ant: List[Tuple[float, str]] = [(-cue_dur, 'cue onset'), (0.0, 'cue offset / ANT onset')]
        if not np.isnan(isi_med):
            vlines_ant.append((isi_med, 'median target onset'))
        if not np.isnan(fb_med):
            vlines_ant.append((fb_med, 'median feedback onset'))
        vlines_ant.append((6.0 - cue_dur, 'Wu: +6s from cue onset'))
        vlines_ant.append((6.0, '+6s from ANT onset'))

        vlines_fb: List[Tuple[float, str]] = [(0.0, 'feedback onset'), (6.0, '+6s')]

        # ANTICIPATION — PSC
        ant_psc_R = weighted_mean_and_sem(families["ANT_PSC"]["Reward"])
        ant_psc_N = weighted_mean_and_sem(families["ANT_PSC"]["Neutral"])
        # Optional confound-regressed overlays (computed across runs that had confounds).
        ant_psc_Rc = weighted_mean_and_sem(families["ANT_PSC_CONF"]["Reward"])
        ant_psc_Nc = weighted_mean_and_sem(families["ANT_PSC_CONF"]["Neutral"])
        plot_two_conditions(
            time_axis, ant_psc_R, "Reward", ant_psc_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "anticipation_psc.png",
            vlines=vlines_ant,
            overlayA=ant_psc_Rc,
            overlayB=ant_psc_Nc,
            overlay_label="conf-reg",
        )
        # ANTICIPATION — Z
        ant_z_R = weighted_mean_and_sem(families["ANT_Z"]["Reward"])
        ant_z_N = weighted_mean_and_sem(families["ANT_Z"]["Neutral"])
        ant_z_Rc = weighted_mean_and_sem(families["ANT_Z_CONF"]["Reward"])
        ant_z_Nc = weighted_mean_and_sem(families["ANT_Z_CONF"]["Neutral"])
        plot_two_conditions(
            time_axis, ant_z_R, "Reward", ant_z_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "anticipation_z.png",
        )
        # FEEDBACK (valence) — PSC
        fb_psc_P = weighted_mean_and_sem(families["FB_PSC"]["Positive"])
        fb_psc_N = weighted_mean_and_sem(families["FB_PSC"]["Negative"])
        fb_psc_Pc = weighted_mean_and_sem(families["FB_PSC_CONF"]["Positive"])
        fb_psc_Nc = weighted_mean_and_sem(families["FB_PSC_CONF"]["Negative"])
        plot_two_conditions(
            time_axis, fb_psc_P, "Feedback +", fb_psc_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "feedback_psc.png",
        )
        # FEEDBACK (valence) — Z
        fb_z_P = weighted_mean_and_sem(families["FB_Z"]["Positive"])
        fb_z_N = weighted_mean_and_sem(families["FB_Z"]["Negative"])
        fb_z_Pc = weighted_mean_and_sem(families["FB_Z_CONF"]["Positive"])
        fb_z_Nc = weighted_mean_and_sem(families["FB_Z_CONF"]["Negative"])
        plot_two_conditions(
            time_axis, fb_z_P, "Feedback +", fb_z_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "feedback_z.png",
        )
        summary = (
            f"Subject {subject} — Echo: {echo}\n"
            f"ANT  Reward n={ant_psc_R[2]}, Neutral n={ant_psc_N[2]}\n"
            f"FB   Positive n={fb_psc_P[2]}, Negative n={fb_psc_N[2]}\n"
        )
        (subj_out / echo / "summary.txt").write_text(summary)


# ---------------------------------- Driver ------------------------------------

def _fmt(x: float | None) -> str:
    if x is None or not np.isfinite(x):
        return ""
    return f"{x:.6f}"

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
    tr_values: List[float] = []

    rows_for_tp: List[List[str]] = []

    for feat in feat_paths:
        res = process_one_feat(feat)
        if res is None:
            continue
        save_run_plots(res)
        results.append(res)
        tr_values.append(float(res.tr_used))
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

    # TR sanity check across all processed runs
    if tr_values:
        tr_unique = sorted(set([round(v, 6) for v in tr_values]))
        if len(tr_unique) > 1:
            print(f"[WARN] Multiple TR values observed across runs: {tr_unique}")
        else:
            print(f"[INFO] TR consistent across runs: {tr_unique[0]:.6f} s")

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
    tsv_path = SUMMARY_DIR / "summary_at_4thTR_mid-unsmoothed-new.tsv"

    with open(tsv_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in rows_for_tp:
            f.write("\t".join(row) + "\n")


    print(f"Done. Outputs in: {OUT_TC_DIR}\n  - Discrete summaries: {tsv_path}")


if __name__ == "__main__":
    main()
