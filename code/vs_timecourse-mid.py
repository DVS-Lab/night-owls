#!/usr/bin/env python3
"""
VS time‑course extraction, plotting, and discrete 4th‑TR summaries (plus Wu‑style ANT@+6s)
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
    • ANTICIPATION: Reward vs Neutral (delay‑onset locked; cue‑offset)
    • FEEDBACK (valence): Positive vs Negative (feedback‑locked, pooling across
      incentive conditions)
- Produces per‑run plots and subject‑level aggregated plots split by echo
  (single‑echo vs multi‑echo) for PSC and Z (four figures per subject per family).
- Additionally, writes **discrete 4th‑TR‑after‑onset summaries** (no interpolation), plus an interpolated Wu‑style +6.0 s point estimate for ANT
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
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-unsmoothed_interpolated"
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
# Timing constants
# - In your MID, cue duration is fixed at 0.75 s; the ISI is jittered 1.5–3.0 s.
# - To align "anticipation" with Wu et al. (delay onset), we treat ANT t=0 as cue offset.
CUE_DUR_S   = 0.750  # seconds
WU_LAG_S    = 6.0    # seconds after ANT delay onset (Wu-style NAcc peak lag)

# Vertical reference lines (sec from event)
# ANT plots are delay-onset locked: 0 = cue offset; target occurs 1.5–3.0 s later.
VERT_LINES_ANT      = [0.0, 1.5, 3.0, 6.0]
# FB plots are feedback-onset locked: 0 = feedback onset; 6 s is a helpful hemodynamic reference.
VERT_LINES_FEEDBACK = [0.0, 6.0]
# Default if not otherwise specified
VERT_LINES_DEFAULT  = [0.0, 6.0]

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
    # Discrete summaries: 4th-TR (no interpolation) + ANT@+6s (interpolated)
    tp_means_psc: Dict[str, float | None]
    tp_means_z:   Dict[str, float | None]
    tp_counts:    Dict[str, int]

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


def vals_at_offset(ts: np.ndarray, onsets_s: np.ndarray, tr: float, offset_s: float) -> np.ndarray:
    """Sample ts at (onset + offset_s) using linear interpolation.

    Uses left/right NaN padding so events too close to the run edges won't contribute
    to the mean (handled via np.nanmean).
    """
    if onsets_s.size == 0:
        return np.empty((0,), dtype=float)
    t_series = np.arange(ts.size) * tr
    t_abs = onsets_s + float(offset_s)
    return np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)



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
                        title: str, ylabel: str, out_png: Path,
                        vlines: List[float] | None = None) -> None:
    mA, sA, nA = condA
    mB, sB, nB = condB
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(time_axis, mA, label=f"{labelA} (n={nA})")
    ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)
    ax.plot(time_axis, mB, label=f"{labelB} (n={nB})")
    ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)
    for v in (VERT_LINES_DEFAULT if vlines is None else vlines):
        ax.axvline(v, ls=":", lw=1)
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
    T = ts_raw.size

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)
    ant_R = load_ev(ev_dir / "_anticipation_reward.txt")
    ant_N = load_ev(ev_dir / "_anticipation_neutral.txt")

    # Wu-style alignment for anticipation (NAcc): lock to delay onset (cue offset)
    # Your EV onsets are cue-onset; shift by cue duration so t=0 is cue offset.
    ant_R = ant_R + CUE_DUR_S
    ant_N = ant_N + CUE_DUR_S

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
        if arr.size == 0:
            return None
        if not np.any(~np.isnan(arr)):
            return None
        return float(np.nanmean(arr))

    def vals_at_tp(ts: np.ndarray, onsets: np.ndarray) -> np.ndarray:
        idx = fourth_tr_indices(onsets, tr, T)
        return ts[idx]

    # Anticipation (Reward/Neutral)
    ant_R_psc_tp = vals_at_tp(ts_psc, ant_R)
    ant_N_psc_tp = vals_at_tp(ts_psc, ant_N)
    ant_R_z_tp   = vals_at_tp(ts_z,   ant_R)
    ant_N_z_tp   = vals_at_tp(ts_z,   ant_N)

    # Wu-style point estimate for NAcc/VS: sample exactly +6.0 s after delay onset
    # using linear interpolation (avoids TR rounding).
    ant_R_psc_wu6 = vals_at_offset(ts_psc, ant_R, tr, WU_LAG_S)
    ant_N_psc_wu6 = vals_at_offset(ts_psc, ant_N, tr, WU_LAG_S)
    ant_R_z_wu6   = vals_at_offset(ts_z,   ant_R, tr, WU_LAG_S)
    ant_N_z_wu6   = vals_at_offset(ts_z,   ant_N, tr, WU_LAG_S)

    def n_valid(arr: np.ndarray) -> int:
        return int(np.sum(~np.isnan(arr))) if arr.size else 0

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
        "ANT_REWARD_WU6":  m(ant_R_psc_wu6),
        "ANT_NEUTRAL_WU6": m(ant_N_psc_wu6),
        "FB_POS_REWARD":  m(fb_PR_psc_tp),
        "FB_NEG_REWARD":  m(fb_NR_psc_tp),
        "FB_POS_NEUTRAL": m(fb_PN_psc_tp),
        "FB_NEG_NEUTRAL": m(fb_NN_psc_tp),
    }
    tp_means_z = {
        "ANT_REWARD":  m(ant_R_z_tp),
        "ANT_NEUTRAL": m(ant_N_z_tp),
        "ANT_REWARD_WU6":  m(ant_R_z_wu6),
        "ANT_NEUTRAL_WU6": m(ant_N_z_wu6),
        "FB_POS_REWARD":  m(fb_PR_z_tp),
        "FB_NEG_REWARD":  m(fb_NR_z_tp),
        "FB_POS_NEUTRAL": m(fb_PN_z_tp),
        "FB_NEG_NEUTRAL": m(fb_NN_z_tp),
    }
    tp_counts = {
        "ANT_REWARD":  int(ant_R_psc_tp.size),
        "ANT_NEUTRAL": int(ant_N_psc_tp.size),
        "ANT_REWARD_WU6":  n_valid(ant_R_psc_wu6),
        "ANT_NEUTRAL_WU6": n_valid(ant_N_psc_wu6),
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
    )


def save_run_plots(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

    # Anticipation: PSC and Z
    plot_two_conditions(
        res.time_axis,
        res.ant_psc["Reward"], "Reward",
        res.ant_psc["Neutral"], "Neutral",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, PSC)",
        ylabel="% signal change (PSC)",
        out_png=run_out / "anticipation_psc.png",
        vlines=VERT_LINES_ANT,
    )
    plot_two_conditions(
        res.time_axis,
        res.ant_z["Reward"], "Reward",
        res.ant_z["Neutral"], "Neutral",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, Z)",
        ylabel="Z (SD units)",
        out_png=run_out / "anticipation_z.png",
        vlines=VERT_LINES_ANT,
    )

    # Feedback (valence): PSC and Z
    plot_two_conditions(
        res.time_axis,
        res.fb_psc["Positive"], "Feedback +",
        res.fb_psc["Negative"], "Feedback −",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, PSC)",
        ylabel="% signal change (PSC)",
        out_png=run_out / "feedback_psc.png",
        vlines=VERT_LINES_FEEDBACK,
    )
    plot_two_conditions(
        res.time_axis,
        res.fb_z["Positive"], "Feedback +",
        res.fb_z["Negative"], "Feedback −",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, Z)",
        ylabel="Z (SD units)",
        out_png=run_out / "feedback_z.png",
        vlines=VERT_LINES_FEEDBACK,
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
        "ANT_REWARD_PSC","ANT_NEUTRAL_PSC","ANT_REWARD_PSC_WU6","ANT_NEUTRAL_PSC_WU6",
        "FB_POS_REWARD_PSC","FB_NEG_REWARD_PSC","FB_POS_NEUTRAL_PSC","FB_NEG_NEUTRAL_PSC",
        "ANT_REWARD_Z","ANT_NEUTRAL_Z","ANT_REWARD_Z_WU6","ANT_NEUTRAL_Z_WU6",
        "FB_POS_REWARD_Z","FB_NEG_REWARD_Z","FB_POS_NEUTRAL_Z","FB_NEG_NEUTRAL_Z",
        "N_ANT_REWARD","N_ANT_NEUTRAL","N_ANT_REWARD_WU6","N_ANT_NEUTRAL_WU6","N_FB_POS_REWARD","N_FB_NEG_REWARD","N_FB_POS_NEUTRAL","N_FB_NEG_NEUTRAL"
    ]
    row = [
        res.sub, res.ses, res.run, res.echo,
        _fmt(res.tp_means_psc.get("ANT_REWARD")),
        _fmt(res.tp_means_psc.get("ANT_NEUTRAL")),
        _fmt(res.tp_means_psc.get("ANT_REWARD_WU6")),
        _fmt(res.tp_means_psc.get("ANT_NEUTRAL_WU6")),
        _fmt(res.tp_means_psc.get("FB_POS_REWARD")),
        _fmt(res.tp_means_psc.get("FB_NEG_REWARD")),
        _fmt(res.tp_means_psc.get("FB_POS_NEUTRAL")),
        _fmt(res.tp_means_psc.get("FB_NEG_NEUTRAL")),
        _fmt(res.tp_means_z.get("ANT_REWARD")),
        _fmt(res.tp_means_z.get("ANT_NEUTRAL")),
        _fmt(res.tp_means_z.get("ANT_REWARD_WU6")),
        _fmt(res.tp_means_z.get("ANT_NEUTRAL_WU6")),
        _fmt(res.tp_means_z.get("FB_POS_REWARD")),
        _fmt(res.tp_means_z.get("FB_NEG_REWARD")),
        _fmt(res.tp_means_z.get("FB_POS_NEUTRAL")),
        _fmt(res.tp_means_z.get("FB_NEG_NEUTRAL")),
        str(res.tp_counts.get("ANT_REWARD",0)),
        str(res.tp_counts.get("ANT_NEUTRAL",0)),
        str(res.tp_counts.get("ANT_REWARD_WU6",0)),
        str(res.tp_counts.get("ANT_NEUTRAL_WU6",0)),
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

    for res in results:
        if res.sub != subject:
            continue
        if time_axis is None:
            time_axis = res.time_axis
        echo = res.echo
        pool.setdefault(echo, {
            "ANT_PSC": {"Reward": [], "Neutral": []},
            "ANT_Z":   {"Reward": [], "Neutral": []},
            "FB_PSC":  {"Positive": [], "Negative": []},
            "FB_Z":    {"Positive": [], "Negative": []},
        })
        pool[echo]["ANT_PSC"]["Reward"].append((res.ant_psc["Reward"][0], res.ant_psc["Reward"][2]))
        pool[echo]["ANT_PSC"]["Neutral"].append((res.ant_psc["Neutral"][0], res.ant_psc["Neutral"][2]))
        pool[echo]["ANT_Z"]["Reward"].append((res.ant_z["Reward"][0], res.ant_z["Reward"][2]))
        pool[echo]["ANT_Z"]["Neutral"].append((res.ant_z["Neutral"][0], res.ant_z["Neutral"][2]))
        pool[echo]["FB_PSC"]["Positive"].append((res.fb_psc["Positive"][0], res.fb_psc["Positive"][2]))
        pool[echo]["FB_PSC"]["Negative"].append((res.fb_psc["Negative"][0], res.fb_psc["Negative"][2]))
        pool[echo]["FB_Z"]["Positive"].append((res.fb_z["Positive"][0], res.fb_z["Positive"][2]))
        pool[echo]["FB_Z"]["Negative"].append((res.fb_z["Negative"][0], res.fb_z["Negative"][2]))

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
        # ANTICIPATION — PSC
        ant_psc_R = weighted_mean_and_sem(families["ANT_PSC"]["Reward"])
        ant_psc_N = weighted_mean_and_sem(families["ANT_PSC"]["Neutral"])
        plot_two_conditions(
            time_axis, ant_psc_R, "Reward", ant_psc_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "anticipation_psc.png",
            vlines=VERT_LINES_ANT,
        )
        # ANTICIPATION — Z
        ant_z_R = weighted_mean_and_sem(families["ANT_Z"]["Reward"])
        ant_z_N = weighted_mean_and_sem(families["ANT_Z"]["Neutral"])
        plot_two_conditions(
            time_axis, ant_z_R, "Reward", ant_z_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "anticipation_z.png",
            vlines=VERT_LINES_ANT,
        )
        # FEEDBACK (valence) — PSC
        fb_psc_P = weighted_mean_and_sem(families["FB_PSC"]["Positive"])
        fb_psc_N = weighted_mean_and_sem(families["FB_PSC"]["Negative"])
        plot_two_conditions(
            time_axis, fb_psc_P, "Feedback +", fb_psc_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "feedback_psc.png",
            vlines=VERT_LINES_FEEDBACK,
        )
        # FEEDBACK (valence) — Z
        fb_z_P = weighted_mean_and_sem(families["FB_Z"]["Positive"])
        fb_z_N = weighted_mean_and_sem(families["FB_Z"]["Negative"])
        plot_two_conditions(
            time_axis, fb_z_P, "Feedback +", fb_z_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "feedback_z.png",
            vlines=VERT_LINES_FEEDBACK,
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
            _fmt(res.tp_means_psc.get("ANT_REWARD_WU6")),
            _fmt(res.tp_means_psc.get("ANT_NEUTRAL_WU6")),
            _fmt(res.tp_means_psc.get("FB_POS_REWARD")),
            _fmt(res.tp_means_psc.get("FB_NEG_REWARD")),
            _fmt(res.tp_means_psc.get("FB_POS_NEUTRAL")),
            _fmt(res.tp_means_psc.get("FB_NEG_NEUTRAL")),
            _fmt(res.tp_means_z.get("ANT_REWARD")),
            _fmt(res.tp_means_z.get("ANT_NEUTRAL")),
            _fmt(res.tp_means_z.get("ANT_REWARD_WU6")),
            _fmt(res.tp_means_z.get("ANT_NEUTRAL_WU6")),
            _fmt(res.tp_means_z.get("FB_POS_REWARD")),
            _fmt(res.tp_means_z.get("FB_NEG_REWARD")),
            _fmt(res.tp_means_z.get("FB_POS_NEUTRAL")),
            _fmt(res.tp_means_z.get("FB_NEG_NEUTRAL")),
            str(res.tp_counts.get("ANT_REWARD",0)),
            str(res.tp_counts.get("ANT_NEUTRAL",0)),
            str(res.tp_counts.get("ANT_REWARD_WU6",0)),
            str(res.tp_counts.get("ANT_NEUTRAL_WU6",0)),
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
        "ANT_REWARD_PSC","ANT_NEUTRAL_PSC","ANT_REWARD_PSC_WU6","ANT_NEUTRAL_PSC_WU6",
        "FB_POS_REWARD_PSC","FB_NEG_REWARD_PSC","FB_POS_NEUTRAL_PSC","FB_NEG_NEUTRAL_PSC",
        "ANT_REWARD_Z","ANT_NEUTRAL_Z","ANT_REWARD_Z_WU6","ANT_NEUTRAL_Z_WU6",
        "FB_POS_REWARD_Z","FB_NEG_REWARD_Z","FB_POS_NEUTRAL_Z","FB_NEG_NEUTRAL_Z",
        "N_ANT_REWARD","N_ANT_NEUTRAL","N_ANT_REWARD_WU6","N_ANT_NEUTRAL_WU6","N_FB_POS_REWARD","N_FB_NEG_REWARD","N_FB_POS_NEUTRAL","N_FB_NEG_NEUTRAL"
    ]
    tsv_path = SUMMARY_DIR / "summary_at_6s_mid-unsmoothed_interpolated.tsv"

    with open(tsv_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in rows_for_tp:
            f.write("\t".join(row) + "\n")


    print(f"Done. Outputs in: {OUT_TC_DIR}\n  - Discrete summaries: {tsv_path}")


if __name__ == "__main__":
    main()
