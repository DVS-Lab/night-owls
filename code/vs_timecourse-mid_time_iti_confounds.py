#!/usr/bin/env python3
"""
VS time-course extraction, plotting, and *time-based* summaries (seconds everywhere)
-------------------------------------------------------------------------------
Adds (v3)
- Subject-level ANT plots stratified by ITI bins (short / medium / long).
- Keeps the prior ITI-threshold option (MIN_ITI_S) as an optional additional filter.
- Confound regression still supported (confoundevs.txt).
- Robust to empty bins / empty conditions (no plotting crashes).

What "ITI" means here
- For each anticipation cue onset on trial t:
      ITI(t) = cue_onset(t) - feedback_offset(t-1)
  where feedback_offset = feedback_onset + feedback_duration.
- ITIs are inferred per run from the same EV timing files used by FEAT.

Important
- ITI binning is intended mainly for *subject-level* inspection.
  Run-level ANT plots are disabled by default to avoid sparse, misleading figures.

Outputs
  <rootdir>/derivatives/extractions/timecourses-mid-unsmoothed/
    ├── runs/
    └── subjects/
  <rootdir>/derivatives/extractions/summary_at_times_mid-unsmoothed.tsv
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

OUT_TC_DIR  = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-unsmoothed"
SUMMARY_DIR = ROOT_DIR / "derivatives" / "extractions"

EV_BASE = FSL_DERIV / "EVFiles"
VS_MNI = MASKS_DIR / "space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-unsmoothed.txt"

# ----------------------------- Analysis parameters -----------------------------
TASK = "mid"

# Peri-event window (seconds)
TMIN = -4.0
TMAX = 16.0

# Plot sampling grid (seconds). Independent of TR.
DT = 0.10  # 100 ms grid

# Timing constants
CUE_DUR_S = 0.750  # cue offset = cue onset + 0.75

# Summary times (seconds from event anchor)
SUMMARY_OFFSETS_S = [6.0]

# Vertical reference markers (seconds from event anchor)
VERT_LINES_ANT = [0.0, 1.5, 3.0, 6.0]
VERT_LINES_FB  = [0.0, 6.0]

# ----------------------------- Optional controls ----------------------
# Additional ITI threshold for ANT trials (after binning). Set to None to disable.
MIN_ITI_S: Optional[float] = None  # e.g., 8.0 if you want to drop very short ITIs entirely

# Confound regression from confoundevs.txt
APPLY_CONFOUND_REGRESSION = True

# Plotting
PLOT_RUN_LEVEL_FB = True
PLOT_RUN_LEVEL_ANT = False  # keep OFF; subject-level ITI bin plots are the goal

# ITI binning configuration (subject-level)
# Two common strategies:
#   - "fixed": bins use fixed cutpoints (recommended if you want interpretability across subjects)
#   - "tertiles": bins are subject-specific tertiles (recommended if some subjects have skewed ITI distributions)
ITI_BIN_METHOD = "fixed"  # "fixed" or "tertiles"

# If ITI_BIN_METHOD == "fixed", these edges define [short, medium, long]
# Example: short < 4s, medium 4-7s, long >= 7s
ITI_FIXED_EDGES = (4.0, 7.0)

ITI_BIN_LABELS = ("short", "medium", "long")

# Minimum number of trials per (bin x condition) to draw a plot.
# If fewer, we still compute summaries, but the curves can be very unstable.
MIN_TRIALS_TO_PLOT = 10

# Plot style
plt.rcParams.update({
    "figure.dpi": 120,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 12,
})

# ------------------------------ Data structures --------------------------------
@dataclass
class RunResult:
    feat_path: Path
    sub: str
    ses: str
    run: str
    echo: str
    time_axis: np.ndarray

    # Overall ANT/FB curves (not binned)
    ant_psc: Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    ant_z:   Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_psc:  Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_z:    Dict[str, Tuple[np.ndarray, np.ndarray, int]]

    # ITIs (cue-onset anchored) for ANT events, for binning downstream
    iti_reward: np.ndarray
    iti_neutral: np.ndarray

    # Binned ANT windows (delay-onset anchored; after +CUE_DUR_S shift)
    # Keys: bin -> condition -> W array (n_trials, n_time)
    ant_bins_psc_W: Dict[str, Dict[str, np.ndarray]]
    ant_bins_z_W:   Dict[str, Dict[str, np.ndarray]]

    # Point estimates and counts
    points_psc: Dict[str, float | None]
    points_z:   Dict[str, float | None]
    points_n:   Dict[str, int]

    # QC counts
    n_ant_total: int

# ------------------------------ Helper functions --------------------------------
def parse_sub_ses_run_from_feat(feat: Path) -> Tuple[str, str, str]:
    m_sub = re.search(r"sub-(\d+)", feat.as_posix())
    m_ses = re.search(r"ses-(\d+)", feat.as_posix())
    m_run = re.search(r"run-(\d+)", feat.as_posix())
    if not (m_sub and m_ses and m_run):
        raise ValueError(f"Could not parse sub/ses/run from: {feat}")
    return m_sub.group(1), m_ses.group(1), m_run.group(1)


def echo_from_feat(feat: Path) -> str:
    s = feat.as_posix().lower()
    if "tedana" in s or "multi-echo" in s:
        return "multi-echo"
    if "fmriprep" in s or "single-echo" in s:
        return "single-echo"
    return "unknown-echo"


def get_ev_dir(sub: str, ses: str, run: str) -> Path:
    return EV_BASE / f"sub-{sub}" / f"ses-{ses}" / TASK / f"run-{run}"


def load_ev(ev_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load FSL EV (onset, duration, amplitude). Returns (onsets, durations)."""
    if not ev_path.exists():
        return np.array([], dtype=float), np.array([], dtype=float)
    try:
        arr = np.loadtxt(ev_path, ndmin=2)
    except Exception:
        return np.array([], dtype=float), np.array([], dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    return np.asarray(arr[:, 0], dtype=float), np.asarray(arr[:, 1], dtype=float)


def load_vs_timeseries(feat: Path) -> Tuple[np.ndarray, float]:
    func_path = feat / "filtered_func_data.nii.gz"
    if not func_path.exists():
        raise FileNotFoundError(f"Missing: {func_path}")
    if not VS_MNI.exists():
        raise FileNotFoundError(f"Missing mask: {VS_MNI}")

    img = nib.load(str(func_path))
    tr = float(img.header.get_zooms()[3])
    data = img.get_fdata(dtype=np.float32)

    mask_img = nib.load(str(VS_MNI))
    mask = mask_img.get_fdata().astype(bool)
    if mask.shape != data.shape[:3]:
        raise ValueError(f"Mask shape {mask.shape} != data shape {data.shape[:3]} for {func_path}")

    ts = data[mask, :].mean(axis=0)
    return ts.astype(float), tr


def load_confound_matrix(feat: Path, n_tp: int) -> Optional[np.ndarray]:
    for p in [feat / "confoundevs.txt", feat / "confoundEVs.txt", feat / "confounds.txt"]:
        if p.exists():
            try:
                X = np.loadtxt(p, ndmin=2)
            except Exception:
                return None
            if X.shape[0] != n_tp:
                return None
            return X.astype(float)
    return None


def regress_confounds(ts: np.ndarray, X: np.ndarray) -> np.ndarray:
    if X.size == 0:
        return ts
    Xc = X - np.nanmean(X, axis=0, keepdims=True)
    good = np.isfinite(ts) & np.all(np.isfinite(Xc), axis=1)
    if np.sum(good) < max(10, Xc.shape[1] + 2):
        return ts
    y = ts[good]
    Xg = Xc[good, :]
    A = np.column_stack([np.ones(Xg.shape[0]), Xg])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    conf_part = Xc @ beta[1:]
    return ts - conf_part


def to_psc(ts: np.ndarray) -> np.ndarray:
    mu = float(np.nanmean(ts))
    if mu == 0 or not np.isfinite(mu):
        return np.zeros_like(ts)
    return ((ts / mu) - 1.0) * 100.0


def to_z(ts: np.ndarray) -> np.ndarray:
    mu = float(np.nanmean(ts))
    sd = float(np.nanstd(ts, ddof=1))
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(ts)
    return (ts - mu) / sd


def build_time_axis(tmin: float, tmax: float, dt: float) -> np.ndarray:
    n = int(np.floor((tmax - tmin) / dt)) + 1
    return (tmin + np.arange(n) * dt).astype(float)


def sample_windows(ts: np.ndarray, onsets_s: np.ndarray, tr: float, t_axis: np.ndarray) -> np.ndarray:
    """Return (n_trials, n_time) sampled by linear interpolation."""
    if onsets_s.size == 0:
        return np.empty((0, t_axis.size), dtype=float)
    t_series = np.arange(ts.size) * tr
    W = np.empty((onsets_s.size, t_axis.size), dtype=float)
    W[:] = np.nan
    for i, onset in enumerate(onsets_s):
        t_abs = onset + t_axis
        W[i, :] = np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)
    return W


def mean_and_sem(W: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Always returns vectors of length W.shape[1] (even if 0 trials), so plotting never crashes.
    SEM is computed safely:
      - if n<2 at a timepoint, SEM=0 at that timepoint.
    """
    if W.ndim != 2:
        raise ValueError("W must be 2D (n_trials, n_time).")
    n_time = W.shape[1]
    if W.shape[0] == 0:
        return np.full(n_time, np.nan), np.full(n_time, np.nan), 0

    mean = np.nanmean(W, axis=0)
    n = np.sum(np.isfinite(W), axis=0)

    sd = np.full(n_time, 0.0)
    for j in range(n_time):
        if n[j] >= 2:
            sd[j] = float(np.nanstd(W[:, j], ddof=1))
        else:
            sd[j] = 0.0

    sem = sd / np.sqrt(np.maximum(n, 1))
    n_trials = int(np.sum(~np.all(np.isnan(W), axis=1)))
    return mean, sem, n_trials


def vals_at_offset(ts: np.ndarray, onsets_s: np.ndarray, tr: float, offset_s: float) -> np.ndarray:
    if onsets_s.size == 0:
        return np.array([], dtype=float)
    t_series = np.arange(ts.size) * tr
    t_abs = onsets_s + float(offset_s)
    return np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)


def compute_itis_from_feedback_offsets(cue_onsets: np.ndarray, fb_offsets_sorted: np.ndarray) -> np.ndarray:
    """
    For each cue onset, find the most recent feedback offset strictly before it and compute ITI.
    ITI = cue_onset - prev_feedback_offset
    If none exists, ITI is NaN.
    """
    if cue_onsets.size == 0:
        return np.array([], dtype=float)
    if fb_offsets_sorted.size == 0:
        return np.full_like(cue_onsets, np.nan, dtype=float)

    idx = np.searchsorted(fb_offsets_sorted, cue_onsets, side="left") - 1
    itis = np.full_like(cue_onsets, np.nan, dtype=float)
    valid = idx >= 0
    itis[valid] = cue_onsets[valid] - fb_offsets_sorted[idx[valid]]
    return itis


def bin_edges_for_subject(itis_all: np.ndarray) -> Tuple[float, float]:
    """
    Returns (edge1, edge2) such that:
      short < edge1
      medium [edge1, edge2)
      long >= edge2
    """
    itis = itis_all[np.isfinite(itis_all)]
    if itis.size < 30:
        # too few to estimate stable tertiles; fall back to fixed
        return ITI_FIXED_EDGES

    if ITI_BIN_METHOD == "tertiles":
        e1, e2 = np.quantile(itis, [1/3, 2/3])
        # Guard against degenerate quantiles
        if not np.isfinite(e1) or not np.isfinite(e2) or e2 <= e1:
            return ITI_FIXED_EDGES
        return float(e1), float(e2)

    # fixed
    return ITI_FIXED_EDGES


def assign_bins(itis: np.ndarray, edge1: float, edge2: float) -> Dict[str, np.ndarray]:
    """Return dict bin_label -> boolean mask for that bin (NaNs go to none)."""
    out: Dict[str, np.ndarray] = {}
    finite = np.isfinite(itis)
    short = finite & (itis < edge1)
    medium = finite & (itis >= edge1) & (itis < edge2)
    long = finite & (itis >= edge2)
    out[ITI_BIN_LABELS[0]] = short
    out[ITI_BIN_LABELS[1]] = medium
    out[ITI_BIN_LABELS[2]] = long
    return out


def plot_two_conditions(time_axis: np.ndarray,
                        condA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
                        condB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
                        title: str, ylabel: str, out_png: Path,
                        vlines: List[float]) -> None:
    mA, sA, nA = condA
    mB, sB, nB = condB
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(time_axis, mA, label=f"{labelA} (n={nA})")
    ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)
    ax.plot(time_axis, mB, label=f"{labelB} (n={nB})")
    ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)
    for v in vlines:
        ax.axvline(v, ls=":", lw=1)
    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def _fmt(x: float | None) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return f"{x:.6f}"

# ------------------------------- Core pipeline ----------------------------------
def process_one_feat(feat: Path, time_axis: np.ndarray) -> RunResult | None:
    try:
        sub, ses, run = parse_sub_ses_run_from_feat(feat)
        echo = echo_from_feat(feat)
        ts_raw, tr = load_vs_timeseries(feat)
    except Exception as e:
        print(f"[WARN] Skipping FEAT due to error: {feat}\n  -> {e}")
        return None

    if APPLY_CONFOUND_REGRESSION:
        X = load_confound_matrix(feat, n_tp=ts_raw.size)
        if X is not None:
            ts_raw = regress_confounds(ts_raw, X)

    ts_psc = to_psc(ts_raw)
    ts_z   = to_z(ts_raw)

    ev_dir = get_ev_dir(sub, ses, run)

    # ANT cue onsets (cue-onset anchored)
    ant_R_on, _ = load_ev(ev_dir / "_anticipation_reward.txt")
    ant_N_on, _ = load_ev(ev_dir / "_anticipation_neutral.txt")

    # Feedback events (for offsets / ITI inference)
    fb_PR_on, fb_PR_dur = load_ev(ev_dir / "_feedback_positive_reward.txt")
    fb_NR_on, fb_NR_dur = load_ev(ev_dir / "_feedback_negative_reward.txt")
    fb_PN_on, fb_PN_dur = load_ev(ev_dir / "_feedback_positive_neutral.txt")
    fb_NN_on, fb_NN_dur = load_ev(ev_dir / "_feedback_negative_neutral.txt")

    fb_pos = np.sort(np.concatenate([fb_PR_on, fb_PN_on]))
    fb_neg = np.sort(np.concatenate([fb_NR_on, fb_NN_on]))

    fb_all_on  = np.concatenate([fb_PR_on, fb_NR_on, fb_PN_on, fb_NN_on])
    fb_all_dur = np.concatenate([fb_PR_dur, fb_NR_dur, fb_PN_dur, fb_NN_dur])
    fb_offsets = np.sort(fb_all_on + fb_all_dur)

    iti_R = compute_itis_from_feedback_offsets(ant_R_on, fb_offsets)
    iti_N = compute_itis_from_feedback_offsets(ant_N_on, fb_offsets)

    # Optional additional ITI threshold (after computing ITIs)
    if MIN_ITI_S is not None:
        keep_R = np.isfinite(iti_R) & (iti_R >= MIN_ITI_S)
        keep_N = np.isfinite(iti_N) & (iti_N >= MIN_ITI_S)
        ant_R_on, iti_R = ant_R_on[keep_R], iti_R[keep_R]
        ant_N_on, iti_N = ant_N_on[keep_N], iti_N[keep_N]

    n_ant_total = int(ant_R_on.size + ant_N_on.size)

    # Wu-style ANT anchoring (delay onset = cue offset)
    ant_R_delay = ant_R_on + CUE_DUR_S
    ant_N_delay = ant_N_on + CUE_DUR_S

    # Overall (not binned) ANT/FB windows
    R_psc_w = sample_windows(ts_psc, ant_R_delay, tr, time_axis)
    N_psc_w = sample_windows(ts_psc, ant_N_delay, tr, time_axis)
    R_z_w   = sample_windows(ts_z,   ant_R_delay, tr, time_axis)
    N_z_w   = sample_windows(ts_z,   ant_N_delay, tr, time_axis)

    P_psc_w = sample_windows(ts_psc, fb_pos, tr, time_axis)
    G_psc_w = sample_windows(ts_psc, fb_neg, tr, time_axis)
    P_z_w   = sample_windows(ts_z,   fb_pos, tr, time_axis)
    G_z_w   = sample_windows(ts_z,   fb_neg, tr, time_axis)

    ant_psc = {"Reward": mean_and_sem(R_psc_w), "Neutral": mean_and_sem(N_psc_w)}
    ant_z   = {"Reward": mean_and_sem(R_z_w),   "Neutral": mean_and_sem(N_z_w)}
    fb_psc  = {"Positive": mean_and_sem(P_psc_w), "Negative": mean_and_sem(G_psc_w)}
    fb_z    = {"Positive": mean_and_sem(P_z_w),   "Negative": mean_and_sem(G_z_w)}

    # Binned ANT windows (store trial-level W for subject-level concatenation)
    ant_bins_psc_W: Dict[str, Dict[str, np.ndarray]] = {b: {"Reward": np.empty((0, time_axis.size)),
                                                           "Neutral": np.empty((0, time_axis.size))}
                                                       for b in ITI_BIN_LABELS}
    ant_bins_z_W:   Dict[str, Dict[str, np.ndarray]] = {b: {"Reward": np.empty((0, time_axis.size)),
                                                           "Neutral": np.empty((0, time_axis.size))}
                                                       for b in ITI_BIN_LABELS}

    # We *do not* choose edges here (subject-level edges are better). We store ITIs + trial windows.
    # We'll assign bins at subject-aggregation time.

    # -------------------- Time-based point estimates (interpolated) --------------------
    def m(arr: np.ndarray) -> float | None:
        if arr.size == 0 or not np.any(np.isfinite(arr)):
            return None
        return float(np.nanmean(arr))

    def n_valid(arr: np.ndarray) -> int:
        return int(np.sum(np.isfinite(arr))) if arr.size else 0

    points_psc: Dict[str, float | None] = {}
    points_z:   Dict[str, float | None] = {}
    points_n:   Dict[str, int] = {}

    for off in SUMMARY_OFFSETS_S:
        tag = f"{off:.3f}".rstrip("0").rstrip(".")
        ant_R_psc = vals_at_offset(ts_psc, ant_R_delay, tr, off)
        ant_N_psc = vals_at_offset(ts_psc, ant_N_delay, tr, off)
        ant_R_z   = vals_at_offset(ts_z,   ant_R_delay, tr, off)
        ant_N_z   = vals_at_offset(ts_z,   ant_N_delay, tr, off)

        points_psc[f"ANT_REWARD_t{tag}"]  = m(ant_R_psc)
        points_psc[f"ANT_NEUTRAL_t{tag}"] = m(ant_N_psc)
        points_z[f"ANT_REWARD_t{tag}"]    = m(ant_R_z)
        points_z[f"ANT_NEUTRAL_t{tag}"]   = m(ant_N_z)
        points_n[f"N_ANT_REWARD_t{tag}"]  = n_valid(ant_R_psc)
        points_n[f"N_ANT_NEUTRAL_t{tag}"] = n_valid(ant_N_psc)

        for name, on in [
            ("FB_POS_REWARD", fb_PR_on),
            ("FB_NEG_REWARD", fb_NR_on),
            ("FB_POS_NEUTRAL", fb_PN_on),
            ("FB_NEG_NEUTRAL", fb_NN_on),
        ]:
            v_psc = vals_at_offset(ts_psc, on, tr, off)
            v_z   = vals_at_offset(ts_z,   on, tr, off)
            points_psc[f"{name}_t{tag}"] = m(v_psc)
            points_z[f"{name}_t{tag}"]   = m(v_z)
            points_n[f"N_{name}_t{tag}"] = n_valid(v_psc)

    # Store the trial-level windows so we can bin at the subject-level
    ant_bins_psc_W["__ALL__"] = {"Reward": R_psc_w, "Neutral": N_psc_w}
    ant_bins_z_W["__ALL__"]   = {"Reward": R_z_w,   "Neutral": N_z_w}

    return RunResult(
        feat_path=feat,
        sub=sub, ses=ses, run=run, echo=echo,
        time_axis=time_axis,
        ant_psc=ant_psc, ant_z=ant_z,
        fb_psc=fb_psc, fb_z=fb_z,
        iti_reward=iti_R,
        iti_neutral=iti_N,
        ant_bins_psc_W=ant_bins_psc_W,
        ant_bins_z_W=ant_bins_z_W,
        points_psc=points_psc,
        points_z=points_z,
        points_n=points_n,
        n_ant_total=n_ant_total,
    )


def save_run_plots(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

    if PLOT_RUN_LEVEL_ANT:
        plot_two_conditions(
            res.time_axis,
            res.ant_psc["Reward"], "Reward",
            res.ant_psc["Neutral"], "Neutral",
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=run_out / "anticipation_psc.png",
            vlines=VERT_LINES_ANT,
        )

    if PLOT_RUN_LEVEL_FB:
        plot_two_conditions(
            res.time_axis,
            res.fb_psc["Positive"], "Feedback +",
            res.fb_psc["Negative"], "Feedback −",
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=run_out / "feedback_psc.png",
            vlines=VERT_LINES_FB,
        )
        plot_two_conditions(
            res.time_axis,
            res.fb_z["Positive"], "Feedback +",
            res.fb_z["Negative"], "Feedback −",
            title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, Z)",
            ylabel="Z (SD units)",
            out_png=run_out / "feedback_z.png",
            vlines=VERT_LINES_FB,
        )


def concat_windows(w_list: List[np.ndarray], n_time: int) -> np.ndarray:
    if not w_list:
        return np.empty((0, n_time), dtype=float)
    return np.concatenate(w_list, axis=0)


def aggregate_subject(results: List[RunResult], subject: str, time_axis: np.ndarray) -> Dict[str, Dict[str, Dict[str, int]]]:
    """
    Build subject-level plots and return bin counts for summary.
    Returns counts[echo][bin][condition] = n_trials
    """
    subj_out = OUT_TC_DIR / "subjects" / f"sub-{subject}"
    n_time = time_axis.size

    # Organize by echo
    by_echo: Dict[str, List[RunResult]] = {}
    for r in results:
        if r.sub == subject:
            by_echo.setdefault(r.echo, []).append(r)

    counts_out: Dict[str, Dict[str, Dict[str, int]]] = {}

    for echo, runs in by_echo.items():
        # Collect all ITIs for this subject+echo to define bin edges
        itis_all = np.concatenate([r.iti_reward for r in runs] + [r.iti_neutral for r in runs])
        edge1, edge2 = bin_edges_for_subject(itis_all)

        # Collect trial windows + ITIs per condition (these are aligned to delay onset)
        R_W_psc = []
        N_W_psc = []
        R_W_z   = []
        N_W_z   = []

        R_iti = []
        N_iti = []

        for r in runs:
            # "__ALL__" contains the full trial windows for this run
            R_W_psc.append(r.ant_bins_psc_W["__ALL__"]["Reward"])
            N_W_psc.append(r.ant_bins_psc_W["__ALL__"]["Neutral"])
            R_W_z.append(r.ant_bins_z_W["__ALL__"]["Reward"])
            N_W_z.append(r.ant_bins_z_W["__ALL__"]["Neutral"])
            R_iti.append(r.iti_reward)
            N_iti.append(r.iti_neutral)

        R_W_psc_all = concat_windows(R_W_psc, n_time)
        N_W_psc_all = concat_windows(N_W_psc, n_time)
        R_W_z_all   = concat_windows(R_W_z, n_time)
        N_W_z_all   = concat_windows(N_W_z, n_time)

        R_iti_all = np.concatenate(R_iti) if R_iti else np.array([], dtype=float)
        N_iti_all = np.concatenate(N_iti) if N_iti else np.array([], dtype=float)

        # Assign bins separately by condition using the same edges
        bins_R = assign_bins(R_iti_all, edge1, edge2)
        bins_N = assign_bins(N_iti_all, edge1, edge2)

        counts_out.setdefault(echo, {})
        for b in ITI_BIN_LABELS:
            # Slice windows by bin masks (must align: ITIs and W rows correspond to the same events)
            mR = bins_R[b] if bins_R[b].size == R_W_psc_all.shape[0] else np.zeros(R_W_psc_all.shape[0], dtype=bool)
            mN = bins_N[b] if bins_N[b].size == N_W_psc_all.shape[0] else np.zeros(N_W_psc_all.shape[0], dtype=bool)

            W_R_psc = R_W_psc_all[mR, :] if R_W_psc_all.size else np.empty((0, n_time))
            W_N_psc = N_W_psc_all[mN, :] if N_W_psc_all.size else np.empty((0, n_time))
            W_R_z   = R_W_z_all[mR, :]   if R_W_z_all.size else np.empty((0, n_time))
            W_N_z   = N_W_z_all[mN, :]   if N_W_z_all.size else np.empty((0, n_time))

            cR = W_R_psc.shape[0]
            cN = W_N_psc.shape[0]
            counts_out[echo].setdefault(b, {"Reward": int(cR), "Neutral": int(cN)})

            # Build mean/sem
            R_psc = mean_and_sem(W_R_psc)
            N_psc = mean_and_sem(W_N_psc)
            R_z   = mean_and_sem(W_R_z)
            N_z   = mean_and_sem(W_N_z)

            # Only plot if there's a minimally stable number of trials in both conditions
            if min(cR, cN) < MIN_TRIALS_TO_PLOT:
                continue

            plot_two_conditions(
                time_axis, R_psc, "Reward", N_psc, "Neutral",
                title=(f"VS — subject {subject} [{echo}] (ANT, PSC) | ITI bin: {b} "
                       f"({ITI_BIN_METHOD}; edges={edge1:.2f},{edge2:.2f})"),
                ylabel="% signal change (PSC)",
                out_png=subj_out / echo / f"anticipation_psc_iti-{b}.png",
                vlines=VERT_LINES_ANT,
            )
            plot_two_conditions(
                time_axis, R_z, "Reward", N_z, "Neutral",
                title=(f"VS — subject {subject} [{echo}] (ANT, Z) | ITI bin: {b} "
                       f"({ITI_BIN_METHOD}; edges={edge1:.2f},{edge2:.2f})"),
                ylabel="Z (SD units)",
                out_png=subj_out / echo / f"anticipation_z_iti-{b}.png",
                vlines=VERT_LINES_ANT,
            )

        # Also save the overall (all ITIs) subject-level curves (not binned)
        R_all = mean_and_sem(R_W_psc_all)
        N_all = mean_and_sem(N_W_psc_all)
        plot_two_conditions(
            time_axis, R_all, "Reward", N_all, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, PSC) | all ITIs",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "anticipation_psc_allITIs.png",
            vlines=VERT_LINES_ANT,
        )

    return counts_out


def main() -> None:
    if not FEAT_LIST_PATH.exists():
        print(f"[ERROR] FEAT list not found: {FEAT_LIST_PATH}")
        return

    feat_paths: List[Path] = []
    for line in FEAT_LIST_PATH.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        feat_paths.append(Path(line))

    if not feat_paths:
        print("[ERROR] No FEAT paths found in list.")
        return

    OUT_TC_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    time_axis = build_time_axis(TMIN, TMAX, DT)

    results: List[RunResult] = []
    echo_counts: Dict[str, int] = {}

    for feat in feat_paths:
        res = process_one_feat(feat, time_axis)
        if res is None:
            continue
        save_run_plots(res)
        results.append(res)
        echo_counts[res.echo] = echo_counts.get(res.echo, 0) + 1

    print("Processed runs by echo:", echo_counts)
    print(f"ITI binning method: {ITI_BIN_METHOD}. Fixed edges (if used): {ITI_FIXED_EDGES}.")
    if MIN_ITI_S is not None:
        print(f"Additional ANT ITI threshold active: ITI >= {MIN_ITI_S} s (applied before binning).")

    if not results:
        print("[ERROR] No runs processed successfully.")
        return

    # Aggregate by subject + write out bin counts
    subjects = sorted({r.sub for r in results})
    subject_bin_counts: Dict[str, Dict[str, Dict[str, Dict[str, int]]]] = {}
    for sub in subjects:
        subject_bin_counts[sub] = aggregate_subject(results, sub, time_axis)

    # Write summary table (run-level rows) + include run-level ANT total counts
    header = ["sub", "ses", "run", "echo", "n_ant_total"]
    keys_psc = sorted({k for r in results for k in r.points_psc.keys()})
    keys_z   = sorted({k for r in results for k in r.points_z.keys()})
    keys_n   = sorted({k for r in results for k in r.points_n.keys()})
    header += keys_psc + keys_z + keys_n

    # Add subject-level bin counts columns (will repeat per run; still useful for quick lookup)
    # Columns: SUBJ_ANT_N_<echo>_<bin>_<cond>
    subj_cols: List[str] = []
    echos = sorted({r.echo for r in results})
    for e in echos:
        for b in ITI_BIN_LABELS:
            for c in ("Reward", "Neutral"):
                subj_cols.append(f"SUBJ_ANT_N_{e}_{b}_{c}")
    header += subj_cols

    out_tsv = SUMMARY_DIR / "summary_at_times_mid-unsmoothed.tsv"
    with open(out_tsv, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in results:
            row = [r.sub, r.ses, r.run, r.echo, str(r.n_ant_total)]
            row += [_fmt(r.points_psc.get(k)) for k in keys_psc]
            row += [_fmt(r.points_z.get(k)) for k in keys_z]
            row += [str(r.points_n.get(k, 0)) for k in keys_n]

            # subject-level bin counts (may be missing if nothing plotted/available)
            sc = subject_bin_counts.get(r.sub, {})
            for col in subj_cols:
                # parse the column
                # SUBJ_ANT_N_<echo>_<bin>_<cond>
                parts = col.split("_")
                e = parts[3]
                b = parts[4]
                c = parts[5]
                v = sc.get(e, {}).get(b, {}).get(c, 0)
                row.append(str(v))

            f.write("\t".join(row) + "\n")

    print(f"Done.\n  Plots: {OUT_TC_DIR}\n  Summaries: {out_tsv}")
    print("Subject-level ITI-bin plots are saved under: derivatives/extractions/timecourses-mid-unsmoothed/subjects/")


if __name__ == "__main__":
    main()
