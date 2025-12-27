#!/usr/bin/env python3
"""
VS time-course extraction, plotting, and *time-based* summaries (seconds everywhere)
-------------------------------------------------------------------------------
Why this version exists
- Avoids "Nth TR" thinking in summaries and QC markers.
- Plots and point-estimates are defined in **seconds**, then sampled from the
  TR-sampled time series via linear interpolation.

Key alignment choices (Wu comparability without re-running FEAT)
- ANTICIPATION (Reward vs Neutral):
    * Your EV1/EV2 onsets are cue-onset; durations are cue+ISI.
    * For Wu-style NAcc/VS "anticipation" (delay onset), we shift ANT onsets by:
          +CUE_DUR_S
      so t=0 in ANT plots/metrics corresponds to cue offset / delay onset.
- FEEDBACK (Positive vs Negative):
    * Uses feedback-onset EVs as-is; t=0 is feedback onset.

What this script outputs
- Per-run and per-subject (echo-split) plots for:
    ANT (Reward vs Neutral) and FB (Positive vs Negative), both in PSC and Z.
- A TSV summary with point-estimates sampled at explicit times (e.g., +6.0 s),
  using interpolation so the requested time is hit exactly.

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
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses"
SUMMARY_DIR = ROOT_DIR / "derivatives" / "extractions"

# Keep your existing conventions from the original script
FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-unsmoothed.txt"
VS_MNI         = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"

# ----------------------------- Analysis parameters -----------------------------
TASK = "mid"

# Peri-event window (seconds)
TMIN = -4.0
TMAX = 16.0

# Plot sampling grid (seconds). This is independent of TR.
# Under the hood, values are sampled from the TR time series via interpolation.
DT = 0.10  # 100 ms grid is a nice QC resolution without being huge

# Timing constants for your MID
CUE_DUR_S = 0.750  # fixed cue duration in your task

# If design.fsf TR disagrees with NIfTI header TR, you can optionally prefer design.fsf.
PREFER_TR_FROM_DESIGN = False

# Point-estimate times (seconds from event anchor). Add more if you want.
# - 6.0 s is commonly used as an HRF-lag point estimate, and is used in Wu for NAcc (delay onset + 6 s).
SUMMARY_OFFSETS_S = [6.0]

# Vertical reference markers (seconds from event anchor)
# ANT is delay-onset locked (after shifting by +0.75 s), so target occurs 1.5–3.0 s later (ISI).
VERT_LINES_ANT = [0.0, 1.5, 3.0, 6.0]
# FB is feedback-onset locked
VERT_LINES_FB  = [0.0, 6.0]

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
    echo: str  # 'single-echo' or 'multi-echo'
    time_axis: np.ndarray
    # TR sanity-check (header vs design.fsf)
    tr_header: float
    tr_design: float | None
    tr_used: float
    # Timing diagnostics (per-trial, seconds)
    isi_s: np.ndarray  # target-onset relative to cue offset, from ANT EV durations
    fb_delay_s: np.ndarray  # feedback-onset relative to cue offset (nearest following feedback)
    # Each dict maps condition label -> (mean, sem, n_trials)
    ant_psc: Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    ant_z:   Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_psc:  Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_z:    Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    # Time-based point estimates
    points_psc: Dict[str, float | None]
    points_z:   Dict[str, float | None]
    points_n:   Dict[str, int]

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
    # Heuristic: tweak to match your naming if needed
    if "tedana" in s or "multi-echo" in s or "me" in s:
        return "multi-echo"
    return "single-echo"


def get_ev_dir(sub: str, ses: str, run: str) -> Path:
    # <root>/derivatives/fsl/EVFiles/sub-XXX/ses-YY/mid/run-ZZ/
    return FSL_DERIV / "EVFiles" / f"sub-{sub}" / f"ses-{ses}" / TASK / f"run-{run}"


def load_ev_3col(ev_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load a standard 3-column FSL EV file: onset(s), duration(s), weight.

    Returns empty arrays if the file is missing or empty.
    """
    if not ev_path.exists():
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
    try:
        arr = np.loadtxt(ev_path, ndmin=2)
    except Exception:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
    if arr.shape[1] < 3:
        raise ValueError(f"EV file does not have 3 columns: {ev_path}")
    on = np.asarray(arr[:, 0], dtype=float)
    dur = np.asarray(arr[:, 1], dtype=float)
    w = np.asarray(arr[:, 2], dtype=float)
    return on, dur, w


def load_ev_onsets(ev_path: Path) -> np.ndarray:
    """Convenience: load only onsets from a 3-col EV file."""
    on, _, _ = load_ev_3col(ev_path)
    return on


def load_tr_from_design_fsf(feat: Path) -> float | None:
    """Extract TR (seconds) from FEAT's design.fsf if available."""
    fsf = feat / "design.fsf"
    if not fsf.exists():
        return None
    try:
        txt = fsf.read_text(errors="ignore")
    except Exception:
        return None
    m = re.search(r"set\s+fmri\(tr\)\s+([0-9]*\.?[0-9]+)", txt)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None
def load_vs_timeseries(feat: Path) -> Tuple[np.ndarray, float]:
    """Load filtered_func_data and extract mean within VS/NAcc mask."""
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


def to_psc(ts: np.ndarray) -> np.ndarray:
    mu = float(np.mean(ts))
    if mu == 0:
        return np.zeros_like(ts)
    return ((ts / mu) - 1.0) * 100.0


def to_z(ts: np.ndarray) -> np.ndarray:
    mu = float(np.mean(ts))
    sd = float(np.std(ts, ddof=1))
    if sd == 0:
        return np.zeros_like(ts)
    return (ts - mu) / sd


def build_time_axis(tmin: float, tmax: float, dt: float) -> np.ndarray:
    n = int(np.floor((tmax - tmin) / dt)) + 1
    return (tmin + np.arange(n) * dt).astype(float)


def sample_windows(ts: np.ndarray, onsets_s: np.ndarray, tr: float,
                   t_axis: np.ndarray) -> np.ndarray:
    """Return array of shape (n_trials, n_time) sampled by linear interpolation."""
    if onsets_s.size == 0:
        return np.empty((0, t_axis.size), dtype=float)
    t_series = np.arange(ts.size) * tr
    windows = np.empty((onsets_s.size, t_axis.size), dtype=float)
    windows[:] = np.nan
    for i, onset in enumerate(onsets_s):
        t_abs = onset + t_axis
        windows[i, :] = np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)
    return windows


def mean_and_sem(w: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """Return (mean, SEM, n_trials) across trials, NaN-safe."""
    if w.size == 0:
        return np.array([]), np.array([]), 0
    mean = np.nanmean(w, axis=0)
    n = np.sum(np.isfinite(w), axis=0)
    sem = np.nanstd(w, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))
    n_trials = int(np.sum(~np.all(np.isnan(w), axis=1)))
    return mean, sem, n_trials


def vals_at_offset(ts: np.ndarray, onsets_s: np.ndarray, tr: float, offset_s: float) -> np.ndarray:
    """Sample ts at (onset + offset_s) using linear interpolation."""
    if onsets_s.size == 0:
        return np.array([], dtype=float)
    t_series = np.arange(ts.size) * tr
    t_abs = onsets_s + float(offset_s)
    return np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)


def plot_two_conditions(
    time_axis: np.ndarray,
    condA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
    condB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
    title: str, ylabel: str, out_png: Path,
    vlines: List[Tuple[float, str]] | List[float] | None = None,
    spans: List[Tuple[float, float, str]] | None = None,
) -> None:
    """Plot two conditions with optional vertical markers and shaded windows.

    vlines can be:
      - None
      - List[float]
      - List[(x, label)]
    spans is a list of (x0, x1, label) shaded intervals.
    """
    mA, sA, nA = condA
    mB, sB, nB = condB

    fig, ax = plt.subplots(figsize=(8, 4.8))

    # Optional shaded windows first so lines/curves sit on top
    if spans:
        for x0, x1, lab in spans:
            if lab:
                ax.axvspan(x0, x1, alpha=0.10, label=lab)
            else:
                ax.axvspan(x0, x1, alpha=0.10)

    ax.plot(time_axis, mA, label=f"{labelA} (n={nA})")
    ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)

    ax.plot(time_axis, mB, label=f"{labelB} (n={nB})")
    ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)

    # Vertical reference lines + (optional) labels
    if vlines:
        # Normalize to [(x, label), ...]
        norm: List[Tuple[float, str]] = []
        for v in vlines:
            if isinstance(v, (tuple, list)) and len(v) >= 1:
                x = float(v[0])
                lab = str(v[1]) if len(v) >= 2 else ""
            else:
                x = float(v)
                lab = ""
            if np.isfinite(x):
                norm.append((x, lab))

        # Draw lines
        for x, _ in norm:
            ax.axvline(x, ls=":", lw=1)

        # Annotate labels after ylim is set by data
        y_top = ax.get_ylim()[1]
        for x, lab in norm:
            if lab:
                ax.text(x, y_top, lab, rotation=90, va="top", ha="right", fontsize=8)

    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
# ------------------------------- Core pipeline ----------------------------------

def process_one_feat(feat: Path, time_axis: np.ndarray) -> RunResult | None:
    try:
        sub, ses, run = parse_sub_ses_run_from_feat(feat)
        echo = echo_from_feat(feat)
        ts_raw, tr = load_vs_timeseries(feat)
    except Exception as e:
        print(f"[WARN] Skipping FEAT due to error: {feat}\n  -> {e}")
        return None

    # TR sanity-check (header vs design.fsf)
    tr_header = float(tr)
    tr_design = load_tr_from_design_fsf(feat)
    tr_used = tr_header
    if tr_design is not None and abs(tr_design - tr_header) > 1e-3:
        print(f"[WARN] TR mismatch for {feat}: header={tr_header:.6f}s, design.fsf={tr_design:.6f}s")
    if tr_design is not None and PREFER_TR_FROM_DESIGN:
        tr_used = float(tr_design)
    tr = tr_used

    ts_psc = to_psc(ts_raw)
    ts_z   = to_z(ts_raw)

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)

    # Anticipation EVs (cue onset + duration). We shift onsets so t=0 is cue offset
    # to match Wu et al.'s "anticipatory fixation" / delay-onset interpretation.
    ant_R_on, ant_R_dur, _ = load_ev_3col(ev_dir / "_anticipation_reward.txt")
    ant_N_on, ant_N_dur, _ = load_ev_3col(ev_dir / "_anticipation_neutral.txt")

    # Wu-style ANT anchoring (delay onset): shift cue-onset to cue-offset
    ant_R = ant_R_on + CUE_DUR_S
    ant_N = ant_N_on + CUE_DUR_S

    # Timing diagnostics:
    # - target-onset relative to cue-offset ~ (ANT duration - cue duration)
    isi_s = np.concatenate([ant_R_dur, ant_N_dur]) - CUE_DUR_S
    isi_s = isi_s[np.isfinite(isi_s)]

    # Feedback (pooled by valence for plotting)
    fb_pos = np.sort(np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_positive_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_positive_neutral.txt"),
    ]))
    fb_neg = np.sort(np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_negative_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_negative_neutral.txt"),
    ]))

    # Feedback timing diagnostic: for each cue-offset, find the first subsequent feedback onset
    cueoffs = np.sort(np.concatenate([ant_R, ant_N]))
    fb_all = np.sort(np.concatenate([fb_pos, fb_neg]))
    if cueoffs.size == 0 or fb_all.size == 0:
        fb_delay_s = np.array([], dtype=float)
    else:
        idx_next = np.searchsorted(fb_all, cueoffs, side="right")
        valid = idx_next < fb_all.size
        fb_delay_s = fb_all[idx_next[valid]] - cueoffs[valid]
    # Windows (interpolated)
    R_psc_w = sample_windows(ts_psc, ant_R, tr, time_axis)
    N_psc_w = sample_windows(ts_psc, ant_N, tr, time_axis)
    R_z_w   = sample_windows(ts_z,   ant_R, tr, time_axis)
    N_z_w   = sample_windows(ts_z,   ant_N, tr, time_axis)

    P_psc_w = sample_windows(ts_psc, fb_pos, tr, time_axis)
    G_psc_w = sample_windows(ts_psc, fb_neg, tr, time_axis)
    P_z_w   = sample_windows(ts_z,   fb_pos, tr, time_axis)
    G_z_w   = sample_windows(ts_z,   fb_neg, tr, time_axis)

    ant_psc = {"Reward": mean_and_sem(R_psc_w), "Neutral": mean_and_sem(N_psc_w)}
    ant_z   = {"Reward": mean_and_sem(R_z_w),   "Neutral": mean_and_sem(N_z_w)}
    fb_psc  = {"Positive": mean_and_sem(P_psc_w), "Negative": mean_and_sem(G_psc_w)}
    fb_z    = {"Positive": mean_and_sem(P_z_w),   "Negative": mean_and_sem(G_z_w)}

    # -------------------- Time-based point estimates (interpolated) --------------------
    def m(arr: np.ndarray) -> float | None:
        if arr.size == 0:
            return None
        if not np.any(np.isfinite(arr)):
            return None
        return float(np.nanmean(arr))

    def n_valid(arr: np.ndarray) -> int:
        return int(np.sum(np.isfinite(arr))) if arr.size else 0

    points_psc: Dict[str, float | None] = {}
    points_z:   Dict[str, float | None] = {}
    points_n:   Dict[str, int] = {}

    for off in SUMMARY_OFFSETS_S:
        tag = f"{off:.3f}".rstrip("0").rstrip(".")  # pretty label like "6" or "6.5"
        # ANT (delay-onset locked)
        ant_R_psc = vals_at_offset(ts_psc, ant_R, tr, off)
        ant_N_psc = vals_at_offset(ts_psc, ant_N, tr, off)
        ant_R_z   = vals_at_offset(ts_z,   ant_R, tr, off)
        ant_N_z   = vals_at_offset(ts_z,   ant_N, tr, off)

        points_psc[f"ANT_REWARD_t{tag}"]  = m(ant_R_psc)
        points_psc[f"ANT_NEUTRAL_t{tag}"] = m(ant_N_psc)
        points_z[f"ANT_REWARD_t{tag}"]    = m(ant_R_z)
        points_z[f"ANT_NEUTRAL_t{tag}"]   = m(ant_N_z)
        points_n[f"N_ANT_REWARD_t{tag}"]  = n_valid(ant_R_psc)
        points_n[f"N_ANT_NEUTRAL_t{tag}"] = n_valid(ant_N_psc)

        # FB (feedback-onset locked)
        for name, on in [
            ("FB_POS_REWARD", fb_PR),
            ("FB_NEG_REWARD", fb_NR),
            ("FB_POS_NEUTRAL", fb_PN),
            ("FB_NEG_NEUTRAL", fb_NN),
        ]:
            v_psc = vals_at_offset(ts_psc, on, tr, off)
            v_z   = vals_at_offset(ts_z,   on, tr, off)
            points_psc[f"{name}_t{tag}"] = m(v_psc)
            points_z[f"{name}_t{tag}"]   = m(v_z)
            points_n[f"N_{name}_t{tag}"] = n_valid(v_psc)

    return RunResult(
        feat_path=feat,
        sub=sub, ses=ses, run=run, echo=echo,
        time_axis=time_axis,
        tr_header=tr_header, tr_design=tr_design, tr_used=tr_used,
        isi_s=isi_s, fb_delay_s=fb_delay_s,
        ant_psc=ant_psc, ant_z=ant_z,
        fb_psc=fb_psc, fb_z=fb_z,
        points_psc=points_psc,
        points_z=points_z,
        points_n=points_n,
    )


def save_run_plots(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

    def _nanmedian(a: np.ndarray) -> float:
        if a is None or a.size == 0:
            return float("nan")
        return float(np.nanmedian(a))

    # ANT markers (t=0 is cue offset / delay onset)
    med_isi = _nanmedian(res.isi_s)           # ~ target onset
    med_fb = _nanmedian(res.fb_delay_s)       # ~ feedback onset
    ant_vlines: List[Tuple[float, str]] = [
        (-CUE_DUR_S, "Cue onset"),
        (0.0, "Cue offset"),
        (med_isi, "Target (median)"),
        (med_fb, "Feedback (median)"),
        (6.0, "+6s"),
    ]
    ant_spans = [(0.0, 2.0, "Wu ANT window (0–2s)")]

    # FB markers (t=0 is feedback onset)
    fb_vlines: List[Tuple[float, str]] = [
        (0.0, "Feedback onset"),
        (6.0, "+6s"),
    ]

    plot_two_conditions(
        res.time_axis,
        res.ant_psc["Reward"], "Reward",
        res.ant_psc["Neutral"], "Neutral",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, PSC)",
        ylabel="% signal change (PSC)",
        out_png=run_out / "anticipation_psc.png",
        vlines=ant_vlines,
        spans=ant_spans,
    )
    plot_two_conditions(
        res.time_axis,
        res.ant_z["Reward"], "Reward",
        res.ant_z["Neutral"], "Neutral",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (ANT, Z)",
        ylabel="Z (SD units)",
        out_png=run_out / "anticipation_z.png",
        vlines=ant_vlines,
        spans=ant_spans,
    )

    plot_two_conditions(
        res.time_axis,
        res.fb_psc["Positive"], "Feedback +",
        res.fb_psc["Negative"], "Feedback −",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, PSC)",
        ylabel="% signal change (PSC)",
        out_png=run_out / "feedback_psc.png",
        vlines=fb_vlines,
    )
    plot_two_conditions(
        res.time_axis,
        res.fb_z["Positive"], "Feedback +",
        res.fb_z["Negative"], "Feedback −",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (FB, Z)",
        ylabel="Z (SD units)",
        out_png=run_out / "feedback_z.png",
        vlines=fb_vlines,
    )

def aggregate_subject(results: List[RunResult], subject: str, time_axis: np.ndarray) -> None:
    subs = [r for r in results if r.sub == subject]
    if not subs:
        return

    # Group by echo (multi-echo vs single-echo, etc.)
    by_echo: Dict[str, List[RunResult]] = {}
    for r in subs:
        by_echo.setdefault(r.echo, []).append(r)

    def _nanmedian(a: np.ndarray) -> float:
        if a is None or a.size == 0:
            return float("nan")
        return float(np.nanmedian(a))

    for echo, runs in by_echo.items():
        # Pooled timing diagnostics across all runs/sessions for this subject/echo
        isi_all = np.concatenate([r.isi_s for r in runs if r.isi_s is not None and r.isi_s.size > 0]) if runs else np.array([], dtype=float)
        fbdelay_all = np.concatenate([r.fb_delay_s for r in runs if r.fb_delay_s is not None and r.fb_delay_s.size > 0]) if runs else np.array([], dtype=float)

        med_isi = _nanmedian(isi_all)
        med_fb = _nanmedian(fbdelay_all)

        ant_vlines: List[Tuple[float, str]] = [
            (-CUE_DUR_S, "Cue onset"),
            (0.0, "Cue offset"),
            (med_isi, "Target (median)"),
            (med_fb, "Feedback (median)"),
            (6.0, "+6s"),
        ]
        ant_spans = [(0.0, 2.0, "Wu ANT window (0–2s)")]

        fb_vlines: List[Tuple[float, str]] = [
            (0.0, "Feedback onset"),
            (6.0, "+6s"),
        ]

        # Weighted aggregation across runs: weights = number of events contributing to each mean curve
        ant_psc_R = weighted_mean_and_sem([(r.ant_psc["Reward"][0], r.ant_psc["Reward"][2]) for r in runs])
        ant_psc_N = weighted_mean_and_sem([(r.ant_psc["Neutral"][0], r.ant_psc["Neutral"][2]) for r in runs])
        ant_z_R = weighted_mean_and_sem([(r.ant_z["Reward"][0], r.ant_z["Reward"][2]) for r in runs])
        ant_z_N = weighted_mean_and_sem([(r.ant_z["Neutral"][0], r.ant_z["Neutral"][2]) for r in runs])

        fb_psc_P = weighted_mean_and_sem([(r.fb_psc["Positive"][0], r.fb_psc["Positive"][2]) for r in runs])
        fb_psc_N = weighted_mean_and_sem([(r.fb_psc["Negative"][0], r.fb_psc["Negative"][2]) for r in runs])
        fb_z_P = weighted_mean_and_sem([(r.fb_z["Positive"][0], r.fb_z["Positive"][2]) for r in runs])
        fb_z_N = weighted_mean_and_sem([(r.fb_z["Negative"][0], r.fb_z["Negative"][2]) for r in runs])

        out_dir = OUT_TC_DIR / "subjects" / f"sub-{subject}" / echo
        plot_two_conditions(
            time_axis,
            ant_psc_R, "Reward",
            ant_psc_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=out_dir / "anticipation_psc.png",
            vlines=ant_vlines,
            spans=ant_spans,
        )
        plot_two_conditions(
            time_axis,
            ant_z_R, "Reward",
            ant_z_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, Z)",
            ylabel="Z (SD units)",
            out_png=out_dir / "anticipation_z.png",
            vlines=ant_vlines,
            spans=ant_spans,
        )

        plot_two_conditions(
            time_axis,
            fb_psc_P, "Feedback +",
            fb_psc_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=out_dir / "feedback_psc.png",
            vlines=fb_vlines,
        )
        plot_two_conditions(
            time_axis,
            fb_z_P, "Feedback +",
            fb_z_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, Z)",
            ylabel="Z (SD units)",
            out_png=out_dir / "feedback_z.png",
            vlines=fb_vlines,
        )

def _fmt(x: float | None) -> str:
    if x is None:
        return ""
    if np.isnan(x):
        return ""
    return f"{x:.6f}"


def main() -> None:
    if not FEAT_LIST_PATH.exists():
        print(f"[WARN] FEAT list not found: {FEAT_LIST_PATH}. Nothing to do.")
        return

    feat_paths: List[Path] = []
    for line in FEAT_LIST_PATH.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        feat_paths.append(Path(line))

    if not feat_paths:
        print("[WARN] No FEAT paths found in list. Nothing to do.")
        return

    OUT_TC_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    # Shared time axis (seconds)
    time_axis = build_time_axis(TMIN, TMAX, DT)

    # Process all FEATs
    results: List[RunResult] = []
    for feat in feat_paths:
        res = process_one_feat(feat, time_axis)
        if res is None:
            continue
        save_run_plots(res)
        results.append(res)

    if not results:
        print("[WARN] No runs processed successfully. Nothing to summarize.")
        return

    # Aggregate by subject (split by echo)
    subjects = sorted({r.sub for r in results})
    for sub in subjects:
        aggregate_subject(results, sub, time_axis)

    # --------------------------- Write summary table ---------------------------
    # Dynamic header based on requested offsets
    header = ["sub", "ses", "run", "echo", "tr_header", "tr_design", "tr_used", "med_isi", "med_fbdelay", "n_isi", "n_fbdelay"]
    # Put PSC columns first, then Z, then Ns
    keys_psc = sorted({k for r in results for k in r.points_psc.keys()})
    keys_z   = sorted({k for r in results for k in r.points_z.keys()})
    keys_n   = sorted({k for r in results for k in r.points_n.keys()})

    header += keys_psc + keys_z + keys_n

    out_tsv = SUMMARY_DIR / "summary_at_times_mid-unsmoothed.tsv"
    with open(out_tsv, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in results:
            med_isi = float(np.nanmedian(r.isi_s)) if r.isi_s is not None and r.isi_s.size > 0 else float('nan')
            med_fb = float(np.nanmedian(r.fb_delay_s)) if r.fb_delay_s is not None and r.fb_delay_s.size > 0 else float('nan')
            row = [
                r.sub, r.ses, r.run, r.echo,
                _fmt(r.tr_header), _fmt(r.tr_design), _fmt(r.tr_used),
                _fmt(med_isi), _fmt(med_fb),
                str(int(r.isi_s.size if r.isi_s is not None else 0)),
                str(int(r.fb_delay_s.size if r.fb_delay_s is not None else 0)),
            ]
            row += [_fmt(r.points_psc.get(k)) for k in keys_psc]
            row += [_fmt(r.points_z.get(k)) for k in keys_z]
            row += [str(r.points_n.get(k, 0)) for k in keys_n]
            f.write("\t".join(row) + "\n")


    # Quick TR summary (sanity check)
    trs = np.array([r.tr_header for r in results], dtype=float)
    if trs.size > 0:
        vals, counts = np.unique(np.round(trs, 6), return_counts=True)
        pairs = ", ".join([f"{v:.6f} ({c})" for v, c in zip(vals, counts)])
        print(f"TR (header) values across runs: {pairs}")
        # Also report any large header-vs-design mismatches
        mism = [r for r in results if r.tr_design is not None and abs(r.tr_design - r.tr_header) > 1e-3]
        if mism:
            print(f"[WARN] {len(mism)} run(s) have TR mismatch > 1e-3 between header and design.fsf (see summary TSV).")

    print(f"Done.\n  Plots: {OUT_TC_DIR}\n  Time-based summaries: {out_tsv}")


if __name__ == "__main__":
    main()
