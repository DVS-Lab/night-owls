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
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-unsmoothed_interpolated-time"
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


def load_ev_onsets(ev_path: Path) -> np.ndarray:
    """Load FSL EV (3 cols) and return onsets. Empty file -> empty array."""
    if not ev_path.exists():
        raise FileNotFoundError(f"Missing EV: {ev_path}")
    arr = np.loadtxt(ev_path, ndmin=2)
    if arr.size == 0:
        return np.array([], dtype=float)
    return np.asarray(arr[:, 0], dtype=float)


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

# ------------------------------- Core pipeline ----------------------------------

def process_one_feat(feat: Path, time_axis: np.ndarray) -> RunResult | None:
    try:
        sub, ses, run = parse_sub_ses_run_from_feat(feat)
        echo = echo_from_feat(feat)
        ts_raw, tr = load_vs_timeseries(feat)
    except Exception as e:
        print(f"[WARN] Skipping FEAT due to error: {feat}\n  -> {e}")
        return None

    ts_psc = to_psc(ts_raw)
    ts_z   = to_z(ts_raw)

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)
    ant_R = load_ev_onsets(ev_dir / "_anticipation_reward.txt")
    ant_N = load_ev_onsets(ev_dir / "_anticipation_neutral.txt")

    # Wu-style ANT anchoring (delay onset): shift cue-onset to cue-offset
    ant_R = ant_R + CUE_DUR_S
    ant_N = ant_N + CUE_DUR_S

    # Feedback (pooled by valence for plotting)
    fb_pos = np.sort(np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_positive_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_positive_neutral.txt"),
    ]))
    fb_neg = np.sort(np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_negative_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_negative_neutral.txt"),
    ]))

    # Feedback (separate conditions for summary table)
    fb_PR = load_ev_onsets(ev_dir / "_feedback_positive_reward.txt")
    fb_NR = load_ev_onsets(ev_dir / "_feedback_negative_reward.txt")
    fb_PN = load_ev_onsets(ev_dir / "_feedback_positive_neutral.txt")
    fb_NN = load_ev_onsets(ev_dir / "_feedback_negative_neutral.txt")

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
        ant_psc=ant_psc, ant_z=ant_z,
        fb_psc=fb_psc, fb_z=fb_z,
        points_psc=points_psc,
        points_z=points_z,
        points_n=points_n,
    )


def save_run_plots(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

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


def aggregate_subject(results: List[RunResult], subject: str, time_axis: np.ndarray) -> None:
    pool: Dict[str, Dict[str, Dict[str, List[Tuple[np.ndarray, int]]]]] = {}

    for res in results:
        if res.sub != subject:
            continue
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

    subj_out = OUT_TC_DIR / "subjects" / f"sub-{subject}"

    def weighted_mean_and_sem(items: List[Tuple[np.ndarray, int]]) -> Tuple[np.ndarray, np.ndarray, int]:
        if not items:
            return np.array([]), np.array([]), 0
        means = np.stack([m for m, _ in items], axis=0)
        weights = np.asarray([n for _, n in items], dtype=float)
        if np.sum(weights) == 0:
            w = np.ones_like(weights) / weights.size
        else:
            w = weights / np.sum(weights)
        wmean = np.nansum(means * w[:, None], axis=0)
        sd = np.nanstd(means, axis=0, ddof=1)
        sem = sd / np.sqrt(max(1, means.shape[0]))
        return wmean, sem, int(np.sum(weights))

    for echo, families in pool.items():
        ant_psc_R = weighted_mean_and_sem(families["ANT_PSC"]["Reward"])
        ant_psc_N = weighted_mean_and_sem(families["ANT_PSC"]["Neutral"])
        plot_two_conditions(
            time_axis, ant_psc_R, "Reward", ant_psc_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "anticipation_psc.png",
            vlines=VERT_LINES_ANT,
        )

        ant_z_R = weighted_mean_and_sem(families["ANT_Z"]["Reward"])
        ant_z_N = weighted_mean_and_sem(families["ANT_Z"]["Neutral"])
        plot_two_conditions(
            time_axis, ant_z_R, "Reward", ant_z_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "anticipation_z.png",
            vlines=VERT_LINES_ANT,
        )

        fb_psc_P = weighted_mean_and_sem(families["FB_PSC"]["Positive"])
        fb_psc_N = weighted_mean_and_sem(families["FB_PSC"]["Negative"])
        plot_two_conditions(
            time_axis, fb_psc_P, "Feedback +", fb_psc_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "feedback_psc.png",
            vlines=VERT_LINES_FB,
        )

        fb_z_P = weighted_mean_and_sem(families["FB_Z"]["Positive"])
        fb_z_N = weighted_mean_and_sem(families["FB_Z"]["Negative"])
        plot_two_conditions(
            time_axis, fb_z_P, "Feedback +", fb_z_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "feedback_z.png",
            vlines=VERT_LINES_FB,
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
    header = ["sub", "ses", "run", "echo"]
    # Put PSC columns first, then Z, then Ns
    keys_psc = sorted({k for r in results for k in r.points_psc.keys()})
    keys_z   = sorted({k for r in results for k in r.points_z.keys()})
    keys_n   = sorted({k for r in results for k in r.points_n.keys()})

    header += keys_psc + keys_z + keys_n

    out_tsv = SUMMARY_DIR / "summary_at_times_mid-unsmoothed.tsv"
    with open(out_tsv, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in results:
            row = [r.sub, r.ses, r.run, r.echo]
            row += [_fmt(r.points_psc.get(k)) for k in keys_psc]
            row += [_fmt(r.points_z.get(k)) for k in keys_z]
            row += [str(r.points_n.get(k, 0)) for k in keys_n]
            f.write("\t".join(row) + "\n")

    print(f"Done.\n  Plots: {OUT_TC_DIR}\n  Time-based summaries: {out_tsv}")


if __name__ == "__main__":
    main()
