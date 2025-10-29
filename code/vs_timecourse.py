#!/usr/bin/env python3
"""
VS time‑course extraction and plotting (anticipation + feedback), echo‑split
-----------------------------------------------------------------------------
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
  (so fractional onsets are handled) for two analysis families:
    • ANTICIPATION: Reward vs Neutral (cue‑locked)
    • FEEDBACK (valence): Positive vs Negative (feedback‑locked, pooling across
      incentive conditions)
- Produces per‑run plots and subject‑level aggregated plots split by echo
  (single‑echo vs multi‑echo) for PSC and Z (four figures per subject per family).
- Saves simple text summaries including the number of trials included.

Outputs
  <rootdir>/derivatives/extractions/timecourses/
    ├── runs/  (per FEAT)
    └── subjects/  (per subject, echo‑split)

Dependencies: nibabel, numpy, matplotlib
"""
from __future__ import annotations

import os
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
OUT_DIR    = ROOT_DIR / "derivatives" / "extractions" / "timecourses"
EV_BASE    = FSL_DERIV / "EVFiles"  # EVs are arranged by sub/ses/task/run below

VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

# Input list of FEAT directories (one per line, absolute). You provided this.
FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths.txt"

# --------------------------- Analysis parameters ------------------------------
TASK         = "mid"
TR_HARDCODE  = 1.615  # sec (we also check header and warn on mismatch)
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
    if "single-echo" in feat.name:
        return "single-echo"
    if "multi-echo" in feat.name:
        return "multi-echo"
    # Fallback to searching full path
    if "single-echo" in feat.as_posix():
        return "single-echo"
    if "multi-echo" in feat.as_posix():
        return "multi-echo"
    return "unknown-echo"


def get_ev_dir(sub: str, ses: str, run: str) -> Path:
    return EV_BASE / f"sub-{sub}" / f"ses-{ses}" / TASK / f"run-{run}"


def load_ev(ev_path: Path) -> np.ndarray:
    """Load a 3‑column EV (onset, duration, amplitude) -> onsets (seconds).
    The time‑course uses only onsets; durations are ignored for event‑locking.
    Returns a 1D numpy array of onsets in seconds. Missing file -> empty array.
    """
    if not ev_path.exists():
        return np.array([], dtype=float)
    try:
        arr = np.loadtxt(ev_path, ndmin=2)
    except Exception:
        return np.array([], dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float)
    # First column is onset in seconds
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
    mask_path = VS_MNI
    if not mask_path.exists():
        raise FileNotFoundError(f"Missing VS mask at {mask_path}")
    img = nib.load(str(img_path))
    mask_img = nib.load(str(mask_path))

    data = img.get_fdata()  # X×Y×Z×T
    mask = mask_img.get_fdata() > 0
    # Basic sanity checks on spatial dims
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
    """Return array of shape (n_trials, n_time) sampled by linear interpolation.
    Out‑of‑range samples are set to NaN and later dropped from stats.
    """
    if onsets_s.size == 0:
        return np.empty((0, 0))
    t_series = np.arange(ts.size) * tr
    t_axis = build_time_axis(tmin, tmax, tr)
    windows = []
    for onset in onsets_s:
        t_abs = onset + t_axis
        # Linear interpolation; out‑of‑bounds -> NaN
        vals = np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)
        windows.append(vals)
    return np.vstack(windows)


def mean_and_sem(windows: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    if windows.size == 0:
        return np.array([]), np.array([]), 0
    # Drop timepoints with all NaN across trials
    valid = ~np.all(~np.isfinite(windows), axis=0)
    w = np.where(valid, windows, np.nan)
    mean = np.nanmean(w, axis=0)
    n = np.sum(np.isfinite(w), axis=0)
    sem = np.nanstd(w, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))
    # For reporting, use the number of complete trials (with at least one finite)
    n_trials = int(np.sum(~np.all(np.isnan(w), axis=1)))
    return mean, sem, n_trials


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

    # Load EVs
    ev_dir = get_ev_dir(sub, ses, run)
    ant_R = load_ev(ev_dir / "_anticipation_reward.txt")
    ant_N = load_ev(ev_dir / "_anticipation_neutral.txt")
    # Feedback: pool across incentive (reward/neutral) within valence
    fb_pos = np.sort(np.concatenate([
        load_ev(ev_dir / "_feedback_positive_reward.txt"),
        load_ev(ev_dir / "_feedback_positive_neutral.txt"),
    ]))
    fb_neg = np.sort(np.concatenate([
        load_ev(ev_dir / "_feedback_negative_reward.txt"),
        load_ev(ev_dir / "_feedback_negative_neutral.txt"),
    ]))

    t_axis = build_time_axis(TMIN, TMAX, tr)

    # ANTICIPATION windows
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

    # FEEDBACK (valence) windows
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

    return RunResult(
        feat_path=feat,
        sub=sub, ses=ses, run=run, echo=echo,
        time_axis=t_axis,
        ant_psc=ant_psc, ant_z=ant_z,
        fb_psc=fb_psc, fb_z=fb_z,
    )


def save_run_plots(res: RunResult) -> None:
    run_out = OUT_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

    # Anticipation: PSC and Z
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
        ylabel="Z (SD units)",
        out_png=run_out / "anticipation_z.png",
    )

    # Feedback (valence): PSC and Z
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


# -------------------------- Subject‑level aggregation --------------------------

def aggregate_subject(results: List[RunResult], subject: str) -> None:
    # Organize windows by echo and family/condition/metric
    buckets: Dict[str, Dict[str, Dict[str, List[np.ndarray]]]] = {}
    # Structure: buckets[echo][family_metric_condition] -> list of (mean across trials?)
    # We keep per‑trial windows by reconstructing from mean/sem is impossible.
    # Instead, we will average across runs using the means, and combine SEM via
    # a simple pooled SEM approximation using trial counts. To do that, we need
    # to retain mean vectors and counts.

    # Instead of reconstructing trials, pool means with weights by n_trials.
    # We will store for each echo/family/metric/condition: list of (mean_vec, n)
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

    subj_out = OUT_DIR / "subjects" / f"sub-{subject}"

    def weighted_mean_and_sem(items: List[Tuple[np.ndarray, int]]) -> Tuple[np.ndarray, np.ndarray, int]:
        if len(items) == 0:
            return np.array([]), np.array([]), 0
        means = np.vstack([m for (m, n) in items])
        ns = np.array([max(n, 0) for (m, n) in items], dtype=float)
        # Avoid zero division
        w = ns / np.clip(ns.sum(), 1.0, None)
        mean = np.nansum(means * w[:, None], axis=0)
        # Approximate SEM across runs weighted by trial counts
        # First compute run‑wise variance estimates from SEM if we had them;
        # since we do not, approximate by variance of means across runs and then
        # divide by sqrt(number of runs). This is conservative.
        sem = np.nanstd(means, axis=0, ddof=1) / np.sqrt(max(len(items), 1))
        n_total = int(ns.sum())
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
        )
        # ANTICIPATION — Z
        ant_z_R = weighted_mean_and_sem(families["ANT_Z"]["Reward"])
        ant_z_N = weighted_mean_and_sem(families["ANT_Z"]["Neutral"])
        plot_two_conditions(
            time_axis, ant_z_R, "Reward", ant_z_N, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "anticipation_z.png",
        )
        # FEEDBACK (valence) — PSC
        fb_psc_P = weighted_mean_and_sem(families["FB_PSC"]["Positive"])
        fb_psc_N = weighted_mean_and_sem(families["FB_PSC"]["Negative"])
        plot_two_conditions(
            time_axis, fb_psc_P, "Feedback +", fb_psc_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=subj_out / echo / "feedback_psc.png",
        )
        # FEEDBACK (valence) — Z
        fb_z_P = weighted_mean_and_sem(families["FB_Z"]["Positive"])
        fb_z_N = weighted_mean_and_sem(families["FB_Z"]["Negative"])
        plot_two_conditions(
            time_axis, fb_z_P, "Feedback +", fb_z_N, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, Z)",
            ylabel="Z (SD units)",
            out_png=subj_out / echo / "feedback_z.png",
        )
        # Summary text
        summary = (
            f"Subject {subject} — Echo: {echo}\n"
            f"ANT  Reward n={ant_psc_R[2]}, Neutral n={ant_psc_N[2]}\n"
            f"FB   Positive n={fb_psc_P[2]}, Negative n={fb_psc_N[2]}\n"
        )
        (subj_out / echo / "summary.txt").write_text(summary)


# ---------------------------------- Driver ------------------------------------

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
    for feat in feat_paths:
        res = process_one_feat(feat)
        if res is None:
            continue
        save_run_plots(res)
        results.append(res)

    # Aggregate by subject (split by echo)
    subjects = sorted({r.sub for r in results})
    for sub in subjects:
        aggregate_subject(results, sub)

    print(f"Done. Outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
