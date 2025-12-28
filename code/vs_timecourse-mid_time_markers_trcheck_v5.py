#!/usr/bin/env python3
"""
VS time-course extraction for MID with timing markers + TR sanity checks (v5.2)

Fixes:
- Defines weighted_mean_and_sem (was missing in your failing copy).
- Loads per-condition feedback EVs (POS/NEG × REW/NEU) so point estimates work.
- Robust to empty conditions/runs (won't crash; will skip plots when needed).

What this script does
- ANT plots are anchored at *cue offset* (delay onset): ANT_onset = cue_onset + 0.75s.
- Adds vertical markers:
    cue onset  = -0.75s
    cue offset = 0.00s (the anchor)
    target (median)  = median(ANT_duration - 0.75)
    feedback (median)= median(next feedback onset - cue offset) using nearest subsequent feedback
    +6s marker
- Shades Wu-analogous anticipation window: 0–2s after cue offset.
- Checks TR consistency:
    TR from NIfTI header vs TR in design.fsf (fmri(tr)); warns if mismatched.
    Optionally prefer design.fsf TR for sampling (PREFER_TR_FROM_DESIGN).

Inputs expected (same conventions as your earlier scripts)
- FEAT paths listed in: code/feat_paths-unsmoothed.txt (relative to this script)
  If absent, falls back to glob under derivatives/fsl/sub-*/ses-*/L1_*.feat
- EVs in: derivatives/fsl/EVFiles/sub-<SUB>/ses-<SES>/mid/run-<RUN>/
    _anticipation_reward.txt
    _anticipation_neutral.txt
    _feedback_positive_reward.txt
    _feedback_negative_reward.txt
    _feedback_positive_neutral.txt
    _feedback_negative_neutral.txt
- VS mask in project masks directory (MNI space):
    masks/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz

Outputs
- derivatives/extractions/timecourses-mid-unsmoothed_markers/
    subjects/... (subject-level aggregated)
- derivatives/extractions/summary_mid_markers.tsv
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt


# --------------------------- Paths (project-relative) ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
EV_BASE    = FSL_DERIV / "EVFiles"
MASKS_DIR  = ROOT_DIR / "masks"

VS_MNI = MASKS_DIR / "space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-unsmoothed.txt"

OUT_DIR   = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-unsmoothed_markers"
SUMMARY_TSV = ROOT_DIR / "derivatives" / "extractions" / "summary_mid_markers.tsv"


# ----------------------------- Analysis parameters -----------------------------
TASK = "mid"

# Peri-event window (seconds)
TMIN = -4.0
TMAX = 16.0

# Sampling grid (seconds)
DT = 0.10  # 100 ms

# Cue duration
CUE_DUR_S = 0.75

# Baseline window in peri-event coordinates (seconds)
BASELINE_WIN = (-4.0, -2.0)

# Timepoints to sample for “point estimates” (seconds after the anchor)
POINT_SAMPLES_S = [0.0, 1.5, 3.0, 6.0]

# Markers / windows
WU_ANT_WINDOW = (0.0, 2.0)  # Wu-style 2s anticipation window after cue offset

# Plot toggles
PLOT_RUN_LEVEL = False  # subject-level is usually enough

# TR preference
PREFER_TR_FROM_DESIGN = False  # if True, use design.fsf TR for sampling when available

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

    tr_header: float
    tr_design: Optional[float]
    tr_used: float

    time_axis: np.ndarray

    # Timing diagnostics (seconds; relative to cue offset anchor)
    isi_s: np.ndarray          # target onset relative to cue offset: ANT_duration - 0.75
    fb_delay_s: np.ndarray     # next feedback onset - cue_offset

    ant_psc: Dict[str, Tuple[np.ndarray, np.ndarray, int]]  # Reward/Neutral
    ant_z:   Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_psc:  Dict[str, Tuple[np.ndarray, np.ndarray, int]]  # Positive/Negative (pooled)
    fb_z:    Dict[str, Tuple[np.ndarray, np.ndarray, int]]

    # Point estimates for TSV (includes the 4 feedback EVs separately)
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
    if "tedana" in s or "multi-echo" in s:
        return "multi-echo"
    if "fmriprep" in s or "single-echo" in s:
        return "single-echo"
    return "unknown-echo"


def get_ev_dir(sub: str, ses: str, run: str) -> Path:
    return EV_BASE / f"sub-{sub}" / f"ses-{ses}" / TASK / f"run-{run}"


def load_ev(ev_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load FSL EV 3-column file: onset duration amplitude. Returns (onsets, durations)."""
    if not ev_path.exists():
        return np.array([], dtype=float), np.array([], dtype=float)
    try:
        arr = np.loadtxt(ev_path, ndmin=2)
    except Exception:
        return np.array([], dtype=float), np.array([], dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    return np.asarray(arr[:, 0], dtype=float), np.asarray(arr[:, 1], dtype=float)


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
    mu = float(np.nanmean(ts))
    if mu == 0 or not np.isfinite(mu):
        return np.full_like(ts, np.nan, dtype=float)
    return ((ts / mu) - 1.0) * 100.0


def to_z(ts: np.ndarray) -> np.ndarray:
    mu = float(np.nanmean(ts))
    sd = float(np.nanstd(ts, ddof=1))
    if sd == 0 or not np.isfinite(sd):
        return np.full_like(ts, np.nan, dtype=float)
    return (ts - mu) / sd


def build_time_axis(tmin: float, tmax: float, dt: float) -> np.ndarray:
    n = int(np.floor((tmax - tmin) / dt)) + 1
    return (tmin + np.arange(n) * dt).astype(float)


def sample_windows(ts: np.ndarray, onsets_s: np.ndarray, tr: float, t_axis: np.ndarray) -> np.ndarray:
    """Return (n_events, n_time) sampled by linear interpolation."""
    if onsets_s.size == 0:
        return np.empty((0, t_axis.size), dtype=float)
    t_series = np.arange(ts.size) * tr
    W = np.empty((onsets_s.size, t_axis.size), dtype=float)
    W[:] = np.nan
    for i, onset in enumerate(onsets_s):
        t_abs = onset + t_axis
        W[i, :] = np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)
    return W


def baseline_subtract(W: np.ndarray, t_axis: np.ndarray, base_win: Tuple[float, float]) -> np.ndarray:
    if W.shape[0] == 0:
        return W
    b0, b1 = base_win
    bmask = (t_axis >= b0) & (t_axis <= b1)
    if not np.any(bmask):
        return W
    base = np.nanmean(W[:, bmask], axis=1, keepdims=True)
    return W - base


def mean_and_sem(W: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """Always returns correct-length vectors; SEM safe for small n."""
    n_time = W.shape[1]
    if W.shape[0] == 0:
        return np.full(n_time, np.nan), np.full(n_time, np.nan), 0

    mean = np.nanmean(W, axis=0)
    n = np.sum(np.isfinite(W), axis=0)

    sd = np.zeros(n_time, dtype=float)
    for j in range(n_time):
        if n[j] >= 2:
            sd[j] = float(np.nanstd(W[:, j], ddof=1))
        else:
            sd[j] = 0.0
    sem = sd / np.sqrt(np.maximum(n, 1))
    n_events = int(np.sum(~np.all(np.isnan(W), axis=1)))
    return mean, sem, n_events


def weighted_mean_and_sem(curves: List[Tuple[np.ndarray, np.ndarray, int]]) -> Tuple[np.ndarray, np.ndarray, int]:
    """Combine run-level curves into subject-level (weights = n_events)."""
    filtered = [(m, s, n) for (m, s, n) in curves if n and n > 0 and np.any(np.isfinite(m))]
    if not filtered:
        raise ValueError("No valid curves to combine (all n=0).")

    means = np.stack([m for m, _, _ in filtered], axis=0)
    weights = np.array([n for _, _, n in filtered], dtype=float)
    total_n = int(np.sum(weights))

    w = weights / np.sum(weights)
    wmean = np.sum(means * w[:, None], axis=0)

    if means.shape[0] == 1:
        wsem = np.zeros_like(wmean)
    else:
        wvar = np.sum(w[:, None] * (means - wmean[None, :]) ** 2, axis=0)
        wstd = np.sqrt(wvar)
        wsem = wstd / np.sqrt(means.shape[0])

    return wmean, wsem, total_n


def vals_at_offset(ts: np.ndarray, onsets_s: np.ndarray, tr: float, offset_s: float) -> np.ndarray:
    if onsets_s.size == 0:
        return np.array([], dtype=float)
    t_series = np.arange(ts.size) * tr
    t_abs = onsets_s + float(offset_s)
    return np.interp(t_abs, t_series, ts, left=np.nan, right=np.nan)


def nanmedian_or_nan(x: np.ndarray) -> float:
    if x is None or x.size == 0:
        return float("nan")
    return float(np.nanmedian(x))


def plot_two_conditions(time_axis: np.ndarray,
                        condA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
                        condB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
                        title: str, ylabel: str, out_png: Path,
                        vlines: List[Tuple[float, str]] | None = None,
                        spans: List[Tuple[float, float, str]] | None = None) -> None:
    mA, sA, nA = condA
    mB, sB, nB = condB

    fig, ax = plt.subplots(figsize=(8, 4.8))

    if spans:
        for x0, x1, lab in spans:
            ax.axvspan(x0, x1, alpha=0.10, label=lab if lab else None)

    ax.plot(time_axis, mA, label=f"{labelA} (n={nA})")
    ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)

    ax.plot(time_axis, mB, label=f"{labelB} (n={nB})")
    ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)

    if vlines:
        for x, _ in vlines:
            if np.isfinite(x):
                ax.axvline(x, ls=":", lw=1)

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


def _fmt(x: float | None) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return f"{x:.6f}"


# ------------------------------- Core pipeline ----------------------------------
def process_one_feat(feat: Path, time_axis: np.ndarray) -> Optional[RunResult]:
    try:
        sub, ses, run = parse_sub_ses_run_from_feat(feat)
        echo = echo_from_feat(feat)
        ts_raw, tr_header = load_vs_timeseries(feat)
        tr_design = load_tr_from_design_fsf(feat)
        tr_used = tr_header
        if tr_design is not None and abs(tr_design - tr_header) > 1e-3:
            print(f"[WARN] TR mismatch: {feat}\n       header={tr_header:.6f} design.fsf={tr_design:.6f}")
        if tr_design is not None and PREFER_TR_FROM_DESIGN:
            tr_used = float(tr_design)
    except Exception as e:
        print(f"[WARN] Skipping FEAT due to error: {feat}\n  -> {e}")
        return None

    ts_psc = to_psc(ts_raw)
    ts_z   = to_z(ts_raw)

    ev_dir = get_ev_dir(sub, ses, run)

    # Anticipation cue onsets and durations (cue onset)
    ant_R_on, ant_R_dur = load_ev(ev_dir / "_anticipation_reward.txt")
    ant_N_on, ant_N_dur = load_ev(ev_dir / "_anticipation_neutral.txt")

    # Anchor ANT at cue offset (Wu-like delay onset)
    ant_R = ant_R_on + CUE_DUR_S
    ant_N = ant_N_on + CUE_DUR_S

    # ISI proxy: target onset relative to cue offset
    isi_s = np.concatenate([ant_R_dur, ant_N_dur]) - CUE_DUR_S
    isi_s = isi_s[np.isfinite(isi_s)]

    # Feedback EVs (per condition)
    fb_PR_on, _ = load_ev(ev_dir / "_feedback_positive_reward.txt")
    fb_NR_on, _ = load_ev(ev_dir / "_feedback_negative_reward.txt")
    fb_PN_on, _ = load_ev(ev_dir / "_feedback_positive_neutral.txt")
    fb_NN_on, _ = load_ev(ev_dir / "_feedback_negative_neutral.txt")

    fb_pos = np.sort(np.concatenate([fb_PR_on, fb_PN_on]))
    fb_neg = np.sort(np.concatenate([fb_NR_on, fb_NN_on]))
    fb_all = np.sort(np.concatenate([fb_pos, fb_neg]))

    # Feedback delay relative to each cue offset: next subsequent feedback onset
    cueoffs = np.sort(np.concatenate([ant_R, ant_N]))
    if cueoffs.size == 0 or fb_all.size == 0:
        fb_delay_s = np.array([], dtype=float)
    else:
        idx_next = np.searchsorted(fb_all, cueoffs, side="right")
        valid = idx_next < fb_all.size
        fb_delay_s = fb_all[idx_next[valid]] - cueoffs[valid]

    # Windows + baseline subtract
    R_psc_w = baseline_subtract(sample_windows(ts_psc, ant_R, tr_used, time_axis), time_axis, BASELINE_WIN)
    N_psc_w = baseline_subtract(sample_windows(ts_psc, ant_N, tr_used, time_axis), time_axis, BASELINE_WIN)

    P_psc_w = baseline_subtract(sample_windows(ts_psc, fb_pos, tr_used, time_axis), time_axis, BASELINE_WIN)
    G_psc_w = baseline_subtract(sample_windows(ts_psc, fb_neg, tr_used, time_axis), time_axis, BASELINE_WIN)

    R_z_w = baseline_subtract(sample_windows(ts_z, ant_R, tr_used, time_axis), time_axis, BASELINE_WIN)
    N_z_w = baseline_subtract(sample_windows(ts_z, ant_N, tr_used, time_axis), time_axis, BASELINE_WIN)

    P_z_w = baseline_subtract(sample_windows(ts_z, fb_pos, tr_used, time_axis), time_axis, BASELINE_WIN)
    G_z_w = baseline_subtract(sample_windows(ts_z, fb_neg, tr_used, time_axis), time_axis, BASELINE_WIN)

    ant_psc = {"Reward": mean_and_sem(R_psc_w), "Neutral": mean_and_sem(N_psc_w)}
    ant_z   = {"Reward": mean_and_sem(R_z_w),   "Neutral": mean_and_sem(N_z_w)}
    fb_psc  = {"Positive": mean_and_sem(P_psc_w), "Negative": mean_and_sem(G_psc_w)}
    fb_z    = {"Positive": mean_and_sem(P_z_w),   "Negative": mean_and_sem(G_z_w)}

    # Point estimates (sample raw ts, not baseline-subtracted curves)
    points_psc: Dict[str, float | None] = {}
    points_z: Dict[str, float | None] = {}
    points_n: Dict[str, int] = {}

    def m(arr: np.ndarray) -> float | None:
        if arr.size == 0 or not np.any(np.isfinite(arr)):
            return None
        return float(np.nanmean(arr))

    def n_valid(arr: np.ndarray) -> int:
        return int(np.sum(np.isfinite(arr))) if arr.size else 0

    for off in POINT_SAMPLES_S:
        tag = f"{off:.3f}".rstrip("0").rstrip(".")
        # ANT (reward/neutral)
        v = vals_at_offset(ts_psc, ant_R, tr_used, off)
        points_psc[f"ANT_REWARD_t{tag}"] = m(v); points_n[f"N_ANT_REWARD_t{tag}"] = n_valid(v)
        v = vals_at_offset(ts_psc, ant_N, tr_used, off)
        points_psc[f"ANT_NEUTRAL_t{tag}"] = m(v); points_n[f"N_ANT_NEUTRAL_t{tag}"] = n_valid(v)

        v = vals_at_offset(ts_z, ant_R, tr_used, off)
        points_z[f"ANT_REWARD_t{tag}"] = m(v)
        v = vals_at_offset(ts_z, ant_N, tr_used, off)
        points_z[f"ANT_NEUTRAL_t{tag}"] = m(v)

        # FB per EV (pos/neg × rew/neu)
        for name, on in [
            ("FB_POS_REWARD", fb_PR_on),
            ("FB_NEG_REWARD", fb_NR_on),
            ("FB_POS_NEUTRAL", fb_PN_on),
            ("FB_NEG_NEUTRAL", fb_NN_on),
        ]:
            v = vals_at_offset(ts_psc, on, tr_used, off)
            points_psc[f"{name}_t{tag}"] = m(v); points_n[f"N_{name}_t{tag}"] = n_valid(v)
            vz = vals_at_offset(ts_z, on, tr_used, off)
            points_z[f"{name}_t{tag}"] = m(vz)

        # FB pooled
        v = vals_at_offset(ts_psc, fb_pos, tr_used, off)
        points_psc[f"FB_POS_t{tag}"] = m(v); points_n[f"N_FB_POS_t{tag}"] = n_valid(v)
        v = vals_at_offset(ts_psc, fb_neg, tr_used, off)
        points_psc[f"FB_NEG_t{tag}"] = m(v); points_n[f"N_FB_NEG_t{tag}"] = n_valid(v)

        vz = vals_at_offset(ts_z, fb_pos, tr_used, off)
        points_z[f"FB_POS_t{tag}"] = m(vz)
        vz = vals_at_offset(ts_z, fb_neg, tr_used, off)
        points_z[f"FB_NEG_t{tag}"] = m(vz)

    return RunResult(
        feat_path=feat, sub=sub, ses=ses, run=run, echo=echo,
        tr_header=tr_header, tr_design=tr_design, tr_used=tr_used,
        time_axis=time_axis,
        isi_s=isi_s,
        fb_delay_s=fb_delay_s,
        ant_psc=ant_psc, ant_z=ant_z,
        fb_psc=fb_psc, fb_z=fb_z,
        points_psc=points_psc,
        points_z=points_z,
        points_n=points_n,
    )


def markers_for_ant(isi_s: np.ndarray, fb_delay_s: np.ndarray) -> Tuple[List[Tuple[float, str]], List[Tuple[float, float, str]]]:
    med_isi = nanmedian_or_nan(isi_s)
    med_fb = nanmedian_or_nan(fb_delay_s)
    vlines = [
        (-CUE_DUR_S, "Cue onset"),
        (0.0, "Cue offset"),
        (med_isi, "Target (median)"),
        (med_fb, "Feedback (median)"),
        (6.0, "+6s"),
    ]
    spans = [(WU_ANT_WINDOW[0], WU_ANT_WINDOW[1], "Wu ANT window (0–2s)")]
    return vlines, spans


def list_feat_paths() -> List[Path]:
    if FEAT_LIST_PATH.exists():
        feats: List[Path] = []
        for line in FEAT_LIST_PATH.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            feats.append(Path(line))
        return feats
    return sorted(FSL_DERIV.glob("sub-*/ses-*/L1_*.feat"))


def aggregate_subject(results: List[RunResult], subject: str, time_axis: np.ndarray) -> None:
    subj_runs = [r for r in results if r.sub == subject]
    if not subj_runs:
        return

    by_echo: Dict[str, List[RunResult]] = {}
    for r in subj_runs:
        by_echo.setdefault(r.echo, []).append(r)

    for echo, runs in by_echo.items():
        isi_all = np.concatenate([r.isi_s for r in runs if r.isi_s.size > 0]) if runs else np.array([], dtype=float)
        fb_all = np.concatenate([r.fb_delay_s for r in runs if r.fb_delay_s.size > 0]) if runs else np.array([], dtype=float)
        vlines, spans = markers_for_ant(isi_all, fb_all)

        try:
            antR = weighted_mean_and_sem([r.ant_psc["Reward"] for r in runs])
            antN = weighted_mean_and_sem([r.ant_psc["Neutral"] for r in runs])
            fbP = weighted_mean_and_sem([r.fb_psc["Positive"] for r in runs])
            fbG = weighted_mean_and_sem([r.fb_psc["Negative"] for r in runs])
        except ValueError:
            continue

        out_dir = OUT_DIR / "subjects" / f"sub-{subject}" / echo
        plot_two_conditions(
            time_axis,
            antR, "Reward",
            antN, "Neutral",
            title=f"VS — subject {subject} [{echo}] (ANT, PSC)",
            ylabel="% signal change (PSC)",
            out_png=out_dir / "anticipation_psc.png",
            vlines=vlines,
            spans=spans,
        )
        plot_two_conditions(
            time_axis,
            fbP, "Feedback +",
            fbG, "Feedback −",
            title=f"VS — subject {subject} [{echo}] (FB, PSC)",
            ylabel="% signal change (PSC)",
            out_png=out_dir / "feedback_psc.png",
            vlines=[(0.0, "Feedback onset"), (6.0, "+6s")],
            spans=None,
        )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_TSV.parent.mkdir(parents=True, exist_ok=True)

    time_axis = build_time_axis(TMIN, TMAX, DT)

    feats = list_feat_paths()
    if not feats:
        print("[ERROR] No FEATs found.")
        return

    results: List[RunResult] = []
    tr_vals = []
    tr_mismatch = 0

    for feat in feats:
        res = process_one_feat(feat, time_axis)
        if res is None:
            continue
        results.append(res)
        tr_vals.append(res.tr_used)
        if res.tr_design is not None and abs(res.tr_design - res.tr_header) > 1e-3:
            tr_mismatch += 1

    if not results:
        print("[ERROR] No runs processed successfully.")
        return

    # Subject-level plots
    for sub in sorted({r.sub for r in results}):
        aggregate_subject(results, sub, time_axis)

    # TR summary
    if tr_vals:
        vals, counts = np.unique(np.round(np.array(tr_vals, dtype=float), 6), return_counts=True)
        s = ", ".join([f"{v:.6f} ({c})" for v, c in zip(vals, counts)])
        print(f"TR used distribution: {s}")
    if tr_mismatch:
        print(f"[WARN] TR mismatch (design.fsf vs NIfTI header) in {tr_mismatch} run(s).")

    # Write summary TSV (run-level)
    base_cols = ["sub", "ses", "run", "echo", "tr_header", "tr_design", "tr_used",
                 "med_isi", "med_fbdelay", "n_isi", "n_fbdelay"]
    keys_psc = sorted({k for r in results for k in r.points_psc.keys()})
    keys_z = sorted({k for r in results for k in r.points_z.keys()})
    keys_n = sorted({k for r in results for k in r.points_n.keys()})
    header = base_cols + keys_psc + keys_z + keys_n

    with open(SUMMARY_TSV, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in results:
            med_isi = nanmedian_or_nan(r.isi_s)
            med_fb = nanmedian_or_nan(r.fb_delay_s)
            row = [
                r.sub, r.ses, r.run, r.echo,
                f"{r.tr_header:.6f}",
                "" if r.tr_design is None else f"{r.tr_design:.6f}",
                f"{r.tr_used:.6f}",
                "" if np.isnan(med_isi) else f"{med_isi:.6f}",
                "" if np.isnan(med_fb) else f"{med_fb:.6f}",
                str(int(r.isi_s.size)),
                str(int(r.fb_delay_s.size)),
            ]
            row += [_fmt(r.points_psc.get(k)) for k in keys_psc]
            row += [_fmt(r.points_z.get(k)) for k in keys_z]
            row += [str(r.points_n.get(k, 0)) for k in keys_n]
            f.write("\t".join(row) + "\n")

    print(f"Done.\n  Plots: {OUT_DIR}\n  Summary: {SUMMARY_TSV}")


if __name__ == "__main__":
    main()
