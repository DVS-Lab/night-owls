#!/usr/bin/env python3
"""
Wu-normalized event-related ROI timecourses for Night-Owls (MID + SharedReward).

What this script does (baked-in; no CLI args):
- Reads a list of L1 FEAT directories from feat_paths_704.txt (or feat_paths_all.txt),
  located in the same folder as this script.
- Extracts an ROI mean time series from filtered_func_data.nii.gz for each FEAT.
- Regresses confounds from confoundevs.txt (if present).
- Converts the residualized ROI time series to Wu-style percent signal change (PSC)
  relative to the run mean: PSC(t) = ((x(t)/mean(x)) - 1) * 100.
- Epochs PSC around events (MID and SharedReward EVs) and produces per-run curves
  (mean +/- SEM across events).
- Produces session- and subject-level QC summaries by averaging run-level curves.

Outputs (under <project_root>/derivatives/extractions/timecourses-unified-wu_noargs_v2/):
- tables/run_curves_long_wu.csv.gz  (long-format per-run curves)
- run_qc/*.png                      (per-run curves)
- session_qc/*.png                  (run-averaged session summaries)
- subject_qc/*.png                  (session-averaged subject summaries)

Assumptions (matches your project layout):
- FEAT dirs: derivatives/fsl/sub-*/ses-*/L1_*.feat/
- EV files: derivatives/fsl/EVFiles/sub-*/ses-*/<task>/run-<run>/*.txt
- ROI mask: masks/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz (in project root)
"""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# ----------------------------- CONFIG (baked-in) -----------------------------

SCRIPT_DIR = Path(__file__).resolve().parent

# One of these must exist next to this script.
PATHLIST_CANDIDATES = [
    SCRIPT_DIR / "feat_paths_704.txt",
    SCRIPT_DIR / "feat_paths_all.txt",
]

MASK_REL = Path("masks/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz")

# Event-locked time axis
DT = 0.1
TMIN = -6.0
TMAX = 16.0

# Expected TR (used only to sanity-check; actual TR read from NIfTI header)
TR_EXPECTED = 1.615

# SharedReward: when computing "next event" marker, ignore gaps larger than this.
SR_MAX_DEC_TO_OUT_GAP_S = 8.0

# Plot styling
DPI = 150

# Robust y-limits so a single outlier run doesn't flatten everything.
# Quantiles computed across all plotted line points.
YLIM_QLO = 0.01
YLIM_QHI = 0.99
YLIM_PAD_FRAC = 0.10

PSC_METHOD = "wu_runmean"


# ----------------------------- helpers -----------------------------

FEAT_NAME_RE = re.compile(r".*task-(mid|sharedreward).*run-(\d+).*", re.IGNORECASE)


def find_pathlist_file() -> Path:
    for p in PATHLIST_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError(
        "No path list found. Expected one of:\n  - "
        + "\n  - ".join(str(p) for p in PATHLIST_CANDIDATES)
    )


def read_feat_paths(pathlist: Path) -> List[Path]:
    out: List[Path] = []
    for line in pathlist.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        out.append(Path(s))
    return out


def infer_project_root_from_feat(feat: Path) -> Path:
    """
    Tries to infer <project_root> from:
      <project_root>/derivatives/fsl/sub-*/ses-*/<feat>.feat
    """
    ps = feat.resolve().parents
    # expected: .../<root>/derivatives/fsl/sub-*/ses-*/<feat>.feat
    if len(ps) >= 6 and ps[2].name.startswith("sub-") and ps[1].name.startswith("ses-"):
        if ps[3].name == "fsl" and ps[4].name == "derivatives":
            return ps[5]
    # fallback: walk upward until we find "derivatives"
    cur = feat.resolve()
    for parent in cur.parents:
        if parent.name == "derivatives":
            return parent.parent
    # last resort: directory containing derivatives
    return feat.resolve().parents[3]


def load_mask(mask_path: Path) -> np.ndarray:
    img = nib.load(str(mask_path))
    dat = img.get_fdata()
    if dat.ndim != 3:
        raise ValueError(f"Expected 3D mask, got {dat.shape}: {mask_path}")
    return dat > 0


def parse_sub_ses_run_task(feat: Path) -> Tuple[str, str, str, str]:
    s = feat.name
    sub = re.search(r"sub-([^_]+)", s)
    ses = re.search(r"ses-([^_]+)", s)
    run = re.search(r"run-([^_]+)", s)
    task = re.search(r"task-(mid|sharedreward)", s)
    if not (sub and ses and run and task):
        # fall back to full path search
        sp = str(feat)
        sub = sub or re.search(r"/sub-([^/_]+)", sp)
        ses = ses or re.search(r"/ses-([^/_]+)", sp)
        run = run or re.search(r"run-([^_./]+)", sp)
        task = task or re.search(r"task-(mid|sharedreward)", sp)
    if not (sub and ses and run and task):
        raise ValueError(f"Could not parse sub/ses/run/task from FEAT name: {feat}")
    return sub.group(1), ses.group(1), run.group(1), task.group(1)


def infer_echo_and_cnfds(feat: Path) -> Tuple[str, str]:
    s = feat.name.lower()
    echo = "multi-echo" if "multi-echo" in s else ("single-echo" if "single-echo" in s else "unknown")
    m = re.search(r"cnfds-([^_]+)", feat.name)
    cnfds = m.group(1) if m else "unknown"
    return echo, cnfds


def regress_confounds(ts: np.ndarray, confound_file: Path) -> np.ndarray:
    """
    Regress columns of confoundevs.txt from the ROI mean time series.
    Returns residual + mean(ts) to preserve overall scaling (then PSC happens later).
    """
    if ts.size == 0:
        return ts
    if not confound_file.exists():
        return ts
    try:
        X = np.loadtxt(confound_file)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        if X.shape[0] != ts.shape[0]:
            return ts
        keep = np.nanstd(X, axis=0) > 0
        if not np.any(keep):
            return ts
        X = X[:, keep]
        X = np.column_stack([X, np.ones(X.shape[0])])
        beta = np.linalg.lstsq(X, ts, rcond=None)[0]
        resid = ts - (X @ beta)
        return resid + np.nanmean(ts)
    except Exception:
        return ts


def roi_mean_ts(func_nii: Path, mask_bool: np.ndarray) -> Tuple[np.ndarray, float]:
    img = nib.load(str(func_nii))
    data = img.get_fdata()
    if data.ndim != 4:
        raise ValueError(f"Expected 4D func data, got {data.shape} in {func_nii}")

    if mask_bool.shape != data.shape[:3]:
        raise ValueError(
            f"Mask grid mismatch:\n  func: {data.shape[:3]} {func_nii}\n  mask: {mask_bool.shape}"
        )

    ts = np.nanmean(data[mask_bool], axis=0)
    tr = float(img.header.get_zooms()[3])
    if abs(tr - TR_EXPECTED) > 0.2:
        # Not fatal—just warn
        print(f"[WARN] TR={tr:.4f}s differs from expected {TR_EXPECTED:.4f}s for {func_nii}")
    return ts.astype(float), tr


def to_wu_psc(ts: np.ndarray) -> np.ndarray:
    """
    Wu-style percent signal change relative to run mean.
    """
    mu = float(np.nanmean(ts))
    if not np.isfinite(mu) or abs(mu) < 1e-12:
        return np.full_like(ts, np.nan, dtype=float)
    return ((ts / mu) - 1.0) * 100.0


def epoch_interpolate(ts: np.ndarray, onsets: np.ndarray, tr: float,
                     tmin: float, tmax: float, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns:
      t_axis (len T)
      wins  (n_events x len T) in the same units as ts (here: PSC)
    """
    if onsets.size == 0:
        t_axis = np.arange(tmin, tmax + 1e-9, dt)
        return t_axis, np.zeros((0, t_axis.size), dtype=float)

    t_axis = np.arange(tmin, tmax + 1e-9, dt)
    t0 = np.arange(ts.size, dtype=float) * tr

    wins = np.full((onsets.size, t_axis.size), np.nan, dtype=float)
    for i, onset in enumerate(onsets.astype(float)):
        samp = onset + t_axis
        good = (samp >= t0[0]) & (samp <= t0[-1])
        if np.any(good):
            wins[i, good] = np.interp(samp[good], t0, ts)
    return t_axis, wins


def mean_sem_nanaware(wins: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    if wins.size == 0:
        return np.array([]), np.array([]), 0
    mu = np.nanmean(wins, axis=0)
    n_t = np.sum(np.isfinite(wins), axis=0)
    sd = np.nanstd(wins, axis=0, ddof=1)
    sem = sd / np.sqrt(np.maximum(n_t, 1))
    return mu, sem, int(wins.shape[0])


def robust_ylim_from_lines(ax) -> None:
    ys = []
    for line in ax.get_lines():
        y = line.get_ydata()
        y = y[np.isfinite(y)]
        if y.size:
            ys.append(y)
    if not ys:
        return
    y = np.concatenate(ys)
    lo = float(np.quantile(y, YLIM_QLO))
    hi = float(np.quantile(y, YLIM_QHI))
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return
    pad = (hi - lo) * YLIM_PAD_FRAC
    ax.set_ylim(lo - pad, hi + pad)


def add_event_markers(ax, event_dur: float, next_onset: float) -> None:
    """
    Adds vertical markers at:
      - 0 s (event onset)
      - event_dur (median duration, if finite)
      - next_onset (median onset to next phase, if finite)
    """
    yl = ax.get_ylim()
    ytop = yl[1]
    ax.axvline(0.0, linewidth=1.0, alpha=0.6)
    ax.text(0.0, ytop, " onset", va="top", ha="left", fontsize=8, alpha=0.8)

    if np.isfinite(event_dur):
        ax.axvline(event_dur, linewidth=1.0, alpha=0.4, linestyle="--")
        ax.text(event_dur, ytop, " end", va="top", ha="left", fontsize=8, alpha=0.8)

    if np.isfinite(next_onset):
        ax.axvline(next_onset, linewidth=1.0, alpha=0.4, linestyle=":")
        ax.text(next_onset, ytop, " next", va="top", ha="left", fontsize=8, alpha=0.8)


def plot_curves_with_sem(out_png: Path, title: str, t: np.ndarray,
                         curves: Dict[str, Tuple[np.ndarray, np.ndarray, int]],
                         event_dur: float, next_onset: float) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8.5, 5.0))
    ax = fig.add_subplot(111)

    for name in sorted(curves.keys()):
        mu, se, n = curves[name]
        if mu.size == 0:
            continue
        ax.plot(t, mu, label=f"{name} (n={n})")
        ax.fill_between(t, mu - se, mu + se, alpha=0.15)

    ax.set_title(title)
    ax.set_xlabel("Time from event onset (s)")
    ax.set_ylabel("PSC (% of run mean; confound-regressed)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, ncol=2)

    robust_ylim_from_lines(ax)
    add_event_markers(ax, event_dur, next_onset)

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)


def safe_slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9\-]+", "-", s).strip("-")


# ----------------------------- EV loading -----------------------------

def load_ev_3col(p: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    FSL 3-column EV: onset(s), duration(s), amplitude.
    We return (onsets, durations). Missing file => empty arrays.
    """
    if not p.exists():
        return np.array([]), np.array([])
    try:
        x = np.loadtxt(p)
        if x.size == 0:
            return np.array([]), np.array([])
        if x.ndim == 1:
            x = x.reshape(1, -1)
        on = x[:, 0].astype(float)
        dur = x[:, 1].astype(float)
        return on, dur
    except Exception:
        return np.array([]), np.array([])


def mid_ev_sets(ev_dir: Path) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """
    MID EV sets:
      - mid_anticipation_primary: Reward vs Neutral
      - mid_feedback_extra: Positive vs Negative (collapsed across cue types)
    """
    ant_r = load_ev_3col(ev_dir / "_anticipation_reward.txt")
    ant_n = load_ev_3col(ev_dir / "_anticipation_neutral.txt")

    fb_pos_r = load_ev_3col(ev_dir / "_feedback_positive_reward.txt")
    fb_neg_r = load_ev_3col(ev_dir / "_feedback_negative_reward.txt")
    fb_pos_n = load_ev_3col(ev_dir / "_feedback_positive_neutral.txt")
    fb_neg_n = load_ev_3col(ev_dir / "_feedback_negative_neutral.txt")

    out: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    out["mid_anticipation_primary"] = {"Reward": ant_r, "Neutral": ant_n}

    pos_on = np.concatenate([fb_pos_r[0], fb_pos_n[0]]) if (fb_pos_r[0].size or fb_pos_n[0].size) else np.array([])
    pos_d  = np.concatenate([fb_pos_r[1], fb_pos_n[1]]) if (fb_pos_r[1].size or fb_pos_n[1].size) else np.array([])
    neg_on = np.concatenate([fb_neg_r[0], fb_neg_n[0]]) if (fb_neg_r[0].size or fb_neg_n[0].size) else np.array([])
    neg_d  = np.concatenate([fb_neg_r[1], fb_neg_n[1]]) if (fb_neg_r[1].size or fb_neg_n[1].size) else np.array([])
    out["mid_feedback_extra"] = {"Positive": (pos_on, pos_d), "Negative": (neg_on, neg_d)}
    return out


def sr_ev_sets(ev_dir: Path) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """
    SharedReward EV sets:
      - sr_outcome_primary: Reward / Neutral / Punish (pooled across partner)
      - sr_decision_partner_extra: Computer / Stranger
      - sr_decision_by_future_outcome_extra: Reward / Neutral / Punish (decision locked)
    """
    out: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}

    # decision
    d_comp = load_ev_3col(ev_dir / "_guess_computer.txt")
    d_str  = load_ev_3col(ev_dir / "_guess_face.txt")
    out["sr_decision_partner_extra"] = {"Computer": d_comp, "Stranger": d_str}

    # outcomes pooled
    reward_on, reward_d = [], []
    neutral_on, neutral_d = [], []
    punish_on, punish_d = [], []
    onset_to_val: Dict[float, str] = {}

    for p in ev_dir.glob("_outcome_*.txt"):
        ons, durs = load_ev_3col(p)
        nm = p.name.lower()
        if "reward" in nm:
            reward_on.append(ons); reward_d.append(durs)
            for o in ons: onset_to_val[float(o)] = "Reward"
        elif "punish" in nm:
            punish_on.append(ons); punish_d.append(durs)
            for o in ons: onset_to_val[float(o)] = "Punish"
        else:
            neutral_on.append(ons); neutral_d.append(durs)
            for o in ons: onset_to_val[float(o)] = "Neutral"

    def cat(xs: List[np.ndarray]) -> np.ndarray:
        return np.sort(np.concatenate(xs)) if any(x.size for x in xs) else np.array([])

    def catd(xs: List[np.ndarray]) -> np.ndarray:
        return np.concatenate(xs) if any(x.size for x in xs) else np.array([])

    out["sr_outcome_primary"] = {
        "Reward": (cat(reward_on), catd(reward_d)),
        "Neutral": (cat(neutral_on), catd(neutral_d)),
        "Punish": (cat(punish_on), catd(punish_d)),
    }

    # decision split by subsequent outcome valence
    d_all = cat([d_comp[0], d_str[0]])
    all_outs = np.sort(np.array(list(onset_to_val.keys()))) if onset_to_val else np.array([])
    bins = {"Reward": [], "Neutral": [], "Punish": []}
    if d_all.size and all_outs.size:
        for d in d_all:
            j = np.searchsorted(all_outs, d, side="right")
            if j < all_outs.size:
                o = float(all_outs[j])
                if (o - d) < SR_MAX_DEC_TO_OUT_GAP_S:
                    bins[onset_to_val[o]].append(float(d))

    # durations for decision split: use decision durations (same EVs), not outcome
    d_durs = np.concatenate([d_comp[1], d_str[1]]) if (d_comp[1].size or d_str[1].size) else np.array([])
    out["sr_decision_by_future_outcome_extra"] = {
        "Reward": (np.array(sorted(bins["Reward"])), d_durs),
        "Neutral": (np.array(sorted(bins["Neutral"])), d_durs),
        "Punish": (np.array(sorted(bins["Punish"])), d_durs),
    }

    return out


# ----------------------------- marker helpers -----------------------------

def median_next_delta(a: np.ndarray, b: np.ndarray, max_delta: float = 30.0) -> float:
    """
    For each onset in a, find the first onset in b that occurs after it.
    Return median(b - a). NaN if none.
    """
    if a.size == 0 or b.size == 0:
        return float("nan")
    a = np.sort(a.astype(float))
    b = np.sort(b.astype(float))
    deltas = []
    for x in a:
        j = np.searchsorted(b, x, side="right")
        if j < b.size:
            d = float(b[j] - x)
            if d > 0 and d <= max_delta:
                deltas.append(d)
    return float(np.median(deltas)) if deltas else float("nan")


def median_duration(durs: List[np.ndarray]) -> float:
    all_d = np.concatenate([d for d in durs if d.size]) if any(d.size for d in durs) else np.array([])
    if all_d.size == 0:
        return float("nan")
    return float(np.median(all_d))


def compute_markers(task: str, plot_group: str,
                    ev_sets: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]) -> Tuple[float, float]:
    """
    Returns (event_duration_med, next_event_onset_med) in seconds from onset.
    next_event_onset_med is relative to the same onset used for the plot group.
    """
    if plot_group not in ev_sets:
        return float("nan"), float("nan")

    group = ev_sets[plot_group]
    group_onsets = [v[0] for v in group.values()]
    group_durs = [v[1] for v in group.values()]
    dur_med = median_duration(group_durs)

    if task == "mid":
        all_ant = np.sort(np.concatenate([
            ev_sets["mid_anticipation_primary"]["Reward"][0],
            ev_sets["mid_anticipation_primary"]["Neutral"][0],
        ]))
        fb_pos = ev_sets["mid_feedback_extra"]["Positive"][0]
        fb_neg = ev_sets["mid_feedback_extra"]["Negative"][0]
        all_fb = np.sort(np.concatenate([fb_pos, fb_neg])) if (fb_pos.size or fb_neg.size) else np.array([])

        a = np.sort(np.concatenate(group_onsets)) if any(x.size for x in group_onsets) else np.array([])
        if plot_group == "mid_anticipation_primary":
            next_med = median_next_delta(a, all_fb)
        else:
            next_med = median_next_delta(a, all_ant)
        return dur_med, next_med

    # sharedreward
    all_dec = np.sort(np.concatenate([
        ev_sets["sr_decision_partner_extra"]["Computer"][0],
        ev_sets["sr_decision_partner_extra"]["Stranger"][0],
    ])) if (ev_sets["sr_decision_partner_extra"]["Computer"][0].size or ev_sets["sr_decision_partner_extra"]["Stranger"][0].size) else np.array([])

    out_group = ev_sets["sr_outcome_primary"]
    all_out = np.sort(np.concatenate([
        out_group["Reward"][0], out_group["Neutral"][0], out_group["Punish"][0]
    ])) if any(out_group[k][0].size for k in ["Reward","Neutral","Punish"]) else np.array([])

    a = np.sort(np.concatenate(group_onsets)) if any(x.size for x in group_onsets) else np.array([])
    if plot_group == "sr_outcome_primary":
        next_med = median_next_delta(a, all_dec, max_delta=30.0)
    else:
        next_med = median_next_delta(a, all_out, max_delta=30.0)

    return dur_med, next_med


# ----------------------------- datamodel -----------------------------

@dataclass
class RunCurve:
    sub: str
    ses: str
    run: str
    task: str
    echo: str
    cnfds: str
    psc_method: str
    plot_group: str
    condition: str
    t: np.ndarray
    mean: np.ndarray
    sem: np.ndarray
    n_events: int
    event_dur_med_s: float
    next_onset_med_s: float
    feat_dir: str


# ----------------------------- core processing -----------------------------

def process_one_feat(feat: Path, mask_bool: np.ndarray, ev_base: Path) -> List[RunCurve]:
    sub, ses, run, task = parse_sub_ses_run_task(feat)
    echo, cnfds = infer_echo_and_cnfds(feat)

    func = feat / "filtered_func_data.nii.gz"
    if not func.exists():
        return []

    ts, tr = roi_mean_ts(func, mask_bool)
    ts = regress_confounds(ts, feat / "confoundevs.txt")
    ts = to_wu_psc(ts)

    ev_dir = ev_base / f"sub-{sub}" / f"ses-{ses}" / task / f"run-{run}"

    if task == "mid":
        ev_sets = mid_ev_sets(ev_dir)
    else:
        ev_sets = sr_ev_sets(ev_dir)

    # markers per plot_group
    markers: Dict[str, Tuple[float, float]] = {}
    for pg in ev_sets.keys():
        markers[pg] = compute_markers(task, pg, ev_sets)

    out: List[RunCurve] = []
    for plot_group, condmap in ev_sets.items():
        dur_med, next_med = markers.get(plot_group, (float("nan"), float("nan")))
        for cond, (onsets, _durs) in condmap.items():
            t_axis, wins = epoch_interpolate(ts, onsets, tr, TMIN, TMAX, DT)
            mu, se, n = mean_sem_nanaware(wins)
            out.append(RunCurve(
                sub=sub, ses=ses, run=run, task=task, echo=echo, cnfds=cnfds,
                psc_method=PSC_METHOD,
                plot_group=plot_group, condition=cond,
                t=t_axis, mean=mu, sem=se, n_events=n,
                event_dur_med_s=dur_med, next_onset_med_s=next_med,
                feat_dir=str(feat),
            ))
    return out


# ----------------------------- main -----------------------------

def main() -> None:
    pathlist = find_pathlist_file()
    feat_paths = read_feat_paths(pathlist)
    if not feat_paths:
        raise RuntimeError(f"Path list is empty: {pathlist}")

    root = infer_project_root_from_feat(feat_paths[0])
    mask_path = root / MASK_REL
    if not mask_path.exists():
        raise FileNotFoundError(f"ROI mask not found: {mask_path}")

    fsl_deriv = root / "derivatives" / "fsl"
    ev_base = fsl_deriv / "EVFiles"

    out_base = root / "derivatives" / "extractions" / "timecourses-unified-wu_noargs_v2"
    out_tables = out_base / "tables"
    out_runplots = out_base / "run_qc"
    out_session_qc = out_base / "session_qc"
    out_subject_qc = out_base / "subject_qc"

    out_tables.mkdir(parents=True, exist_ok=True)

    mask_bool = load_mask(mask_path)

    all_curves: List[RunCurve] = []
    n_ok = 0
    for i, feat in enumerate(feat_paths, start=1):
        if not feat.exists():
            print(f"[WARN] Missing FEAT dir, skipping: {feat}")
            continue
        try:
            curves = process_one_feat(feat, mask_bool, ev_base)
            if curves:
                all_curves.extend(curves)
                n_ok += 1
        except Exception as e:
            print(f"[WARN] Failed {feat}: {e}")

        if i % 25 == 0:
            print(f"[INFO] {i}/{len(feat_paths)} FEATs processed...")

    print(f"[INFO] Completed {n_ok}/{len(feat_paths)} FEATs with curves.")

    # ----------------- Write long table (gz) -----------------
    out_csv = out_tables / "run_curves_long_wu.csv.gz"
    with gzip.open(out_csv, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "sub","ses","run","task","echo","cnfds","psc_method",
            "plot_group","condition",
            "time_s","mean","sem","n_events",
            "event_dur_med_s","next_onset_med_s",
            "feat_dir"
        ])
        for c in all_curves:
            if c.mean.size == 0:
                continue
            for ti, mi, si in zip(c.t, c.mean, c.sem):
                w.writerow([
                    c.sub, c.ses, c.run, c.task, c.echo, c.cnfds, c.psc_method,
                    c.plot_group, c.condition,
                    float(ti), float(mi), float(si), int(c.n_events),
                    float(c.event_dur_med_s), float(c.next_onset_med_s),
                    c.feat_dir
                ])
    print(f"[INFO] Wrote: {out_csv}")

# ------------------ COUNT SUMMARY ------------------
try:
    from collections import Counter
    c_task = Counter((c.task for c in all_curves))
    c_pg = Counter(((c.task, c.plot_group) for c in all_curves))
    print("[INFO] COUNT SUMMARY: curves by task:", dict(c_task))
    print("[INFO] COUNT SUMMARY: curves by (task, plot_group):")
    for (t, pg), n in sorted(c_pg.items()):
        print(f"  - {t:12s} {pg:30s} {n}")
except Exception:
    pass


    # ----------------- Run-level QC plots -----------------
    # Group curves by run + plot_group
    run_bucket: Dict[Tuple[str,str,str,str,str,str,str], Dict[str, Tuple[np.ndarray,np.ndarray,int,float,float]]] = {}
    for c in all_curves:
        k = (c.sub, c.ses, c.run, c.task, c.echo, c.cnfds, c.plot_group)
        run_bucket.setdefault(k, {})[c.condition] = (c.mean, c.sem, c.n_events, c.event_dur_med_s, c.next_onset_med_s)

    for (sub, ses, run, task, echo, cnfds, plot_group), condmap in run_bucket.items():
        # use markers from any condition (they're plot_group-level)
        anyc = next(iter(condmap.values()))
        dur, nxt = anyc[3], anyc[4]
        t = all_curves[0].t  # same for all
        curves_dict = {cond: (vals[0], vals[1], vals[2]) for cond, vals in condmap.items()}
        fname = (
            f"runQC_task-{task}_sub-{sub}_ses-{ses}_run-{run}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_"
            f"psc-{PSC_METHOD}_{plot_group}.png"
        )
        out_png = out_runplots / fname
        title = f"RUN QC: {task} {plot_group} — sub-{sub} ses-{ses} run-{run} {echo} cnfds-{cnfds} ({PSC_METHOD})"
        plot_curves_with_sem(out_png, title, t, curves_dict, dur, nxt)

    print(f"[INFO] Wrote run QC under: {out_runplots}")

    # ----------------- Session-level QC (average across runs) -----------------
    sess_bucket: Dict[Tuple[str,str,str,str,str,str,str,str], List[np.ndarray]] = {}
    sess_markers: Dict[Tuple[str,str,str,str,str,str,str], List[Tuple[float,float]]] = {}

    # We'll aggregate within (sub, ses, task, echo, cnfds, plot_group, condition)
    for c in all_curves:
        k = (c.sub, c.ses, c.task, c.echo, c.cnfds, c.psc_method, c.plot_group, c.condition)
        if c.mean.size:
            sess_bucket.setdefault(k, []).append(c.mean)
        k2 = (c.sub, c.ses, c.task, c.echo, c.cnfds, c.psc_method, c.plot_group)
        sess_markers.setdefault(k2, []).append((c.event_dur_med_s, c.next_onset_med_s))

    # Build session plots dictionary: key=(sub,ses,task,echo,cnfds,psc,plot_group) -> cond->(mu,sem,n)
    sess_plots: Dict[Tuple[str,str,str,str,str,str,str], Dict[str, Tuple[np.ndarray,np.ndarray,int]]] = {}
    for (sub, ses, task, echo, cnfds, psc, plot_group, cond), arrs in sess_bucket.items():
        stack = [a for a in arrs if a.size]
        if not stack:
            continue
        mat = np.vstack(stack)
        mu = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0, ddof=1)
        sem = sd / np.sqrt(mat.shape[0])
        sess_plots.setdefault((sub, ses, task, echo, cnfds, psc, plot_group), {})[cond] = (mu, sem, mat.shape[0])

    sess_plot_markers: Dict[Tuple[str,str,str,str,str,str,str], Tuple[float,float]] = {}
    for k2, vals in sess_markers.items():
        durs = np.array([v[0] for v in vals], float)
        nxts = np.array([v[1] for v in vals], float)
        dur = float(np.nanmedian(durs)) if np.any(np.isfinite(durs)) else float("nan")
        nxt = float(np.nanmedian(nxts)) if np.any(np.isfinite(nxts)) else float("nan")
        sess_plot_markers[k2] = (dur, nxt)

    t = all_curves[0].t if all_curves else np.arange(TMIN, TMAX + 1e-9, DT)
    for (sub, ses, task, echo, cnfds, psc, plot_group), condmap in sess_plots.items():
        dur, nxt = sess_plot_markers.get((sub, ses, task, echo, cnfds, psc, plot_group), (float("nan"), float("nan")))
        fname = (
            f"sessionQC_task-{task}_sub-{sub}_ses-{ses}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_"
            f"psc-{psc}_{plot_group}.png"
        )
        out_png = out_session_qc / fname
        title = f"SESSION QC (avg across runs): {task} {plot_group} — sub-{sub} ses-{ses} {echo} cnfds-{cnfds} ({psc})"
        plot_curves_with_sem(out_png, title, t, condmap, dur, nxt)

    print(f"[INFO] Wrote session QC under: {out_session_qc}")

    # ----------------- Subject-level QC (average across sessions) -----------------
    subj_bucket: Dict[Tuple[str,str,str,str,str,str,str], List[np.ndarray]] = {}
    subj_markers: Dict[Tuple[str,str,str,str,str,str], List[Tuple[float,float]]] = {}

    for (sub, ses, task, echo, cnfds, psc, plot_group), condmap in sess_plots.items():
        for cond, (mu, _sem, _n) in condmap.items():
            k = (sub, task, echo, cnfds, psc, plot_group, cond)
            subj_bucket.setdefault(k, []).append(mu)
        k2 = (sub, task, echo, cnfds, psc, plot_group)
        subj_markers.setdefault(k2, []).append(sess_plot_markers.get((sub, ses, task, echo, cnfds, psc, plot_group), (float("nan"), float("nan"))))

    subj_plots: Dict[Tuple[str,str,str,str,str,str], Dict[str, Tuple[np.ndarray,np.ndarray,int]]] = {}
    for (sub, task, echo, cnfds, psc, plot_group, cond), arrs in subj_bucket.items():
        stack = [a for a in arrs if a.size]
        if not stack:
            continue
        mat = np.vstack(stack)
        mu = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0, ddof=1)
        sem = sd / np.sqrt(mat.shape[0])
        subj_plots.setdefault((sub, task, echo, cnfds, psc, plot_group), {})[cond] = (mu, sem, mat.shape[0])

    subj_plot_markers: Dict[Tuple[str,str,str,str,str,str], Tuple[float,float]] = {}
    for k2, vals in subj_markers.items():
        durs = np.array([v[0] for v in vals], float)
        nxts = np.array([v[1] for v in vals], float)
        dur = float(np.nanmedian(durs)) if np.any(np.isfinite(durs)) else float("nan")
        nxt = float(np.nanmedian(nxts)) if np.any(np.isfinite(nxts)) else float("nan")
        subj_plot_markers[k2] = (dur, nxt)

    for (sub, task, echo, cnfds, psc, plot_group), condmap in subj_plots.items():
        dur, nxt = subj_plot_markers.get((sub, task, echo, cnfds, psc, plot_group), (float("nan"), float("nan")))
        fname = (
            f"subjectQC_task-{task}_sub-{sub}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_"
            f"psc-{psc}_{plot_group}.png"
        )
        out_png = out_subject_qc / fname
        title = f"SUBJECT QC (avg across sessions): {task} {plot_group} — sub-{sub} {echo} cnfds-{cnfds} ({psc})"
        plot_curves_with_sem(out_png, title, t, condmap, dur, nxt)

    print(f"[INFO] Wrote subject QC under: {out_subject_qc}")
    print("[INFO] Done.")


if __name__ == "__main__":
    main()
