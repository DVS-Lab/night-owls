#!/usr/bin/env python3
from __future__ import annotations

import re
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# ----------------------------- CONFIG (baked-in) -----------------------------

SCRIPT_DIR = Path(__file__).resolve().parent

PATHLIST_CANDIDATES = [
    SCRIPT_DIR / "feat_paths_704.txt",
    SCRIPT_DIR / "feat_paths_all.txt",
]

MASK_REL = Path("masks/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz")

DT = 0.1
TMIN = -6.0
TMAX = 16.0

TR_EXPECTED = 1.615
TR_TOL = 0.01

BASELINE_NTRS = [2, 3]
SR_MAX_DEC_TO_OUT_GAP_S = 8.0

DPI = 150

# Robust y-limits so a single outlier run doesn't flatten everything.
# (These are quantiles of the plotted lines; adjust if you want more/less zoom.)
YLIM_QLO = 0.01
YLIM_QHI = 0.99
YLIM_PAD_FRAC = 0.10

# ----------------------------- helpers -----------------------------

FEAT_NAME_RE = re.compile(r".*task-(mid|sharedreward).*run-(\d+).*", re.IGNORECASE)


def find_pathlist_file() -> Path:
    for p in PATHLIST_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError(
        "No path-list file found next to the script. Expected one of:\n"
        + "\n".join([f"  - {x.name}" for x in PATHLIST_CANDIDATES])
    )


def read_feat_paths(pathlist: Path) -> List[Path]:
    out: List[Path] = []
    for ln in pathlist.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        out.append(Path(s))
    return out


def infer_project_root_from_feat(feat: Path) -> Path:
    ps = feat.resolve().parents
    # expected: <root>/derivatives/fsl/sub-*/ses-*/<feat>.feat
    if len(ps) >= 6 and ps[2].name.startswith("sub-") and ps[1].name.startswith("ses-"):
        if ps[3].name == "fsl" and ps[4].name == "derivatives":
            return ps[5]
    # fallback: walk upward until we find "derivatives"
    cur = feat.resolve()
    for parent in cur.parents:
        if (parent / "derivatives").exists():
            return parent
    raise RuntimeError(f"Could not infer project root from feat path: {feat}")


def parse_sub_ses_run_task(feat: Path) -> Tuple[str, str, str, str]:
    s = feat.as_posix()
    m = re.search(r"sub-(\d+).*ses-(\d+).*run-(\d+)", s)
    if not m:
        raise ValueError(f"Could not parse sub/ses/run from: {feat}")
    sub, ses, run = m.group(1), m.group(2), m.group(3)

    mn = FEAT_NAME_RE.match(feat.name)
    if mn:
        task = mn.group(1).lower()
    else:
        if "task-mid" in s:
            task = "mid"
        elif "task-sharedreward" in s:
            task = "sharedreward"
        else:
            raise ValueError(f"Could not infer task (mid/sharedreward) from: {feat}")

    return sub, ses, run, task


def infer_echo_and_cnfds(feat: Path) -> Tuple[str, str]:
    name = feat.name
    echo = "multi-echo" if "multi-echo" in name else ("single-echo" if "single-echo" in name else "unknown-echo")
    m = re.search(r"cnfds-([^_]+)", name)
    cnfds = m.group(1) if m else "unknown-cnfds"
    return echo, cnfds


def load_mask(mask_path: Path) -> np.ndarray:
    m = nib.load(str(mask_path)).get_fdata()
    return (m > 0)


def regress_confounds(ts: np.ndarray, confound_file: Path) -> np.ndarray:
    """Regress confound columns (+ intercept) from ts, if present and dimension-matched."""
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
    if abs(tr - TR_EXPECTED) > TR_TOL:
        tr = TR_EXPECTED
    return ts, tr


def load_ev_3col(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads an FSL 3-column EV file: onset, duration, amplitude.
    Returns (onsets, durations). Missing/empty => ([], []).
    If file has 1 column: duration=0. If 2 columns: amplitude ignored.
    """
    if not path.exists():
        return np.array([]), np.array([])
    try:
        arr = np.loadtxt(path, ndmin=2)
    except Exception:
        return np.array([]), np.array([])
    if arr.size == 0:
        return np.array([]), np.array([])
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    on = arr[:, 0].astype(float)
    if arr.shape[1] >= 2:
        dur = arr[:, 1].astype(float)
    else:
        dur = np.zeros_like(on)
    return on, dur


def epoch_interpolate(ts: np.ndarray, onsets: np.ndarray, tr: float,
                      tmin: float, tmax: float, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    n_pts = int(np.round((tmax - tmin) / dt)) + 1
    t_axis = np.linspace(tmin, tmax, n_pts)
    if onsets.size == 0:
        return t_axis, np.empty((0, n_pts))

    t_raw = np.arange(ts.size) * tr
    wins = np.empty((onsets.size, n_pts), dtype=float)
    wins[:] = np.nan
    for i, onset in enumerate(onsets):
        t_sample = onset + t_axis
        wins[i, :] = np.interp(t_sample, t_raw, ts, left=np.nan, right=np.nan)
    return t_axis, wins


def epoch_psc_from_baseline(wins: np.ndarray, t_axis: np.ndarray,
                            base_window: Tuple[float, float]) -> np.ndarray:
    """
    PSC relative to per-epoch baseline: ((x/b) - 1) * 100
    Then subtract PSC baseline mean so baseline is exactly 0.
    """
    if wins.size == 0:
        return wins
    lo, hi = base_window
    idx = (t_axis >= lo) & (t_axis < hi)
    if not np.any(idx):
        return wins * np.nan

    base = np.nanmean(wins[:, idx], axis=1)
    good = np.isfinite(base) & (np.abs(base) > 1e-12)
    out = np.full_like(wins, np.nan)
    out[good, :] = ((wins[good, :] / base[good, None]) - 1.0) * 100.0

    # force exactly zero in baseline window
    bpsc = np.nanmean(out[:, idx], axis=1)
    out = out - bpsc[:, None]
    return out


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
        y = np.asarray(line.get_ydata())
        if y.size:
            y = y[np.isfinite(y)]
            if y.size:
                ys.append(y)
    if not ys:
        return
    yall = np.concatenate(ys, axis=0)
    if yall.size == 0:
        return
    lo = float(np.quantile(yall, YLIM_QLO))
    hi = float(np.quantile(yall, YLIM_QHI))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return
    if lo == hi:
        pad = 0.25 if lo == 0 else abs(lo) * 0.1
        ax.set_ylim(lo - pad, hi + pad)
        return
    pad = (hi - lo) * YLIM_PAD_FRAC
    ax.set_ylim(lo - pad, hi + pad)


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


# ----------------------------- EV logic + markers -----------------------------

def mid_ev_sets(ev_dir: Path) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """
    Returns mapping plot_group -> condition -> (onsets, durations)
    """
    ant_r = load_ev_3col(ev_dir / "_anticipation_reward.txt")
    ant_n = load_ev_3col(ev_dir / "_anticipation_neutral.txt")

    fb_pos_r = load_ev_3col(ev_dir / "_feedback_positive_reward.txt")
    fb_neg_r = load_ev_3col(ev_dir / "_feedback_negative_reward.txt")
    fb_pos_n = load_ev_3col(ev_dir / "_feedback_positive_neutral.txt")
    fb_neg_n = load_ev_3col(ev_dir / "_feedback_negative_neutral.txt")

    out: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    out["mid_anticipation_primary"] = {"Reward": ant_r, "Neutral": ant_n}

    # feedback collapsed across cue types, split by valence
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


def compute_markers(task: str, plot_group: str,
                    ev_sets: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]) -> Tuple[float, float]:
    """
    Returns (event_duration_med, next_event_onset_med) in seconds from onset.
    next_event_onset_med is relative to the *same onset* used for the plot group.
    """
    if plot_group not in ev_sets:
        return float("nan"), float("nan")

    group = ev_sets[plot_group]
    group_onsets = [v[0] for v in group.values()]
    group_durs = [v[1] for v in group.values()]
    dur_med = median_duration(group_durs)

    if task == "mid":
        all_ant = np.sort(np.concatenate([ev_sets["mid_anticipation_primary"]["Reward"][0],
                                          ev_sets["mid_anticipation_primary"]["Neutral"][0]]))
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
    all_dec = np.sort(np.concatenate([ev_sets["sr_decision_partner_extra"]["Computer"][0],
                                      ev_sets["sr_decision_partner_extra"]["Stranger"][0]]))
    out_group = ev_sets["sr_outcome_primary"]
    all_out = np.sort(np.concatenate([out_group["Reward"][0], out_group["Neutral"][0], out_group["Punish"][0]]))

    a = np.sort(np.concatenate(group_onsets)) if any(x.size for x in group_onsets) else np.array([])
    if plot_group == "sr_outcome_primary":
        next_med = median_next_delta(a, all_dec)
    else:
        next_med = median_next_delta(a, all_out)

    return dur_med, next_med


# ----------------------------- data structures -----------------------------

@dataclass
class RunCurve:
    sub: str
    ses: str
    run: str
    task: str
    echo: str
    cnfds: str
    baseline_ntr: int
    plot_group: str
    condition: str
    t: np.ndarray
    mean: np.ndarray
    sem: np.ndarray
    n_events: int
    feat_dir: str
    event_dur_med: float
    next_onset_med: float


# ----------------------------- plotting -----------------------------

def add_event_markers(ax, event_dur: float, next_onset: float) -> None:
    # onset always at 0
    ax.axvline(0.0, color="k", linestyle=":", linewidth=1)

    # duration line
    if np.isfinite(event_dur) and event_dur > 0.05:
        ax.axvline(event_dur, color="0.4", linestyle="--", linewidth=1)
        ax.text(event_dur, ax.get_ylim()[1], " end", va="top", ha="left", fontsize=8, color="0.3")

    # next event onset
    if np.isfinite(next_onset) and next_onset > 0.05:
        ax.axvline(next_onset, color="0.4", linestyle="-.", linewidth=1)
        ax.text(next_onset, ax.get_ylim()[1], " next", va="top", ha="left", fontsize=8, color="0.3")


def plot_curves_with_sem(out_png: Path, title: str, t: np.ndarray,
                         curves: Dict[str, Tuple[np.ndarray, np.ndarray, int]],
                         event_dur: float, next_onset: float) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    for name in sorted(curves.keys()):
        mu, se, n = curves[name]
        if mu.size == 0:
            continue
        ax.plot(t, mu, label=f"{name} (n={n})")
        ax.fill_between(t, mu - se, mu + se, alpha=0.15)

    ax.set_title(title)
    ax.set_xlabel("Time from event onset (s)")
    ax.set_ylabel("PSC (epoch-baseline referenced)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, ncol=2)

    # y-lims first (so marker labels land properly)
    robust_ylim_from_lines(ax)
    add_event_markers(ax, event_dur, next_onset)

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)


def plot_compare_baselines(out_png: Path, title: str, t: np.ndarray,
                           curves_base2: Dict[str, Tuple[np.ndarray, np.ndarray, int]],
                           curves_base3: Dict[str, Tuple[np.ndarray, np.ndarray, int]],
                           markers_base2: Tuple[float, float],
                           markers_base3: Tuple[float, float]) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    for ax, (lbl, curves, (dur, nxt)) in zip(
        axes,
        [
            ("Baseline = 2 TR", curves_base2, markers_base2),
            ("Baseline = 3 TR", curves_base3, markers_base3),
        ],
    ):
        for name in sorted(curves.keys()):
            mu, se, n = curves[name]
            if mu.size == 0:
                continue
            ax.plot(t, mu, label=f"{name} (n={n})")
            ax.fill_between(t, mu - se, mu + se, alpha=0.15)
        ax.set_title(lbl)
        ax.set_xlabel("Time from event onset (s)")
        ax.grid(True, alpha=0.25)

        robust_ylim_from_lines(ax)
        add_event_markers(ax, dur, nxt)

    axes[0].set_ylabel("PSC (epoch-baseline referenced)")
    axes[1].legend(fontsize=8, ncol=1)

    fig.suptitle(title, y=0.98)
    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)


# ----------------------------- core processing -----------------------------

def process_one_feat(feat: Path, mask_bool: np.ndarray, ev_base: Path) -> List[RunCurve]:
    sub, ses, run, task = parse_sub_ses_run_task(feat)
    echo, cnfds = infer_echo_and_cnfds(feat)

    func = feat / "filtered_func_data.nii.gz"
    if not func.exists():
        return []

    ts, tr = roi_mean_ts(func, mask_bool)
    ts = regress_confounds(ts, feat / "confoundevs.txt")

    ev_dir = ev_base / f"sub-{sub}" / f"ses-{ses}" / task / f"run-{run}"

    if task == "mid":
        ev_sets = mid_ev_sets(ev_dir)
    else:
        ev_sets = sr_ev_sets(ev_dir)

    # markers per plot_group (independent of baseline and condition)
    markers: Dict[str, Tuple[float, float]] = {}
    for pg in ev_sets.keys():
        markers[pg] = compute_markers(task, pg, ev_sets)

    out: List[RunCurve] = []
    for baseline_ntr in BASELINE_NTRS:
        base_window = (-baseline_ntr * tr, 0.0)

        for plot_group, condmap in ev_sets.items():
            dur_med, next_med = markers.get(plot_group, (float("nan"), float("nan")))
            for cond, (onsets, _durs) in condmap.items():
                t_axis, wins = epoch_interpolate(ts, onsets, tr, TMIN, TMAX, DT)
                wins_psc = epoch_psc_from_baseline(wins, t_axis, base_window)
                mu, se, n = mean_sem_nanaware(wins_psc)

                out.append(
                    RunCurve(
                        sub=sub, ses=ses, run=run, task=task,
                        echo=echo, cnfds=cnfds,
                        baseline_ntr=baseline_ntr, plot_group=plot_group,
                        condition=cond, t=t_axis, mean=mu, sem=se, n_events=n,
                        feat_dir=str(feat),
                        event_dur_med=dur_med,
                        next_onset_med=next_med,
                    )
                )
    return out


def write_run_curves_long_csv(out_csv: Path, curves: List[RunCurve]) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "sub","ses","run","task","echo","cnfds",
            "baseline_ntr","plot_group","condition",
            "time_s","mean","sem","n_events",
            "event_dur_med_s","next_onset_med_s",
            "feat_dir"
        ])
        for c in curves:
            if c.mean.size == 0:
                continue
            for ti, mi, si in zip(c.t, c.mean, c.sem):
                w.writerow([
                    c.sub, c.ses, c.run, c.task, c.echo, c.cnfds,
                    c.baseline_ntr, c.plot_group, c.condition,
                    f"{ti:.4f}", f"{mi:.8f}", f"{si:.8f}", c.n_events,
                    f"{c.event_dur_med:.4f}" if np.isfinite(c.event_dur_med) else "",
                    f"{c.next_onset_med:.4f}" if np.isfinite(c.next_onset_med) else "",
                    c.feat_dir
                ])


def safe_slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9\-]+", "-", s).strip("-")


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

    # new output root to avoid mixing with your previous run
    out_base = root / "derivatives" / "extractions" / "timecourses-unified-noargs_v2"
    out_tables = out_base / "tables"
    out_runplots = out_base / "run_plots"
    out_session_qc = out_base / "session_qc"
    out_subject_qc = out_base / "subject_qc"
    out_compare_session = out_base / "session_qc_compare_baselines"
    out_compare_subject = out_base / "subject_qc_compare_baselines"

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
            print(f"[WARN] Failed on {feat}: {e}")

        if i % 25 == 0:
            print(f"[INFO] Processed {i}/{len(feat_paths)} paths (ok={n_ok})")

    print(f"[INFO] Finished run processing. ok={n_ok} / total={len(feat_paths)}")

    out_csv = out_tables / "run_curves_long.csv"
    write_run_curves_long_csv(out_csv, all_curves)
    print(f"[INFO] Wrote: {out_csv}")

    # ------------------ RUN PLOTS (flat filenames) ------------------
    run_groups: Dict[Tuple[str,str,str,str,str,str,int,str], Dict[str, RunCurve]] = {}
    for c in all_curves:
        k = (c.sub, c.ses, c.run, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group)
        run_groups.setdefault(k, {})[c.condition] = c

    for (sub, ses, run, task, echo, cnfds, bntr, plot_group), condmap in run_groups.items():
        anyc = next(iter(condmap.values()))
        t = anyc.t
        curves_dict = {cond: (obj.mean, obj.sem, obj.n_events) for cond, obj in condmap.items()}
        dur, nxt = anyc.event_dur_med, anyc.next_onset_med

        fname = (
            f"run_task-{task}_sub-{sub}_ses-{ses}_run-{run}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_"
            f"baseline-{bntr}TR_{plot_group}.png"
        )
        out_png = out_runplots / fname
        title = f"{task} {plot_group} — sub-{sub} ses-{ses} run-{run} {echo} cnfds-{cnfds} (baseline={bntr}TR)"
        plot_curves_with_sem(out_png, title, t, curves_dict, dur, nxt)

    print(f"[INFO] Wrote run plots under: {out_runplots}")

    # ------------------ SESSION QC (summary across runs; includes ALL plot_groups) ------------------
    sess_bucket: Dict[Tuple[str,str,str,str,str,int,str,str], List[np.ndarray]] = {}
    sess_markers: Dict[Tuple[str,str,str,str,str,int,str], List[Tuple[float,float]]] = {}
    sess_taxis: Dict[Tuple[str,str,str,str,str,int,str], np.ndarray] = {}

    for c in all_curves:
        k = (c.sub, c.ses, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group, c.condition)
        sess_bucket.setdefault(k, []).append(c.mean)

        k2 = (c.sub, c.ses, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group)
        sess_taxis[k2] = c.t
        sess_markers.setdefault(k2, []).append((c.event_dur_med, c.next_onset_med))

    sess_plots: Dict[Tuple[str,str,str,str,str,int,str], Dict[str, Tuple[np.ndarray,np.ndarray,int]]] = {}
    for (sub, ses, task, echo, cnfds, bntr, plot_group, cond), arrs in sess_bucket.items():
        stack = [a for a in arrs if a.size]
        if not stack:
            continue
        mat = np.vstack(stack)
        mu = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0, ddof=1)
        sem = sd / np.sqrt(mat.shape[0])
        sess_plots.setdefault((sub, ses, task, echo, cnfds, bntr, plot_group), {})[cond] = (mu, sem, mat.shape[0])

    sess_plot_markers: Dict[Tuple[str,str,str,str,str,int,str], Tuple[float,float]] = {}
    for k2, vals in sess_markers.items():
        durs = np.array([v[0] for v in vals], float)
        nxts = np.array([v[1] for v in vals], float)
        dur_med = float(np.nanmedian(durs)) if np.any(np.isfinite(durs)) else float("nan")
        nxt_med = float(np.nanmedian(nxts)) if np.any(np.isfinite(nxts)) else float("nan")
        sess_plot_markers[k2] = (dur_med, nxt_med)

    for (sub, ses, task, echo, cnfds, bntr, plot_group), condmap in sess_plots.items():
        t = sess_taxis[(sub, ses, task, echo, cnfds, bntr, plot_group)]
        dur, nxt = sess_plot_markers.get((sub, ses, task, echo, cnfds, bntr, plot_group), (float("nan"), float("nan")))

        fname = (
            f"sessionQC_task-{task}_sub-{sub}_ses-{ses}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_"
            f"baseline-{bntr}TR_{plot_group}_summary.png"
        )
        out_png = out_session_qc / fname
        title = f"SESSION QC: {task} {plot_group} — sub-{sub} ses-{ses} {echo} cnfds-{cnfds} (baseline={bntr}TR)"
        plot_curves_with_sem(out_png, title, t, condmap, dur, nxt)

    print(f"[INFO] Wrote session QC under: {out_session_qc}")

    # ------------------ SESSION baseline compare (ALL plot_groups) ------------------
    sess_compare: Dict[Tuple[str,str,str,str,str,str], Dict[int, Dict[str, Tuple[np.ndarray,np.ndarray,int]]]] = {}
    sess_compare_markers: Dict[Tuple[str,str,str,str,str,str], Dict[int, Tuple[float,float]]] = {}
    sess_compare_t: Dict[Tuple[str,str,str,str,str,str], np.ndarray] = {}

    for (sub, ses, task, echo, cnfds, bntr, plot_group), condmap in sess_plots.items():
        k = (sub, ses, task, echo, cnfds, plot_group)
        sess_compare.setdefault(k, {})[bntr] = condmap
        sess_compare_t[k] = sess_taxis[(sub, ses, task, echo, cnfds, bntr, plot_group)]
        sess_compare_markers.setdefault(k, {})[bntr] = sess_plot_markers.get((sub, ses, task, echo, cnfds, bntr, plot_group), (float("nan"), float("nan")))

    for k, by_base in sess_compare.items():
        if 2 not in by_base or 3 not in by_base:
            continue
        sub, ses, task, echo, cnfds, plot_group = k
        t = sess_compare_t[k]
        m2 = sess_compare_markers.get(k, {}).get(2, (float("nan"), float("nan")))
        m3 = sess_compare_markers.get(k, {}).get(3, (float("nan"), float("nan")))

        fname = (
            f"sessionQC_compare_baseline2vs3_task-{task}_sub-{sub}_ses-{ses}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_{plot_group}.png"
        )
        out_png = out_compare_session / fname
        title = f"SESSION QC (baseline compare): {task} {plot_group} — sub-{sub} ses-{ses} {echo} cnfds-{cnfds}"
        plot_compare_baselines(out_png, title, t, by_base[2], by_base[3], m2, m3)

    print(f"[INFO] Wrote session baseline-compare QC under: {out_compare_session}")

    # ------------------ SUBJECT QC (summary across all runs/sessions; ALL plot_groups) ------------------
    subj_bucket: Dict[Tuple[str,str,str,str,int,str,str], List[np.ndarray]] = {}
    subj_markers: Dict[Tuple[str,str,str,str,int,str], List[Tuple[float,float]]] = {}
    subj_taxis: Dict[Tuple[str,str,str,str,int,str], np.ndarray] = {}

    for c in all_curves:
        k = (c.sub, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group, c.condition)
        subj_bucket.setdefault(k, []).append(c.mean)
        k2 = (c.sub, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group)
        subj_taxis[k2] = c.t
        subj_markers.setdefault(k2, []).append((c.event_dur_med, c.next_onset_med))

    subj_plots: Dict[Tuple[str,str,str,str,int,str], Dict[str, Tuple[np.ndarray,np.ndarray,int]]] = {}
    for (sub, task, echo, cnfds, bntr, plot_group, cond), arrs in subj_bucket.items():
        stack = [a for a in arrs if a.size]
        if not stack:
            continue
        mat = np.vstack(stack)
        mu = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0, ddof=1)
        sem = sd / np.sqrt(mat.shape[0])
        subj_plots.setdefault((sub, task, echo, cnfds, bntr, plot_group), {})[cond] = (mu, sem, mat.shape[0])

    subj_plot_markers: Dict[Tuple[str,str,str,str,int,str], Tuple[float,float]] = {}
    for k2, vals in subj_markers.items():
        durs = np.array([v[0] for v in vals], float)
        nxts = np.array([v[1] for v in vals], float)
        dur_med = float(np.nanmedian(durs)) if np.any(np.isfinite(durs)) else float("nan")
        nxt_med = float(np.nanmedian(nxts)) if np.any(np.isfinite(nxts)) else float("nan")
        subj_plot_markers[k2] = (dur_med, nxt_med)

    for (sub, task, echo, cnfds, bntr, plot_group), condmap in subj_plots.items():
        t = subj_taxis[(sub, task, echo, cnfds, bntr, plot_group)]
        dur, nxt = subj_plot_markers.get((sub, task, echo, cnfds, bntr, plot_group), (float("nan"), float("nan")))

        fname = (
            f"subjectQC_task-{task}_sub-{sub}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_"
            f"baseline-{bntr}TR_{plot_group}_summary.png"
        )
        out_png = out_subject_qc / fname
        title = f"SUBJECT QC: {task} {plot_group} — sub-{sub} {echo} cnfds-{cnfds} (baseline={bntr}TR)"
        plot_curves_with_sem(out_png, title, t, condmap, dur, nxt)

    print(f"[INFO] Wrote subject QC under: {out_subject_qc}")

    # ------------------ SUBJECT baseline compare (ALL plot_groups) ------------------
    subj_compare: Dict[Tuple[str,str,str,str,str], Dict[int, Dict[str, Tuple[np.ndarray,np.ndarray,int]]]] = {}
    subj_compare_markers: Dict[Tuple[str,str,str,str,str], Dict[int, Tuple[float,float]]] = {}
    subj_compare_t: Dict[Tuple[str,str,str,str,str], np.ndarray] = {}

    for (sub, task, echo, cnfds, bntr, plot_group), condmap in subj_plots.items():
        k = (sub, task, echo, cnfds, plot_group)
        subj_compare.setdefault(k, {})[bntr] = condmap
        subj_compare_t[k] = subj_taxis[(sub, task, echo, cnfds, bntr, plot_group)]
        subj_compare_markers.setdefault(k, {})[bntr] = subj_plot_markers.get((sub, task, echo, cnfds, bntr, plot_group), (float("nan"), float("nan")))

    for k, by_base in subj_compare.items():
        if 2 not in by_base or 3 not in by_base:
            continue
        sub, task, echo, cnfds, plot_group = k
        t = subj_compare_t[k]
        m2 = subj_compare_markers.get(k, {}).get(2, (float("nan"), float("nan")))
        m3 = subj_compare_markers.get(k, {}).get(3, (float("nan"), float("nan")))

        fname = (
            f"subjectQC_compare_baseline2vs3_task-{task}_sub-{sub}_"
            f"echo-{safe_slug(echo)}_cnfds-{safe_slug(cnfds)}_{plot_group}.png"
        )
        out_png = out_compare_subject / fname
        title = f"SUBJECT QC (baseline compare): {task} {plot_group} — sub-{sub} {echo} cnfds-{cnfds}"
        plot_compare_baselines(out_png, title, t, by_base[2], by_base[3], m2, m3)

    print(f"[INFO] Wrote subject baseline-compare QC under: {out_compare_subject}")
    print("[INFO] Done.")


if __name__ == "__main__":
    main()
