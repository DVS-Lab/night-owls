#!/usr/bin/env python3
"""
Unified peristimulus ROI timecourse extractor + QC plots (MID + SharedReward)
NO COMMAND-LINE ARGS. Everything is baked in.

Inputs
- A single text file of FEAT directories (~704 paths), located in the same folder as this script.
  Accepted filenames (first found wins):
    feat_paths_704.txt
    feat_paths_all.txt
    feat_paths.txt

Assumptions
- Each FEAT dir contains: filtered_func_data.nii.gz
- Optional confounds: confoundevs.txt (regressed from ROI mean time series if present)
- EVs live under: <project_root>/derivatives/fsl/EVFiles/sub-*/ses-*/<task>/run-*/
- Mask lives under: <project_root>/masks/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz
- Task is inferred from FEAT dir name containing "task-mid" or "task-sharedreward".

Baseline correction (two versions, both produced)
- Epoch-based PSC relative to per-epoch baseline:
    PSC(t) = ((x(t)/baseline) - 1) * 100
  where baseline = mean(x(t) for t in [-NTR*TR, 0)).
- Two baselines:
    NTR = 2 and NTR = 3
- Baseline anchoring:
    After PSC, subtract the PSC mean within baseline window so baseline is exactly 0.

QC summaries you asked for (single line + SEM)
- MID primary QC: anticipation Reward vs Neutral (2 conditions)
- SharedReward primary QC: outcome Reward vs Neutral vs Punish (3 conditions)
SEM for subject/session summaries is computed across RUNS (run-level mean curves),
not across trials, so it reflects within-person variability across runs/sessions.

Outputs (created automatically under project root)
<root>/derivatives/extractions/timecourses-unified-noargs/
  tables/run_curves_long.csv
  run_plots/<baseline>/...
  session_qc/<baseline>/...
  subject_qc/<baseline>/...
  session_qc/compare_baselines/...
  subject_qc/compare_baselines/...
"""

from __future__ import annotations

import re
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt


# ----------------------------- CONFIG (baked-in) -----------------------------

SCRIPT_DIR = Path(__file__).resolve().parent

PATHLIST_CANDIDATES = [
    SCRIPT_DIR / "feat_paths_704.txt",
    SCRIPT_DIR / "feat_paths-all.txt"
]

MASK_REL = Path("masks/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz")

# Interpolation grid for peristimulus sampling
DT = 0.1

# Must include the 3-TR baseline: if TR=1.615, 3TR=4.845 sec
TMIN = -6.0
TMAX = 16.0

# Expected TR (the project uses 1.615); if header differs slightly, clamp to this
TR_EXPECTED = 1.615
TR_TOL = 0.01

# Baseline windows (in TRs) to produce automatically
BASELINE_NTRS = [2, 3]

# SharedReward: optional decision->future outcome splitting sanity constraint
SR_MAX_DEC_TO_OUT_GAP_S = 8.0

# Plot DPI
DPI = 150


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
    lines = pathlist.read_text().splitlines()
    out = []
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        out.append(Path(s))
    return out


def infer_project_root_from_feat(feat: Path) -> Path:
    # expects .../<root>/derivatives/fsl/sub-*/ses-*/<feat>.feat
    # parents: feat, ses-*, sub-*, fsl, derivatives, root
    # so root is parents[5] if structure matches
    ps = feat.resolve().parents
    if len(ps) >= 6 and ps[2].name.startswith("sub-") and ps[1].name.startswith("ses-"):
        # parents[3] == fsl, parents[4] == derivatives
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

    # task from name
    mn = FEAT_NAME_RE.match(feat.name)
    if mn:
        task = mn.group(1).lower()
    else:
        # fallback: search in full path
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
    """
    Regress confound columns (+ intercept) from ts.
    Returns residuals + mean(ts) to preserve scale.
    """
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
        X = np.column_stack([X, np.ones(X.shape[0])])  # intercept
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


def load_ev_onsets(path: Path) -> np.ndarray:
    if not path.exists():
        return np.array([])
    try:
        arr = np.loadtxt(path, ndmin=2)
    except Exception:
        return np.array([])
    if arr.size == 0:
        return np.array([])
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0].astype(float)


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

    # force zero baseline
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


def set_data_driven_ylim(ax) -> None:
    # Use current line + fill artists to compute finite bounds
    ys = []
    for line in ax.get_lines():
        y = line.get_ydata()
        if y is not None:
            ys.append(np.asarray(y))
    if not ys:
        return
    yall = np.concatenate([y[np.isfinite(y)] for y in ys if np.any(np.isfinite(y))], axis=0)
    if yall.size == 0:
        return
    lo, hi = float(np.min(yall)), float(np.max(yall))
    if lo == hi:
        pad = 0.25 if lo == 0 else abs(lo) * 0.1
    else:
        pad = (hi - lo) * 0.1
    ax.set_ylim(lo - pad, hi + pad)


# ----------------------------- task EV logic -----------------------------

def evs_mid_primary(ev_dir: Path) -> Dict[str, np.ndarray]:
    # Primary MID QC: anticipation only (2 conditions)
    return {
        "Reward":  np.sort(load_ev_onsets(ev_dir / "_anticipation_reward.txt")),
        "Neutral": np.sort(load_ev_onsets(ev_dir / "_anticipation_neutral.txt")),
    }


def evs_mid_feedback(ev_dir: Path) -> Dict[str, np.ndarray]:
    # Optional extra: feedback collapsed across cue types (2 conditions)
    pos = np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_positive_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_positive_neutral.txt"),
    ])
    neg = np.concatenate([
        load_ev_onsets(ev_dir / "_feedback_negative_reward.txt"),
        load_ev_onsets(ev_dir / "_feedback_negative_neutral.txt"),
    ])
    return {
        "Positive": np.sort(pos) if pos.size else np.array([]),
        "Negative": np.sort(neg) if neg.size else np.array([]),
    }


def evs_sr_outcome_primary(ev_dir: Path) -> Dict[str, np.ndarray]:
    # Primary SR QC: outcome locked, pooled across partner by filename (3 conditions)
    reward, neutral, punish = [], [], []
    for p in ev_dir.glob("_outcome_*.txt"):
        ons = load_ev_onsets(p)
        nm = p.name.lower()
        if "reward" in nm:
            reward.append(ons)
        elif "punish" in nm:
            punish.append(ons)
        else:
            neutral.append(ons)

    def cat(xs):
        return np.sort(np.concatenate(xs)) if xs else np.array([])

    return {
        "Reward":  cat(reward),
        "Neutral": cat(neutral),
        "Punish":  cat(punish),
    }


def evs_sr_decision_partner(ev_dir: Path) -> Dict[str, np.ndarray]:
    # Optional extra: decision locked by partner (2 conditions)
    return {
        "Computer": np.sort(load_ev_onsets(ev_dir / "_guess_computer.txt")),
        "Stranger": np.sort(load_ev_onsets(ev_dir / "_guess_face.txt")),
    }


def evs_sr_decision_by_future_outcome(ev_dir: Path) -> Dict[str, np.ndarray]:
    """
    Optional extra: decision locked, split by subsequent outcome valence (3 conditions)
    Uses the same nearest-future-outcome mapping logic you used before.
    """
    d_all = np.sort(np.concatenate([
        load_ev_onsets(ev_dir / "_guess_computer.txt"),
        load_ev_onsets(ev_dir / "_guess_face.txt"),
    ]))
    out_map = {}
    for p in ev_dir.glob("_outcome_*.txt"):
        ons = load_ev_onsets(p)
        nm = p.name.lower()
        if "reward" in nm:
            val = "Reward"
        elif "punish" in nm:
            val = "Punish"
        else:
            val = "Neutral"
        for o in ons:
            out_map[float(o)] = val

    if d_all.size == 0 or not out_map:
        return {"Reward": np.array([]), "Neutral": np.array([]), "Punish": np.array([])}

    all_outs = np.sort(np.array(list(out_map.keys())))
    bins = {"Reward": [], "Neutral": [], "Punish": []}
    for d in d_all:
        idx = np.searchsorted(all_outs, d)
        if idx < len(all_outs):
            o = all_outs[idx]
            if (o - d) < SR_MAX_DEC_TO_OUT_GAP_S:
                bins[out_map[o]].append(float(d))

    return {k: np.array(sorted(v)) for k, v in bins.items()}


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


# ----------------------------- plotting -----------------------------

def plot_curves_with_sem(out_png: Path, title: str, t: np.ndarray,
                         curves: Dict[str, Tuple[np.ndarray, np.ndarray, int]],
                         vline_label: str = "Onset") -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    for name in sorted(curves.keys()):
        mu, se, n = curves[name]
        if mu.size == 0:
            continue
        ax.plot(t, mu, label=f"{name} (n={n})")
        ax.fill_between(t, mu - se, mu + se, alpha=0.15)

    ax.axvline(0.0, color="k", linestyle=":", linewidth=1)
    ax.set_title(title)
    ax.set_xlabel("Time from event onset (s)")
    ax.set_ylabel("PSC (epoch-baseline referenced)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, ncol=2)

    set_data_driven_ylim(ax)

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)


def plot_compare_baselines(out_png: Path, title: str, t: np.ndarray,
                           curves_base2: Dict[str, Tuple[np.ndarray, np.ndarray, int]],
                           curves_base3: Dict[str, Tuple[np.ndarray, np.ndarray, int]]) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, (lbl, curves) in zip(axes, [("Baseline = 2 TR", curves_base2), ("Baseline = 3 TR", curves_base3)]):
        for name in sorted(curves.keys()):
            mu, se, n = curves[name]
            if mu.size == 0:
                continue
            ax.plot(t, mu, label=f"{name} (n={n})")
            ax.fill_between(t, mu - se, mu + se, alpha=0.15)
        ax.axvline(0.0, color="k", linestyle=":", linewidth=1)
        ax.set_title(lbl)
        ax.set_xlabel("Time from event onset (s)")
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel("PSC (epoch-baseline referenced)")
    axes[1].legend(fontsize=8, ncol=1)

    # data-driven shared y-limits across both panels
    ys = []
    for ax in axes:
        for line in ax.get_lines():
            y = np.asarray(line.get_ydata())
            ys.append(y[np.isfinite(y)])
    if ys:
        yall = np.concatenate([y for y in ys if y.size], axis=0)
        if yall.size:
            lo, hi = float(np.min(yall)), float(np.max(yall))
            pad = (hi - lo) * 0.1 if hi != lo else (0.25 if lo == 0 else abs(lo) * 0.1)
            axes[0].set_ylim(lo - pad, hi + pad)

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

    # Define plot groups
    plot_groups: Dict[str, Dict[str, np.ndarray]] = {}

    if task == "mid":
        plot_groups["mid_anticipation_primary"] = evs_mid_primary(ev_dir)
        plot_groups["mid_feedback_extra"] = evs_mid_feedback(ev_dir)
    else:
        plot_groups["sr_outcome_primary"] = evs_sr_outcome_primary(ev_dir)
        plot_groups["sr_decision_partner_extra"] = evs_sr_decision_partner(ev_dir)
        plot_groups["sr_decision_by_future_outcome_extra"] = evs_sr_decision_by_future_outcome(ev_dir)

    out: List[RunCurve] = []

    for baseline_ntr in BASELINE_NTRS:
        base_window = (-baseline_ntr * tr, 0.0)

        for plot_group, condmap in plot_groups.items():
            # Epoch each condition
            for cond, onsets in condmap.items():
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
            "time_s","mean","sem","n_events","feat_dir"
        ])
        for c in curves:
            if c.mean.size == 0:
                continue
            for ti, mi, si in zip(c.t, c.mean, c.sem):
                w.writerow([
                    c.sub, c.ses, c.run, c.task, c.echo, c.cnfds,
                    c.baseline_ntr, c.plot_group, c.condition,
                    f"{ti:.4f}", f"{mi:.8f}", f"{si:.8f}", c.n_events, c.feat_dir
                ])


def main() -> None:
    pathlist = find_pathlist_file()
    feat_paths = read_feat_paths(pathlist)
    if not feat_paths:
        raise RuntimeError(f"Path list is empty: {pathlist}")

    # Infer project root from first path
    root = infer_project_root_from_feat(feat_paths[0])
    mask_path = root / MASK_REL
    if not mask_path.exists():
        raise FileNotFoundError(f"ROI mask not found: {mask_path}")

    fsl_deriv = root / "derivatives" / "fsl"
    ev_base = fsl_deriv / "EVFiles"

    out_base = root / "derivatives" / "extractions" / "timecourses-unified-noargs"
    out_tables = out_base / "tables"
    out_runplots = out_base / "run_plots"
    out_session_qc = out_base / "session_qc"
    out_subject_qc = out_base / "subject_qc"

    mask_bool = load_mask(mask_path)

    # ---- Process all runs ----
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

    # ---- Write long table ----
    out_csv = out_tables / "run_curves_long.csv"
    write_run_curves_long_csv(out_csv, all_curves)
    print(f"[INFO] Wrote: {out_csv}")

    # ---- RUN-LEVEL PLOTS ----
    # Group: (sub,ses,run,task,echo,cnfds,baseline_ntr,plot_group) -> condition -> curve
    run_groups: Dict[Tuple[str,str,str,str,str,str,int,str], Dict[str, RunCurve]] = {}
    for c in all_curves:
        k = (c.sub, c.ses, c.run, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group)
        run_groups.setdefault(k, {})[c.condition] = c

    for (sub, ses, run, task, echo, cnfds, bntr, plot_group), condmap in run_groups.items():
        # consistent time axis
        anyc = next(iter(condmap.values()))
        t = anyc.t
        curves_dict = {cond: (obj.mean, obj.sem, obj.n_events) for cond, obj in condmap.items()}
        out_png = out_runplots / f"baseline-{bntr}TR" / f"task-{task}" / f"sub-{sub}" / f"ses-{ses}" / f"run-{run}" / echo / f"cnfds-{cnfds}" / f"{plot_group}.png"
        title = f"{task} {plot_group} — sub-{sub} ses-{ses} run-{run} {echo} cnfds-{cnfds} (baseline={bntr}TR)"
        plot_curves_with_sem(out_png, title, t, curves_dict)

    print(f"[INFO] Wrote run plots under: {out_runplots}")

    # ---- SESSION-LEVEL QC SUMMARY (single line + SEM across runs) ----
    # For each session: aggregate across run-level mean curves (not across trials)
    # Group: (sub,ses,task,echo,cnfds,baseline_ntr,plot_group,condition) -> list of run-mean curves
    sess_bucket: Dict[Tuple[str,str,str,str,str,int,str,str], List[np.ndarray]] = {}
    sess_taxis: Dict[Tuple[str,str,str,str,str,int,str], np.ndarray] = {}

    for c in all_curves:
        # run mean curve is c.mean
        k = (c.sub, c.ses, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group, c.condition)
        sess_bucket.setdefault(k, []).append(c.mean)
        k2 = (c.sub, c.ses, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group)
        sess_taxis[k2] = c.t

    # aggregate per-session plot
    # Group plot: (sub,ses,task,echo,cnfds,baseline_ntr,plot_group) -> condition -> (mu, sem, n_runs)
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

    for (sub, ses, task, echo, cnfds, bntr, plot_group), condmap in sess_plots.items():
        t = sess_taxis[(sub, ses, task, echo, cnfds, bntr, plot_group)]
        out_png = out_session_qc / f"baseline-{bntr}TR" / f"task-{task}" / f"sub-{sub}" / f"ses-{ses}" / echo / f"cnfds-{cnfds}" / f"{plot_group}_summary.png"
        title = f"SESSION QC: {task} {plot_group} — sub-{sub} ses-{ses} {echo} cnfds-{cnfds} (baseline={bntr}TR)"
        plot_curves_with_sem(out_png, title, t, condmap)

    # baseline comparison plots for *primary* QC groups only
    primary_groups = set(["mid_anticipation_primary", "sr_outcome_primary"])
    # Group key without baseline: (sub,ses,task,echo,cnfds,plot_group)
    sess_primary_compare: Dict[Tuple[str,str,str,str,str,str], Dict[int, Dict[str, Tuple[np.ndarray,np.ndarray,int]]]] = {}
    sess_primary_t: Dict[Tuple[str,str,str,str,str,str], np.ndarray] = {}

    for (sub, ses, task, echo, cnfds, bntr, plot_group), condmap in sess_plots.items():
        if plot_group not in primary_groups:
            continue
        k = (sub, ses, task, echo, cnfds, plot_group)
        sess_primary_compare.setdefault(k, {})[bntr] = condmap
        sess_primary_t[k] = sess_taxis[(sub, ses, task, echo, cnfds, bntr, plot_group)]

    for k, by_base in sess_primary_compare.items():
        if 2 not in by_base or 3 not in by_base:
            continue
        sub, ses, task, echo, cnfds, plot_group = k
        t = sess_primary_t[k]
        out_png = out_session_qc / "compare_baselines" / f"task-{task}" / f"sub-{sub}" / f"ses-{ses}" / echo / f"cnfds-{cnfds}" / f"{plot_group}_baseline2vs3.png"
        title = f"SESSION QC (baseline compare): {task} {plot_group} — sub-{sub} ses-{ses} {echo} cnfds-{cnfds}"
        plot_compare_baselines(out_png, title, t, by_base[2], by_base[3])

    print(f"[INFO] Wrote session QC under: {out_session_qc}")

    # ---- SUBJECT-LEVEL QC SUMMARY (single line + SEM across runs across sessions) ----
    subj_bucket: Dict[Tuple[str,str,str,str,int,str,str], List[np.ndarray]] = {}
    subj_taxis: Dict[Tuple[str,str,str,str,int,str], np.ndarray] = {}

    for c in all_curves:
        k = (c.sub, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group, c.condition)
        subj_bucket.setdefault(k, []).append(c.mean)
        k2 = (c.sub, c.task, c.echo, c.cnfds, c.baseline_ntr, c.plot_group)
        subj_taxis[k2] = c.t

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

    for (sub, task, echo, cnfds, bntr, plot_group), condmap in subj_plots.items():
        t = subj_taxis[(sub, task, echo, cnfds, bntr, plot_group)]
        out_png = out_subject_qc / f"baseline-{bntr}TR" / f"task-{task}" / f"sub-{sub}" / echo / f"cnfds-{cnfds}" / f"{plot_group}_summary.png"
        title = f"SUBJECT QC: {task} {plot_group} — sub-{sub} {echo} cnfds-{cnfds} (baseline={bntr}TR)"
        plot_curves_with_sem(out_png, title, t, condmap)

    # baseline comparison plots for *primary* QC groups only
    subj_primary_compare: Dict[Tuple[str,str,str,str,str], Dict[int, Dict[str, Tuple[np.ndarray,np.ndarray,int]]]] = {}
    subj_primary_t: Dict[Tuple[str,str,str,str,str], np.ndarray] = {}

    for (sub, task, echo, cnfds, bntr, plot_group), condmap in subj_plots.items():
        if plot_group not in primary_groups:
            continue
        k = (sub, task, echo, cnfds, plot_group)
        subj_primary_compare.setdefault(k, {})[bntr] = condmap
        subj_primary_t[k] = subj_taxis[(sub, task, echo, cnfds, bntr, plot_group)]

    for k, by_base in subj_primary_compare.items():
        if 2 not in by_base or 3 not in by_base:
            continue
        sub, task, echo, cnfds, plot_group = k
        t = subj_primary_t[k]
        out_png = out_subject_qc / "compare_baselines" / f"task-{task}" / f"sub-{sub}" / echo / f"cnfds-{cnfds}" / f"{plot_group}_baseline2vs3.png"
        title = f"SUBJECT QC (baseline compare): {task} {plot_group} — sub-{sub} {echo} cnfds-{cnfds}"
        plot_compare_baselines(out_png, title, t, by_base[2], by_base[3])

    print(f"[INFO] Wrote subject QC under: {out_subject_qc}")
    print("[INFO] Done.")


if __name__ == "__main__":
    main()
