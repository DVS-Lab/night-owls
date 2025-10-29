#!/usr/bin/env python3
"""
SR task (sharedreward) — VS outcome‑locked time‑courses and discrete 4th‑TR extraction
-------------------------------------------------------------------------------------
Assumptions
- Run from <rootdir>/code (this file lives in that directory).
- All data are already in MNI152NLin6Asym space.
- Input FEAT directories are listed one per line in ./feat_paths-sr.txt (absolute
  paths), exactly as provided by the user. The script ignores blank lines and
  lines starting with '#'.
- EV files live at:
    <rootdir>/derivatives/fsl/EVFiles/sub-XXX/ses-YY/sharedreward/run-ZZ/
  with the following filenames (we pool across agent: computer/stranger):
    _outcome_computer_reward.txt
    _outcome_computer_neutral.txt
    _outcome_computer_punish.txt
    _outcome_stranger_reward.txt
    _outcome_stranger_neutral.txt
    _outcome_stranger_punish.txt
- Ventral Striatum (VS) mask lives at:
    <rootdir>/masks/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz

What it does
- Extracts the VS mean time series from each FEAT's filtered_func_data.nii.gz.
- Computes two normalized versions per run:
    • PSC   : ((ts / ts.mean()) - 1) * 100
    • Z     : (ts - ts.mean()) / ts.std(ddof=1)
- Builds event‑locked peri‑outcome windows with linear interpolation to handle
  fractional onsets.
- Produces per‑run plots and subject‑level aggregated plots (split by echo)
  for PSC and Z, showing: Reward, Neutral, Punish.
- Writes **discrete 4th‑TR** summaries (no interpolation) for Reward/Neutral/Punish
  at ~4.845–6.46 s after outcome onset.

Outputs
  <rootdir>/derivatives/extractions/timecourses_sr/
    ├── runs/  (per FEAT plots + per‑run TP TSV)
    └── subjects/  (per subject, echo‑split plots)
  <rootdir>/derivatives/extractions/summary_sr_at_4thTR.tsv  (and .csv)

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
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses_sr"
SUMMARY_DIR= ROOT_DIR / "derivatives" / "extractions"
EV_BASE    = FSL_DERIV / "EVFiles"

VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

# Input list of FEAT directories for SR (one per line, absolute).
FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-sr.txt"

# --------------------------- Analysis parameters ------------------------------
TASK         = "sharedreward"
TR_HARDCODE  = 1.615  # sec; also checked in header and warned if different
TMIN         = -4.0   # sec relative to outcome onset (pre‑event baseline)
TMAX         = 16.0   # sec after outcome onset
K_AFTER      = 3      # 4th TR = 3 TRs after onset
VERT_LINES   = [0.0, 4.0, 6.0]

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
    echo: str
    time_axis: np.ndarray
    out_psc: Dict[str, Tuple[np.ndarray, np.ndarray, int]]  # Reward/Neutral/Punish
    out_z:   Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    # discrete 4th‑TR summaries
    tp_means_psc: Dict[str, float | None]
    tp_means_z:   Dict[str, float | None]
    tp_counts:    Dict[str, int]

# ------------------------------ Helper functions ------------------------------

def parse_sub_ses_run_from_feat(feat: Path) -> Tuple[str, str, str]:
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
    vs_ts = np.nanmean(vox_ts, axis=0)
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
    if onsets.size == 0:
        return np.empty((0,), dtype=int)
    idx = np.floor((onsets + K_AFTER * tr) / tr).astype(int)
    return idx[(idx >= 0) & (idx < n_vols)]


def plot_three(time_axis: np.ndarray,
               A: Tuple[np.ndarray, np.ndarray, int], labelA: str,
               B: Tuple[np.ndarray, np.ndarray, int], labelB: str,
               C: Tuple[np.ndarray, np.ndarray, int], labelC: str,
               title: str, ylabel: str, out_png: Path) -> None:
    mA, sA, nA = A; mB, sB, nB = B; mC, sC, nC = C
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.plot(time_axis, mA, label=f"{labelA} (n={nA})"); ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.25)
    ax.plot(time_axis, mB, label=f"{labelB} (n={nB})"); ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.25)
    ax.plot(time_axis, mC, label=f"{labelC} (n={nC})"); ax.fill_between(time_axis, mC - sC, mC + sC, alpha=0.25)
    for v in VERT_LINES:
        ax.axvline(v, ls=":", lw=1)
    ax.set_xlabel("Time from outcome (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(out_png); plt.close(fig)

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

    # Load EVs (pool across agent)
    ev_dir = get_ev_dir(sub, ses, run)
    rew = np.sort(np.concatenate([
        load_ev(ev_dir / "_outcome_computer_reward.txt"),
        load_ev(ev_dir / "_outcome_stranger_reward.txt"),
    ]))
    neu = np.sort(np.concatenate([
        load_ev(ev_dir / "_outcome_computer_neutral.txt"),
        load_ev(ev_dir / "_outcome_stranger_neutral.txt"),
    ]))
    pun = np.sort(np.concatenate([
        load_ev(ev_dir / "_outcome_computer_punish.txt"),
        load_ev(ev_dir / "_outcome_stranger_punish.txt"),
    ]))

    t_axis = build_time_axis(TMIN, TMAX, tr)

    # Windows for PSC
    R_psc_w = sample_windows(ts_psc, rew, tr, TMIN, TMAX)
    N_psc_w = sample_windows(ts_psc, neu, tr, TMIN, TMAX)
    P_psc_w = sample_windows(ts_psc, pun, tr, TMIN, TMAX)
    # Windows for Z
    R_z_w   = sample_windows(ts_z,   rew, tr, TMIN, TMAX)
    N_z_w   = sample_windows(ts_z,   neu, tr, TMIN, TMAX)
    P_z_w   = sample_windows(ts_z,   pun, tr, TMIN, TMAX)

    out_psc = {
        "Reward":  mean_and_sem(R_psc_w),
        "Neutral": mean_and_sem(N_psc_w),
        "Punish":  mean_and_sem(P_psc_w),
    }
    out_z = {
        "Reward":  mean_and_sem(R_z_w),
        "Neutral": mean_and_sem(N_z_w),
        "Punish":  mean_and_sem(P_z_w),
    }

    # Discrete 4th‑TR (no interpolation)
    def val_at_tp(ts: np.ndarray, onsets: np.ndarray) -> np.ndarray:
        idx = fourth_tr_indices(onsets, tr, T)
        return ts[idx]

    tp_means_psc = {
        "REWARD":  float(np.nanmean(val_at_tp(ts_psc, rew))) if rew.size else None,
        "NEUTRAL": float(np.nanmean(val_at_tp(ts_psc, neu))) if neu.size else None,
        "PUNISH":  float(np.nanmean(val_at_tp(ts_psc, pun))) if pun.size else None,
    }
    tp_means_z = {
        "REWARD":  float(np.nanmean(val_at_tp(ts_z, rew))) if rew.size else None,
        "NEUTRAL": float(np.nanmean(val_at_tp(ts_z, neu))) if neu.size else None,
        "PUNISH":  float(np.nanmean(val_at_tp(ts_z, pun))) if pun.size else None,
    }
    tp_counts = {
        "REWARD":  int(fourth_tr_indices(rew, tr, T).size),
        "NEUTRAL": int(fourth_tr_indices(neu, tr, T).size),
        "PUNISH":  int(fourth_tr_indices(pun, tr, T).size),
    }

    return RunResult(
        feat_path=feat,
        sub=sub, ses=ses, run=run, echo=echo,
        time_axis=t_axis,
        out_psc=out_psc, out_z=out_z,
        tp_means_psc=tp_means_psc,
        tp_means_z=tp_means_z,
        tp_counts=tp_counts,
    )


def save_run_outputs(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo

    # Plots
    plot_three(
        res.time_axis,
        res.out_psc["Reward"],  "Reward",
        res.out_psc["Neutral"], "Neutral",
        res.out_psc["Punish"],  "Punish",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (Outcome, PSC)",
        ylabel="% signal change (PSC)",
        out_png=run_out / "outcome_psc.png",
    )
    plot_three(
        res.time_axis,
        res.out_z["Reward"],  "Reward",
        res.out_z["Neutral"], "Neutral",
        res.out_z["Punish"],  "Punish",
        title=f"VS — sub {res.sub} ses {res.ses} run {res.run} [{res.echo}] (Outcome, Z)",
        ylabel="Z (SD units)",
        out_png=run_out / "outcome_z.png",
    )

    # Summary text
    summary = (
        f"Run: sub-{res.sub} ses-{res.ses} run-{res.run} echo={res.echo}\n"
        f"Outcome n — Reward={res.out_psc['Reward'][2]}, Neutral={res.out_psc['Neutral'][2]}, Punish={res.out_psc['Punish'][2]}\n"
    )
    run_out.mkdir(parents=True, exist_ok=True)
    (run_out / "summary.txt").write_text(summary)

    # Per‑run discrete TSV
    header = [
        "sub","ses","run","echo",
        "REWARD_PSC","NEUTRAL_PSC","PUNISH_PSC",
        "REWARD_Z","NEUTRAL_Z","PUNISH_Z",
        "N_REWARD","N_NEUTRAL","N_PUNISH",
    ]
    row = [
        res.sub, res.ses, res.run, res.echo,
        _fmt(res.tp_means_psc.get("REWARD")),
        _fmt(res.tp_means_psc.get("NEUTRAL")),
        _fmt(res.tp_means_psc.get("PUNISH")),
        _fmt(res.tp_means_z.get("REWARD")),
        _fmt(res.tp_means_z.get("NEUTRAL")),
        _fmt(res.tp_means_z.get("PUNISH")),
        str(res.tp_counts.get("REWARD",0)),
        str(res.tp_counts.get("NEUTRAL",0)),
        str(res.tp_counts.get("PUNISH",0)),
    ]
    tp_path = run_out / "tp_4thTR.tsv"
    with open(tp_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")


# -------------------------- Subject‑level aggregation --------------------------

def aggregate_subject(results: List[RunResult], subject: str) -> None:
    time_axis = None
    pool: Dict[str, Dict[str, Dict[str, List[Tuple[np.ndarray, int]]]]] = {}

    for res in results:
        if res.sub != subject:
            continue
        if time_axis is None:
            time_axis = res.time_axis
        echo = res.echo
        pool.setdefault(echo, {
            "PSC": {"Reward": [], "Neutral": [], "Punish": []},
            "Z":   {"Reward": [], "Neutral": [], "Punish": []},
        })
        for k in ("Reward","Neutral","Punish"):
            pool[echo]["PSC"][k].append((res.out_psc[k][0], res.out_psc[k][2]))
            pool[echo]["Z"][k].append((res.out_z[k][0],   res.out_z[k][2]))

    if time_axis is None:
        print(f"[INFO] No runs aggregated for subject {subject}.")
        return

    subj_out = OUT_TC_DIR / "subjects" / f"sub-{subject}"

    def weighted_msem(items: List[Tuple[np.ndarray, int]]) -> Tuple[np.ndarray, np.ndarray, int]:
        if not items:
            return np.array([]), np.array([]), 0
        means = np.vstack([m for (m, n) in items])
        sem = np.nanstd(means, axis=0, ddof=1) / np.sqrt(max(len(items),1))
        n_total = int(np.sum([max(n,0) for (m,n) in items]))
        mean = np.nanmean(means, axis=0)
        return mean, sem, n_total

    for echo, fam in pool.items():
        psc_R = weighted_msem(fam["PSC"]["Reward"])
        psc_N = weighted_msem(fam["PSC"]["Neutral"])
        psc_P = weighted_msem(fam["PSC"]["Punish"])
        plot_three(time_axis, psc_R, "Reward", psc_N, "Neutral", psc_P, "Punish",
                   title=f"VS — subject {subject} [{echo}] (Outcome, PSC)",
                   ylabel="% signal change (PSC)",
                   out_png=subj_out / echo / "outcome_psc.png")

        z_R = weighted_msem(fam["Z"]["Reward"])
        z_N = weighted_msem(fam["Z"]["Neutral"])
        z_P = weighted_msem(fam["Z"]["Punish"])
        plot_three(time_axis, z_R, "Reward", z_N, "Neutral", z_P, "Punish",
                   title=f"VS — subject {subject} [{echo}] (Outcome, Z)",
                   ylabel="Z (SD units)",
                   out_png=subj_out / echo / "outcome_z.png")

        summary = (
            f"Subject {subject} — Echo: {echo}\n"
            f"Outcome n — Reward={psc_R[2]}, Neutral={psc_N[2]}, Punish={psc_P[2]}\n"
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
        print("Create feat_paths-sr.txt with one FEAT directory per line (absolute paths).")
        return

    feat_paths: List[Path] = []
    for line in FEAT_LIST_PATH.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        feat_paths.append(Path(s))

    if not feat_paths:
        print("No FEATs to process (feat_paths-sr.txt empty). Exiting.")
        return

    results: List[RunResult] = []
    rows_for_tp: List[List[str]] = []

    for feat in feat_paths:
        res = process_one_feat(feat)
        if res is None:
            continue
        save_run_outputs(res)
        results.append(res)
        rows_for_tp.append([
            res.sub, res.ses, res.run, res.echo,
            _fmt(res.tp_means_psc.get("REWARD")),
            _fmt(res.tp_means_psc.get("NEUTRAL")),
            _fmt(res.tp_means_psc.get("PUNISH")),
            _fmt(res.tp_means_z.get("REWARD")),
            _fmt(res.tp_means_z.get("NEUTRAL")),
            _fmt(res.tp_means_z.get("PUNISH")),
            str(res.tp_counts.get("REWARD",0)),
            str(res.tp_counts.get("NEUTRAL",0)),
            str(res.tp_counts.get("PUNISH",0)),
        ])

    # Aggregate by subject (split by echo)
    subjects = sorted({r.sub for r in results})
    for sub in subjects:
        aggregate_subject(results, sub)

    # Write combined spreadsheet (TSV + CSV)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    header = [
        "sub","ses","run","echo",
        "REWARD_PSC","NEUTRAL_PSC","PUNISH_PSC",
        "REWARD_Z","NEUTRAL_Z","PUNISH_Z",
        "N_REWARD","N_NEUTRAL","N_PUNISH",
    ]
    tsv_path = SUMMARY_DIR / "summary_sr_at_4thTR.tsv"
    csv_path = SUMMARY_DIR / "summary_sr_at_4thTR.csv"
    with open(tsv_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in rows_for_tp:
            f.write("\t".join(row) + "\n")
    with open(csv_path, 'w') as f:
        f.write(','.join(header) + "\n")
        for row in rows_for_tp:
            f.write(','.join(row) + "\n")

    print(f"Done. Outputs in: {OUT_TC_DIR}\n  - Discrete summaries: {tsv_path} and {csv_path}")


if __name__ == "__main__":
    main()
