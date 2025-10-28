#!/usr/bin/env python3
"""
vs_timecourse.py — VS‑only peri‑stimulus extraction (PSC + Z), path‑list driven

Run from <rootdir>/code. All paths are resolved relative to this script.

Why this version?
- Avoids missed runs by reading an explicit **feat_paths.txt** list of FEAT directories.
- Computes **two metrics** from the VS ROI:
  1) % signal change (PSC) relative to each run's mean
  2) Z‑scored time series across the entire run (per‑run mean=0, SD=1)
- Adds **error bands (SEM)** for per‑run plots and for **subject‑level** aggregates (pooling trials across runs/sessions).
- Uses a **fixed TR = 1.615 s** for windowing, with a header sanity check (warn if header TR differs by >5 ms).
- Extends the pre‑stimulus baseline to **−4 s** (configurable).
- Hard‑codes EV file names and cue‑onset alignment.

Inputs
- VS mask: masks/space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz (must be on the same grid as functional data)
- FEAT list file: <rootdir>/code/feat_paths.txt (one absolute FEAT dir per line; lines may include a trailing '+' before .feat)
- EVs (hard‑coded):
  derivatives/fsl/EVFiles/sub-<sub>/ses-<ses>/<task>/run-<run>/_anticipation_reward.txt
  derivatives/fsl/EVFiles/sub-<sub>/ses-<ses>/<task>/run-<run>/_anticipation_neutral.txt

Outputs
  derivatives/extractions/vs_timecourse/
   ├── <feat_basename>/csv/vs_t-4_to_16_psc.csv, vs_t-4_to_16_z.csv
   ├── <feat_basename>/figs/vs_t-4_to_16_psc.png, vs_t-4_to_16_z.png
   ├── <feat_basename>/masks/vs_func_assumed_mni.nii.gz
   ├── <feat_basename>/ts/vs_ts.txt (raw), vs_psc.txt, vs_z.txt
   └── _subject_averages/
        subject-<sub>_vs_timecourse_psc.csv/.png
        subject-<sub>_vs_timecourse_z.csv/.png
        subject-<sub>_summary.txt

Notes
- Sessions 01–12 are allowed; a hard‑coded SKIP map filters known bad sessions.
- FEAT list parsing is robust to names like ...cnfds-fmriprep+.feat (the '+' is fine).
- For windowing, events whose windows would extend beyond the run boundaries are dropped (counted in logs).
"""
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import csv
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# -------------------------
# Fixed locations (relative to this script)
# -------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
DERIV_FSL  = ROOT_DIR / "derivatives" / "fsl"
EVFILES    = DERIV_FSL / "EVFiles"
MASKS_DIR  = ROOT_DIR / "masks"
OUT_ROOT   = ROOT_DIR / "derivatives" / "extractions" / "vs_timecourse"
FEAT_LIST  = SCRIPT_DIR / "feat_paths.txt"  # one absolute FEAT path per line

VS_MASK = MASKS_DIR / "space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

# -------------------------
# HARD‑CODED SETTINGS
# -------------------------
TASK       = "mid"            # EV subfolder
ALIGN      = "onset"          # lock to cue onset
TR_FIXED   = 1.615            # seconds (used for indexing)
TR_TOL     = 0.005            # warn if |TR_header − TR_FIXED| > TR_TOL
TMIN       = -4.0             # seconds
TMAX       = 16.0             # seconds
PEAK_LAGS  = [6.0, 4.0]       # vertical reference lines

SESSIONS    = [f"{i:02d}" for i in range(1, 13)]  # 01..12
SKIP_MAP: Dict[str, set] = {
    "101": {"04", "05", "12"},
    "103": {"12"},
}

# tokens seen in FEAT names (for sanity checks only)
SPACE_TOKS   = ["space-mni", "space-MNI152NLin6Asym"]
ECHO_TOKS    = ["single-echo", "multi-echo"]
CONFOUND_TOKS= ["cnfds-fmriprep"]

# -------------------------
# Small helpers
# -------------------------

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def load_nifti(path: Path):
    img = nib.load(str(path))
    return img.get_fdata(), img.affine, img.header


def header_tr(path: Path) -> float:
    """Read TR from NIfTI header; returns seconds (float)."""
    hdr = nib.load(str(path)).header
    z = hdr.get_zooms()
    if len(z) < 4:
        raise ValueError("Functional NIfTI missing time dimension (zoom[3])")
    return float(z[3])


def psc_vs_runmean(ts: np.ndarray) -> np.ndarray:
    mu = float(np.mean(ts))
    if mu == 0:
        return np.zeros_like(ts)
    return 100.0 * (ts - mu) / mu


def zscore_run(ts: np.ndarray) -> np.ndarray:
    mu = float(np.mean(ts))
    sd = float(np.std(ts, ddof=1))
    if sd == 0:
        return np.zeros_like(ts)
    return (ts - mu) / sd


def load_ev_times(ev_path: Path, align: str) -> np.ndarray:
    arr = np.loadtxt(str(ev_path))
    if arr.ndim == 1:
        arr = arr[None, :]
    on = arr[:, 0].astype(float)
    du = arr[:, 1].astype(float)
    wt = arr[:, 2].astype(float)
    keep = wt > 0
    on, du = on[keep], du[keep]
    if align == "onset":
        t = on
    elif align == "offset":
        t = on + du
    elif align == "center":
        t = on + 0.5 * du
    else:
        raise ValueError("align must be onset|offset|center")
    return t


def extract_windows(series: np.ndarray, onsets_s: np.ndarray, tr_s: float, tmin_s: float, tmax_s: float) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """Return (times, windows, n_kept, n_dropped).
    windows shape: [n_kept, n_timepoints]. Trials whose window would cross run edges are dropped.
    """
    n_win = int(np.floor((tmax_s - tmin_s) / tr_s)) + 1
    times = np.arange(n_win) * tr_s + tmin_s
    T = series.shape[0]
    kept = []
    dropped = 0
    for onset in onsets_s:
        idx0  = int(round(onset / tr_s))
        start = idx0 + int(round(tmin_s / tr_s))
        end   = start + n_win
        if start < 0 or end > T:
            dropped += 1
            continue
        kept.append(series[start:end])
    if len(kept) == 0:
        return times, np.empty((0, n_win)), 0, dropped
    return times, np.vstack(kept), len(kept), dropped


def stats_from_windows(wins: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int, np.ndarray, np.ndarray]:
    """Mean, SEM, n, sum, sumsq along axis=0."""
    if wins.size == 0:
        return (np.full(0, np.nan), np.full(0, np.nan), 0, np.array([]), np.array([]))
    n = wins.shape[0]
    s = wins.sum(axis=0)
    ss = (wins * wins).sum(axis=0)
    mean = s / n
    if n > 1:
        var = (ss - (s*s)/n) / (n - 1)
        sem = np.sqrt(np.maximum(var, 0)) / np.sqrt(n)
    else:
        sem = np.full_like(mean, np.nan)
    return mean, sem, n, s, ss


def parse_feat_tokens(feat_dir: Path) -> Tuple[str, str, str, str]:
    """Extract (sub, ses, task, run) from a FEAT directory name without regex backslashes."""
    name = feat_dir.name.replace('.feat','')
    parts = name.split('_')
    tokens: Dict[str, str] = {}
    for p in parts:
        if p.startswith('sub-'): tokens['sub'] = p[4:]
        elif p.startswith('ses-'): tokens['ses'] = p[4:]
        elif p.startswith('task-'): tokens['task'] = p[5:]
        elif p.startswith('run-'): tokens['run'] = p[4:]
    need = ['sub','ses','task','run']
    if not all(k in tokens for k in need):
        raise ValueError(f"Cannot parse sub/ses/task/run from FEAT name: {feat_dir.name}")
    return tokens['sub'], tokens['ses'], tokens['task'], tokens['run']

# -------------------------
# Plotting
# -------------------------

def plot_with_band(times: np.ndarray,
                   mean_a: np.ndarray, sem_a: np.ndarray,
                   mean_b: np.ndarray, sem_b: np.ndarray,
                   label_a: str, label_b: str,
                   title: str, out_png: Path):
    plt.figure(figsize=(7.0, 4.2))
    if mean_a.size:
        plt.plot(times, mean_a, label=label_a)
        if sem_a is not None and sem_a.size:
            plt.fill_between(times, mean_a - sem_a, mean_a + sem_a, alpha=0.25, linewidth=0)
    if mean_b.size:
        plt.plot(times, mean_b, label=label_b)
        if sem_b is not None and sem_b.size:
            plt.fill_between(times, mean_b - sem_b, mean_b + sem_b, alpha=0.25, linewidth=0)
    plt.axvline(0, linestyle='--', linewidth=1)
    for v in PEAK_LAGS:
        plt.axvline(v, linestyle=':', linewidth=1)
    plt.xlabel('Time from cue onset (s)')
    plt.ylabel('% signal change (PSC)' if 'PSC' in title else 'Z score (run‑wise)')
    plt.title(title)
    plt.legend(frameon=False)
    plt.tight_layout(); plt.savefig(out_png, dpi=150); plt.close()

# -------------------------
# Per‑run computation
# -------------------------

def process_feat(feat_dir: Path, vs_mask: Path) -> Dict[str, object]:
    """Process one FEAT; return dict with per‑run stats and aggregation pieces.
    Keys: 'sub','ses','run','feat_base','t','psc_R_mean','psc_R_sem','psc_N_mean','psc_N_sem',
          'z_R_mean','z_R_sem','z_N_mean','z_N_sem','nR','nN',
          'psc_R_sum','psc_R_sumsq','psc_N_sum','psc_N_sumsq',
          'z_R_sum','z_R_sumsq','z_N_sum','z_N_sumsq'
    """
    func = feat_dir / 'filtered_func_data.nii.gz'
    if not func.exists():
        raise FileNotFoundError(f"Missing {func}")

    sub, ses, task, run = parse_feat_tokens(feat_dir)
    if task != TASK:
        raise ValueError(f"FEAT task token '{task}' != expected '{TASK}' in {feat_dir}")
    if ses in SKIP_MAP.get(sub, set()):
        raise RuntimeError(f"Session marked to skip: sub-{sub} ses-{ses}")

    feat_base = feat_dir.name.replace('.feat','')
    outdir = OUT_ROOT / feat_base
    masks_dir = outdir / 'masks'; ts_dir = outdir / 'ts'; csv_dir = outdir / 'csv'; fig_dir = outdir / 'figs'
    for d in (masks_dir, ts_dir, csv_dir, fig_dir):
        ensure_dir(d)

    # Grid check
    func_data, func_aff, func_hdr = load_nifti(func)
    vs_data,   vs_aff,   vs_hdr   = load_nifti(vs_mask)
    if func_data.shape[:3] != vs_data.shape[:3] or not np.allclose(func_aff, vs_aff, atol=1e-3):
        raise RuntimeError(f"Grid mismatch: {feat_dir.name} vs VS mask; pre‑align before running.")

    # Save assumed‑MNI mask copy (for provenance)
    nib.Nifti1Image((vs_data>0).astype(np.uint8), vs_aff, vs_hdr).to_filename(str(masks_dir / 'vs_func_assumed_mni.nii.gz'))

    # Extract TS and metrics
    mask = vs_data > 0
    T = func_data.shape[3]
    ts_raw = func_data[mask].reshape((-1, T)).mean(axis=0)
    np.savetxt(ts_dir / 'vs_ts.txt', ts_raw[np.newaxis, :], fmt='%.6f')

    # TR check (we still **use TR_FIXED** for indexing)
    tr_hdr = header_tr(func)
    if abs(tr_hdr - TR_FIXED) > TR_TOL:
        print(f"[WARN] {feat_dir.name}: header TR={tr_hdr:.6f}s differs from fixed {TR_FIXED:.3f}s")

    vs_psc = psc_vs_runmean(ts_raw); np.savetxt(ts_dir / 'vs_psc.txt', vs_psc[np.newaxis, :], fmt='%.6f')
    vs_z   = zscore_run(ts_raw);     np.savetxt(ts_dir / 'vs_z.txt',   vs_z[np.newaxis, :],   fmt='%.6f')

    # EVs
    ev_dir = EVFILES / f"sub-{sub}" / f"ses-{ses}" / TASK / f"run-{run}"
    ev_reward  = ev_dir / '_anticipation_reward.txt'
    ev_neutral = ev_dir / '_anticipation_neutral.txt'
    if not ev_reward.exists() or not ev_neutral.exists():
        raise FileNotFoundError(f"Missing EVs: {ev_reward} or {ev_neutral}")

    times_R = load_ev_times(ev_reward,  ALIGN)
    times_N = load_ev_times(ev_neutral, ALIGN)

    # Peri windows for PSC
    tvec, wins_R_psc, nR_keep, nR_drop = extract_windows(vs_psc, times_R, TR_FIXED, TMIN, TMAX)
    _,    wins_N_psc, nN_keep, nN_drop = extract_windows(vs_psc, times_N, TR_FIXED, TMIN, TMAX)
    # Peri windows for Z
    _,    wins_R_z,   _, _           = extract_windows(vs_z,   times_R, TR_FIXED, TMIN, TMAX)
    _,    wins_N_z,   _, _           = extract_windows(vs_z,   times_N, TR_FIXED, TMIN, TMAX)
    # Trial retention log
    log_path = outdir / '_log.txt'
    with open(log_path, 'a') as lf:
        lf.write(f"{feat_base}: R kept={nR_keep}, dropped={nR_drop}; N kept={nN_keep}, dropped={nN_drop}\n")
    print(f"{feat_base}: R kept={nR_keep}, dropped={nR_drop}; N kept={nN_keep}, dropped={nN_drop}")

    # Stats
    psc_R_mean, psc_R_sem, nR, psc_R_sum, psc_R_sumsq = stats_from_windows(wins_R_psc)
    psc_N_mean, psc_N_sem, nN, psc_N_sum, psc_N_sumsq = stats_from_windows(wins_N_psc)
    z_R_mean,   z_R_sem,   _,  z_R_sum,   z_R_sumsq   = stats_from_windows(wins_R_z)
    z_N_mean,   z_N_sem,   _,  z_N_sum,   z_N_sumsq   = stats_from_windows(wins_N_z)

    # CSVs
    csv_psc = csv_dir / f"vs_t{int(TMIN)}_to_{int(TMAX)}_psc.csv"
    with open(csv_psc, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_s','Reward_mean','Reward_sem','Neutral_mean','Neutral_sem','nR','nN'])
        for i, t in enumerate(tvec):
            w.writerow([f"{t:.3f}", _fmt(psc_R_mean, i), _fmt(psc_R_sem, i), _fmt(psc_N_mean, i), _fmt(psc_N_sem, i), nR, nN])

    csv_z = csv_dir / f"vs_t{int(TMIN)}_to_{int(TMAX)}_z.csv"
    with open(csv_z, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_s','Reward_mean','Reward_sem','Neutral_mean','Neutral_sem','nR','nN'])
        for i, t in enumerate(tvec):
            w.writerow([f"{t:.3f}", _fmt(z_R_mean, i), _fmt(z_R_sem, i), _fmt(z_N_mean, i), _fmt(z_N_sem, i), nR, nN])

    # Plots
    fig_psc = (OUT_ROOT / feat_base / 'figs' / f"vs_t{int(TMIN)}_to_{int(TMAX)}_psc.png")
    plot_with_band(tvec, psc_R_mean, psc_R_sem, psc_N_mean, psc_N_sem,
                   f'Reward (n={nR})', f'Neutral (n={nN})',
                   f"VS — {feat_base} (PSC)", fig_psc)

    fig_z = (OUT_ROOT / feat_base / 'figs' / f"vs_t{int(TMIN)}_to_{int(TMAX)}_z.png")
    plot_with_band(tvec, z_R_mean, z_R_sem, z_N_mean, z_N_sem,
                   f'Reward (n={nR})', f'Neutral (n={nN})',
                   f"VS — {feat_base} (Z)", fig_z)

    return {
        'sub': sub, 'ses': ses, 'run': run, 'feat_base': feat_base, 't': tvec,
        'psc_R_mean': psc_R_mean, 'psc_R_sem': psc_R_sem, 'psc_N_mean': psc_N_mean, 'psc_N_sem': psc_N_sem,
        'z_R_mean': z_R_mean, 'z_R_sem': z_R_sem, 'z_N_mean': z_N_mean, 'z_N_sem': z_N_sem,
        'nR': nR, 'nN': nN,
        'psc_R_sum': psc_R_sum, 'psc_R_sumsq': psc_R_sumsq,
        'psc_N_sum': psc_N_sum, 'psc_N_sumsq': psc_N_sumsq,
        'z_R_sum': z_R_sum,   'z_R_sumsq':   z_R_sumsq,
        'z_N_sum': z_N_sum,   'z_N_sumsq':   z_N_sumsq,
    }


def _fmt(arr: np.ndarray, i: int) -> str:
    if arr.size == 0 or np.isnan(arr[i]):
        return ''
    return f"{float(arr[i]):.6f}"

# -------------------------
# Subject‑level aggregation
# -------------------------

def init_acc(n_time: int) -> Dict[str, np.ndarray]:
    return {
        'psc_R_sum': np.zeros(n_time), 'psc_R_sumsq': np.zeros(n_time), 'psc_R_cnt': np.zeros(n_time),
        'psc_N_sum': np.zeros(n_time), 'psc_N_sumsq': np.zeros(n_time), 'psc_N_cnt': np.zeros(n_time),
        'z_R_sum':   np.zeros(n_time), 'z_R_sumsq':   np.zeros(n_time), 'z_R_cnt':   np.zeros(n_time),
        'z_N_sum':   np.zeros(n_time), 'z_N_sumsq':   np.zeros(n_time), 'z_N_cnt':   np.zeros(n_time),
    }


def add_run_to_acc(acc: Dict[str, np.ndarray], rd: Dict[str, object]):
    for key in ['psc_R_sum','psc_R_sumsq','psc_N_sum','psc_N_sumsq','z_R_sum','z_R_sumsq','z_N_sum','z_N_sumsq']:
        acc[key] += rd[key]
    acc['psc_R_cnt'] += (rd['psc_R_sum'] == rd['psc_R_sum']) * rd['nR']  # add nR to non‑NaN positions
    acc['psc_N_cnt'] += (rd['psc_N_sum'] == rd['psc_N_sum']) * rd['nN']
    acc['z_R_cnt']   += (rd['z_R_sum']   == rd['z_R_sum'])   * rd['nR']
    acc['z_N_cnt']   += (rd['z_N_sum']   == rd['z_N_sum'])   * rd['nN']


def finalize_acc(acc: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    for prefix in ['psc_R','psc_N','z_R','z_N']:
        s  = acc[f'{prefix}_sum']
        ss = acc[f'{prefix}_sumsq']
        n  = acc[f'{prefix}_cnt']
        mean = np.divide(s, n, out=np.full_like(s, np.nan), where=n>0)
        var  = np.divide(ss - (s*s)/np.maximum(n,1), np.maximum(n-1, 1), out=np.full_like(s, np.nan), where=n>1)
        sem  = np.divide(np.sqrt(np.maximum(var, 0)), np.sqrt(n), out=np.full_like(s, np.nan), where=n>0)
        out[f'{prefix}_mean'] = mean
        out[f'{prefix}_sem']  = sem
    return out

# -------------------------
# FEAT list loading
# -------------------------

def read_feat_list(file: Path) -> List[Path]:
    feats: List[Path] = []
    if not file.exists():
        print(f"[WARN] FEAT list not found: {file}. Nothing to do.")
        return feats
    for line in file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        p = Path(line)
        if not p.name.endswith('.feat'):
            # allow trailing '+' before .feat in the file
            if p.name.endswith('.feat+'):
                p = Path(str(p)[:-1])
        if p.exists():
            feats.append(p)
        else:
            print(f"[WARN] Listed FEAT does not exist: {p}")
    return feats

# -------------------------
# Main
# -------------------------

def main():
    if not VS_MASK.exists():
        raise FileNotFoundError(f"VS mask not found: {VS_MASK}")

    ensure_dir(OUT_ROOT)

    feats = read_feat_list(FEAT_LIST)
    if len(feats) == 0:
        print("No FEATs to process (feat_paths.txt empty or missing). Exiting.")
        return

    # Group by subject
    by_sub: Dict[str, Dict[str, object]] = {}

    for feat_dir in feats:
        try:
            rd = process_feat(feat_dir, VS_MASK)
        except RuntimeError as e:
            print(f"[SKIP] {feat_dir.name}: {e}")
            continue
        except Exception as e:
            print(f"[ERROR] {feat_dir.name}: {e}")
            continue

        sub = rd['sub']
        tvec = rd['t']
        if sub not in by_sub:
            by_sub[sub] = {
                't': tvec,
                'acc': init_acc(len(tvec)),
                'nR_total': 0,
                'nN_total': 0,
            }
        add_run_to_acc(by_sub[sub]['acc'], rd)
        by_sub[sub]['nR_total'] += rd['nR']
        by_sub[sub]['nN_total'] += rd['nN']

    # Subject‑level outputs
    avg_dir = OUT_ROOT / '_subject_averages'; ensure_dir(avg_dir)

    for sub, bundle in sorted(by_sub.items()):
        tvec = bundle['t']; acc = bundle['acc']
        final = finalize_acc(acc)
        nR = int(bundle['nR_total']); nN = int(bundle['nN_total'])

        # CSVs
        with open(avg_dir / f"subject-{sub}_vs_timecourse_psc.csv", 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['time_s','Reward_mean','Reward_sem','Neutral_mean','Neutral_sem','nR','nN'])
            for i, t in enumerate(tvec):
                w.writerow([f"{t:.3f}",
                            _fmt(final['psc_R_mean'], i), _fmt(final['psc_R_sem'], i),
                            _fmt(final['psc_N_mean'], i), _fmt(final['psc_N_sem'], i), nR, nN])
        with open(avg_dir / f"subject-{sub}_vs_timecourse_z.csv", 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['time_s','Reward_mean','Reward_sem','Neutral_mean','Neutral_sem','nR','nN'])
            for i, t in enumerate(tvec):
                w.writerow([f"{t:.3f}",
                            _fmt(final['z_R_mean'], i), _fmt(final['z_R_sem'], i),
                            _fmt(final['z_N_mean'], i), _fmt(final['z_N_sem'], i), nR, nN])

        # Plots
        plot_with_band(tvec, final['psc_R_mean'], final['psc_R_sem'], final['psc_N_mean'], final['psc_N_sem'],
                       f'Reward (n={nR})', f'Neutral (n={nN})',
                       f"VS — subject {sub} (PSC; all sessions/runs)",
                       avg_dir / f"subject-{sub}_vs_timecourse_psc.png")

        plot_with_band(tvec, final['z_R_mean'], final['z_R_sem'], final['z_N_mean'], final['z_N_sem'],
                       f'Reward (n={nR})', f'Neutral (n={nN})',
                       f"VS — subject {sub} (Z; all sessions/runs)",
                       avg_dir / f"subject-{sub}_vs_timecourse_z.png")

        # Summary
        (avg_dir / f"subject-{sub}_summary.txt").write_text(
            f"Trials included — Reward: {nR}, Neutral: {nN} (cue‑aligned, TR={TR_FIXED:.3f}s, window {TMIN}..{TMAX}s)"
        )
        print(f"Subject {sub}: aggregated Reward n={nR}, Neutral n={nN}")


if __name__ == '__main__':
    main()
