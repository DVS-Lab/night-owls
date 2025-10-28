#!/usr/bin/env python3
"""
vs_timecourse.py — VS‑only, single‑purpose, HARD‑CODED batch extractor

Run from <rootdir>/code. All paths are resolved relative to this script.

What it does
- Assumes functional data and mask are already in **space‑MNI152NLin6Asym** on the same grid (no transforms).
- Hard‑codes EV locations for Reward/Neutral anticipation under EVFiles.
- Iterates **all subjects** found in EVFiles, and for **each subject** iterates **sessions 01–12** (skipping a fixed map).
- Within each session, discovers available **runs** from EVFiles and matches the FEAT dir in derivatives/fsl by name tokens.
- Extracts VS mean time series → PSC vs run mean → peristimulus averages (cue‑onset locked) for Reward & Neutral.
- Writes per‑run outputs and a **per‑subject average** curve (weighted by trial counts across all included runs/sessions).

How to run
    python vs_timecourse.py

Outputs
  derivatives/extractions/vs_timecourse/
   ├── <feat_basename>/masks/vs_func_assumed_mni.nii.gz
   ├── <feat_basename>/ts/vs_ts.txt
   ├── <feat_basename>/csv/vs_t-2_to_16.csv
   ├── <feat_basename>/csv/vs_peaks.csv
   ├── <feat_basename>/figs/vs_t-2_to_16.png
   └── _subject_averages/subject-<sub>_vs_timecourse.csv/.png and subject-<sub>_summary.txt

If a grid mismatch is detected between the functional data and the VS mask, the run is skipped with an error message.
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

VS_MASK = MASKS_DIR / "space-MNI152NLin6Asym_desc-VS-Imanova_mask.nii.gz"

# -------------------------
# HARD‑CODED SETTINGS
# -------------------------
TASKS       = ["mid"]
SPACE_TOKS  = ["space-mni", "space-MNI152NLin6Asym"]
ECHO_TOKS   = ["single-echo", "multi-echo"]
CONFOUND_TOKS = ["cnfds-fmriprep"]

ALIGN       = "onset"   # lock to cue onset
TMIN        = -2.0
TMAX        = 16.0
PEAK_LAGS   = [6.0, 4.0]

SESSIONS    = [f"{i:02d}" for i in range(1, 13)]  # 01..12
SKIP_MAP: Dict[str, set] = {
    "101": {"04", "05", "12"},
    "103": {"12"},
}

# -------------------------
# Small helpers
# -------------------------

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def load_nifti(path: Path):
    img = nib.load(str(path))
    return img.get_fdata(), img.affine, img.header


def get_tr(path: Path) -> float:
    hdr = nib.load(str(path)).header
    tr = float(hdr.get_zooms()[3])
    if tr <= 0:
        raise ValueError("Invalid TR in functional data header")
    return tr


def psc_vs_runmean(ts: np.ndarray) -> np.ndarray:
    mu = float(np.mean(ts))
    if mu == 0:
        raise ValueError("Run mean is zero; PSC undefined")
    return 100.0 * (ts - mu) / mu


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


def peri_average(psc: np.ndarray, onsets: np.ndarray, tr: float, tmin: float, tmax: float) -> Tuple[np.ndarray, np.ndarray, int]:
    n_win = int(np.floor((tmax - tmin) / tr)) + 1
    times = np.arange(n_win) * tr + tmin
    T = psc.shape[0]
    segs = []
    for onset in onsets:
        idx0  = int(round(onset / tr))
        start = idx0 + int(round(tmin / tr))
        end   = start + n_win
        if start < 0 or end > T:
            continue
        segs.append(psc[start:end])
    if len(segs) == 0:
        return times, np.full(n_win, np.nan), 0
    return times, np.vstack(segs).mean(axis=0), len(segs)

# -------------------------
# Discovery helpers
# -------------------------

def subjects_from_ev(ev_root: Path) -> List[str]:
    return sorted([p.name.split('-')[-1] for p in ev_root.glob('sub-*') if p.is_dir()])


def runs_from_ev(ev_task_dir: Path) -> List[str]:
    runs: List[str] = []
    for p in sorted(ev_task_dir.glob('run-*')):
        if not p.is_dir():
            continue
        if (p / '_anticipation_reward.txt').exists() and (p / '_anticipation_neutral.txt').exists():
            runs.append(p.name.split('-')[-1])
    return runs


def find_feat(sub: str, ses: str, task: str, run: str) -> Optional[Path]:
    base = DERIV_FSL / f"sub-{sub}" / f"ses-{ses}"
    if not base.exists():
        return None
    for feat in base.rglob('*.feat'):
        name = feat.name
        if f"task-{task}" not in name or f"run-{run}" not in name:
            continue
        if not any(sp in name for sp in SPACE_TOKS):
            continue
        if not any(e in name for e in ECHO_TOKS):
            continue
        if not any(c in name for c in CONFOUND_TOKS):
            continue
        return feat
    return None

# -------------------------
# Plotting
# -------------------------

def plot_twocond(times: np.ndarray, r: np.ndarray, n: np.ndarray, title: str, vlines: List[float], out_png: Path, nR: int, nN: int):
    plt.figure(figsize=(6.5, 4.0))
    if not np.all(np.isnan(r)):
        plt.plot(times, r, label=f'Reward (n={nR})')
    if not np.all(np.isnan(n)):
        plt.plot(times, n, label=f'Neutral (n={nN})')
    plt.axvline(0, linestyle='--', linewidth=1)
    for v in vlines:
        plt.axvline(v, linestyle=':', linewidth=1)
    plt.xlabel('Time from cue onset (s)')
    plt.ylabel('% signal change (PSC)')
    plt.title(title)
    plt.legend(frameon=False)
    plt.tight_layout(); plt.savefig(out_png, dpi=150); plt.close()

# -------------------------
# Core per‑run computation
# -------------------------

def fmt_f(x: float) -> str:
    return '' if np.isnan(x) else f"{float(x):.6f}"


def run_one(feat_dir: Path, sub: str, ses: str, task: str, run: str, vs_mask: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    func = feat_dir / 'filtered_func_data.nii.gz'
    if not func.exists():
        raise FileNotFoundError(f"Missing {func}")

    # per‑run outdir
    feat_base = feat_dir.name.replace('.feat','')
    outdir = OUT_ROOT / feat_base
    masks_dir = outdir / 'masks'
    ts_dir    = outdir / 'ts'
    csv_dir   = outdir / 'csv'
    fig_dir   = outdir / 'figs'
    for d in (masks_dir, ts_dir, csv_dir, fig_dir):
        ensure_dir(d)

    # Load and grid‑check
    func_data, func_aff, func_hdr = load_nifti(func)
    vs_data, vs_aff, vs_hdr = load_nifti(vs_mask)
    if func_data.shape[:3] != vs_data.shape[:3] or not np.allclose(func_aff, vs_aff, atol=1e-3):
        raise RuntimeError(f"Grid mismatch for {feat_dir}: functional and VS mask differ; pre‑align before running.")

    # Save the (assumed) func‑grid VS mask for bookkeeping
    nib.Nifti1Image((vs_data>0).astype(np.uint8), vs_aff, vs_hdr).to_filename(str(masks_dir / 'vs_func_assumed_mni.nii.gz'))

    # Extract TS and PSC
    mask = vs_data > 0
    T = func_data.shape[3]
    ts_raw = func_data[mask].reshape((-1, T)).mean(axis=0)
    np.savetxt(ts_dir / 'vs_ts.txt', ts_raw[np.newaxis, :], fmt='%.6f')
    tr = get_tr(func)
    psc = psc_vs_runmean(ts_raw)

    # EVs (hard‑coded layout)
    ev_dir = EVFILES / f"sub-{sub}" / f"ses-{ses}" / task / f"run-{run}"
    ev_reward  = ev_dir / '_anticipation_reward.txt'
    ev_neutral = ev_dir / '_anticipation_neutral.txt'
    if not ev_reward.exists() or not ev_neutral.exists():
        raise FileNotFoundError(f"Missing EVs for sub-{sub} ses-{ses} {task} run-{run}")

    times_reward  = load_ev_times(ev_reward,  ALIGN)
    times_neutral = load_ev_times(ev_neutral, ALIGN)

    # Peri‑stimulus averages
    n_win = int(np.floor((TMAX - TMIN) / tr)) + 1
    tvec = np.arange(n_win) * tr + TMIN
    # use helper for correctness
    tvec, curve_reward,  nR = peri_average(psc, times_reward,  tr, TMIN, TMAX)
    _,    curve_neutral, nN = peri_average(psc, times_neutral, tr, TMIN, TMAX)

    # CSV time course
    with open(csv_dir / f"vs_t{int(TMIN)}_to_{int(TMAX)}.csv", 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_s','Reward','Neutral'])
        for i, t in enumerate(tvec):
            w.writerow([f"{t:.3f}", fmt_f(curve_reward[i]), fmt_f(curve_neutral[i])])

    # Peaks
    with open(csv_dir / 'vs_peaks.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['roi','lag_s','Reward','Neutral'])
        for lag in PEAK_LAGS:
            idx = int(np.argmin(np.abs(tvec - lag)))
            w.writerow(['VS', lag, fmt_f(curve_reward[idx]), fmt_f(curve_neutral[idx])])

    # Plot
    out_png = fig_dir / f"vs_t{int(TMIN)}_to_{int(TMAX)}.png"
    plot_twocond(tvec, curve_reward, curve_neutral, f"VS — {feat_base}", PEAK_LAGS, out_png, nR, nN)

    return tvec, curve_reward, curve_neutral, nR, nN

# -------------------------
# Subject‑level aggregation
# -------------------------

def init_accumulator(n_time: int) -> Dict[str, np.ndarray]:
    return { 'sum_R': np.zeros(n_time), 'sum_N': np.zeros(n_time), 'cnt_R': np.zeros(n_time), 'cnt_N': np.zeros(n_time) }


def add_to_acc(acc: Dict[str, np.ndarray], r: np.ndarray, n: np.ndarray, nR: int, nN: int):
    if nR > 0 and not np.any(np.isnan(r)):
        acc['sum_R'] += r * nR
        acc['cnt_R'] += (r == r) * nR
    if nN > 0 and not np.any(np.isnan(n)):
        acc['sum_N'] += n * n
        acc['cnt_N'] += (n == n) * nN


def finalize_acc(acc: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    R = np.divide(acc['sum_R'], acc['cnt_R'], out=np.full_like(acc['sum_R'], np.nan), where=acc['cnt_R']>0)
    N = np.divide(acc['sum_N'], acc['cnt_N'], out=np.full_like(acc['sum_N'], np.nan), where=acc['cnt_N']>0)
    return R, N

# -------------------------
# Main (no CLI; everything hard‑coded)
# -------------------------

def main():
    # Verify VS mask
    if not VS_MASK.exists():
        raise FileNotFoundError(f"VS mask not found: {VS_MASK}")

    subjects = subjects_from_ev(EVFILES)

    for sub in subjects:
        acc = None
        tvec_ref = None
        total_R = total_N = 0

        for ses in SESSIONS:
            if ses in SKIP_MAP.get(sub, set()):
                print(f"Skipping sub-{sub} ses-{ses} (per SKIP map)")
                continue
            for task in TASKS:
                ev_task_dir = EVFILES / f"sub-{sub}" / f"ses-{ses}" / task
                if not ev_task_dir.exists():
                    continue
                runs = runs_from_ev(ev_task_dir)
                for run in runs:
                    feat = find_feat(sub, ses, task, run)
                    if feat is None:
                        print(f"No matching FEAT for sub-{sub} ses-{ses} {task} run-{run}; skipping")
                        continue
                    try:
                        tvec, R, N, nR, nN = run_one(feat, sub, ses, task, run, VS_MASK)
                    except Exception as e:
                        print(f"ERROR in {feat}: {e}")
                        continue
                    if acc is None:
                        acc = init_accumulator(len(tvec))
                        tvec_ref = tvec
                    add_to_acc(acc, R, N, nR, nN)
                    total_R += nR; total_N += nN

        # finalize subject average
        avg_dir = OUT_ROOT / '_subject_averages'
        ensure_dir(avg_dir)
        if acc is not None:
            R_subj, N_subj = finalize_acc(acc)
            # CSV
            csv_path = avg_dir / f"subject-{sub}_vs_timecourse.csv"
            with open(csv_path, 'w', newline='') as f:
                w = csv.writer(f); w.writerow(['time_s','Reward','Neutral'])
                for i, t in enumerate(tvec_ref):
                    w.writerow([f"{t:.3f}", fmt_f(R_subj[i]), fmt_f(N_subj[i])])
            # Plot
            png_path = avg_dir / f"subject-{sub}_vs_timecourse.png"
            plot_twocond(tvec_ref, R_subj, N_subj, f"VS — subject {sub} (all sessions/runs)", PEAK_LAGS, png_path, total_R, total_N)
            # Summary
            (avg_dir / f"subject-{sub}_summary.txt").write_text(
                f"Trials included — Reward: {total_R}, Neutral: {total_N}
"
            )
            print(f"Subject {sub}: aggregated Reward n={total_R}, Neutral n={total_N}")
        else:
            print(f"Subject {sub}: nothing processed (all sessions missing or skipped)")


if __name__ == '__main__':
    main()
