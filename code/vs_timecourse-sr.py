#!/usr/bin/env python3
"""
VS time-course: Shared Reward (Decision-Locked Split by Outcome)
----------------------------------------------------------------
Updates:
1. Interpolation: 0.1s smooth curves.
2. Confounds: Regressed out.
3. Analysis Logic: 
   - Locks to DECISION (Guess) onset.
   - Ignores 'Computer' vs 'Stranger'.
   - Splits traces based on SUBSEQUENT Outcome (Reward, Neutral, Punish).
   * This visualizes the full trial trajectory: Decision -> Wait -> Outcome.
"""
from __future__ import annotations
import re
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# --------------------------- Config ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
MASKS_DIR  = ROOT_DIR / "masks"
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses-sr-final"
EV_BASE    = FSL_DERIV / "EVFiles"

VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
FEAT_LIST  = SCRIPT_DIR / "feat_paths-sr.txt"

TR_HARDCODE = 1.615
DT          = 0.1
TMIN, TMAX  = -4.0, 16.0

# --------------------------- Helpers ---------------------------
# (Same helpers as above, duplicated for standalone execution)
def load_ev_3col(path: Path):
    if not path.exists(): return np.array([]), np.array([]), np.array([])
    try: arr = np.loadtxt(path, ndmin=2)
    except: return np.array([]), np.array([]), np.array([])
    if arr.size == 0: return np.array([]), np.array([]), np.array([])
    if arr.ndim == 1: arr = arr.reshape(1, -1)
    if arr.shape[1] == 2: return arr[:,0], arr[:,1], np.ones_like(arr[:,0])
    return arr[:,0], arr[:,1], arr[:,2]

def regress_confounds(ts: np.ndarray, feat_path: Path) -> np.ndarray:
    conf_file = feat_path / "confoundevs.txt"
    if not conf_file.exists(): return ts
    try:
        X = np.loadtxt(conf_file)
        if X.ndim == 1: X = X.reshape(-1, 1)
        if X.shape[0] != ts.shape[0]: return ts
        X = np.column_stack([X, np.ones(X.shape[0])])
        beta = np.linalg.lstsq(X, ts, rcond=None)[0]
        return (ts - X @ beta) + np.mean(ts)
    except: return ts

def get_interpolated_windows(ts: np.ndarray, onsets: np.ndarray, tr: float) -> np.ndarray:
    if onsets.size == 0: return np.empty((0, 0))
    n_pts = int((TMAX - TMIN) / DT) + 1
    t_win = np.linspace(TMIN, TMAX, n_pts)
    t_raw = np.arange(ts.size) * tr
    windows = []
    for onset in onsets:
        val = np.interp(onset + t_win, t_raw, ts, left=np.nan, right=np.nan)
        windows.append(val)
    return np.vstack(windows)

def mean_sem(wins):
    if wins.size == 0: return np.array([]), np.array([]), 0
    mu = np.nanmean(wins, axis=0)
    sem = np.nanstd(wins, axis=0, ddof=1) / np.sqrt(wins.shape[0])
    return mu, sem, wins.shape[0]

# --------------------------- Main Logic ---------------------------
@dataclass
class RunResult:
    sub: str; ses: str; echo: str
    time_axis: np.ndarray
    dec_by_outcome: Dict # {OutcomeType: (Mean, SEM, N)}

def process_feat(feat: Path) -> Optional[RunResult]:
    try:
        m = re.search(r"sub-(\d+).*ses-(\d+).*run-(\d+)", feat.as_posix())
        sub, ses, run = m.groups()
        echo = "multi-echo" if "multi-echo" in feat.name else "single-echo"
    except: return None

    # Load TS
    img_p = feat / "filtered_func_data.nii.gz"
    if not img_p.exists(): return None
    img = nib.load(img_p)
    ts = np.nanmean(img.get_fdata()[nib.load(VS_MNI).get_fdata() > 0], axis=0)
    ts = regress_confounds(ts, feat)
    ts_psc = ((ts / np.mean(ts)) - 1) * 100
    tr = img.header.get_zooms()[3]

    # Load EVs
    ev_d = EV_BASE / f"sub-{sub}" / f"ses-{ses}" / "sharedreward" / f"run-{run}"
    
    # 1. Decisions (Combine Computer + Stranger)
    d_comp, _, _ = load_ev_3col(ev_d / "_guess_computer.txt")
    d_str, _, _  = load_ev_3col(ev_d / "_guess_face.txt")
    d_all = np.sort(np.concatenate([d_comp, d_str]))
    
    # 2. Outcomes (Load all types)
    out_types = {} # {Onset: 'Reward'/'Neutral'/'Punish'}
    for p in ev_d.glob("_outcome_*.txt"):
        ons, _, _ = load_ev_3col(p)
        # Determine valence from filename
        if "reward" in p.name: val = "Reward"
        elif "punish" in p.name: val = "Punish"
        else: val = "Neutral"
        
        for o in ons: out_types[o] = val
        
    all_outcomes = np.sort(np.array(list(out_types.keys())))
    
    # 3. Match Decision -> Subsequent Outcome
    bins = {"Reward": [], "Neutral": [], "Punish": []}
    
    for dec_on in d_all:
        # Find closest outcome occurring AFTER decision
        # searchsorted returns insertion pt. We want index where outcome > dec_on
        idx = np.searchsorted(all_outcomes, dec_on)
        
        if idx < len(all_outcomes):
            out_on = all_outcomes[idx]
            gap = out_on - dec_on
            # Sanity check: Gap should be small (e.g. < 6s)
            if gap < 8.0:
                val = out_types[out_on]
                bins[val].append(dec_on)
                
    # Extract
    t_ax = np.linspace(TMIN, TMAX, int((TMAX-TMIN)/DT)+1)
    results = {}
    for k in bins:
        results[k] = mean_sem(get_interpolated_windows(ts_psc, np.array(bins[k]), tr))
        
    return RunResult(sub, ses, echo, t_ax, results)

def aggregate_and_plot(results: List[RunResult]):
    subs = set(r.sub for r in results)
    for sub in subs:
        sub_res = [r for r in results if r.sub == sub]
        for echo in ["multi-echo", "single-echo"]:
            echo_res = [r for r in sub_res if r.echo == echo]
            if not echo_res: continue
            
            t = echo_res[0].time_axis
            out = OUT_TC_DIR / f"sub-{sub}" / echo
            out.mkdir(parents=True, exist_ok=True)
            
            plt.figure(figsize=(10,6))
            cols = {"Reward": "green", "Neutral": "gray", "Punish": "red"}
            
            for k, col in cols.items():
                # Pool means
                all_means = [r.dec_by_outcome[k][0] for r in echo_res if r.dec_by_outcome[k][2] > 0]
                if not all_means: continue
                stack = np.vstack(all_means)
                mu = np.mean(stack, axis=0)
                sem = np.std(stack, axis=0, ddof=1) / np.sqrt(len(all_means))
                
                plt.plot(t, mu, label=k, color=col)
                plt.fill_between(t, mu-sem, mu+sem, alpha=0.1, color=col)
                
            plt.title(f"Decision-Locked split by Future Outcome - {sub} {echo}")
            plt.xlabel("Time from Decision Onset (s)")
            plt.ylabel("PSC")
            plt.axvline(0, color='k', ls='-', label="Decision")
            plt.axvline(3.5, color='k', ls=':', label="Approx Outcome")
            plt.legend()
            plt.savefig(out / "decision_by_outcome.png"); plt.close()

def main():
    if not FEAT_LIST.exists(): return
    paths = [Path(x.strip()) for x in FEAT_LIST.read_text().splitlines() 
             if x.strip() and not x.startswith("#")]
    results = [r for p in paths if (r := process_feat(p))]
    aggregate_and_plot(results)
    print("Done.")

if __name__ == "__main__":
    main()