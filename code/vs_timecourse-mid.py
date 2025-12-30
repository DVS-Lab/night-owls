#!/usr/bin/env python3
"""
VS time-course: MID Task with "Previous-Trial" Binning
------------------------------------------------------
Updates:
1. Interpolation: Samples at 0.1s resolution (smoothing).
2. Confounds: Regresses out 'confoundevs.txt' if present.
3. Diagnostic Binning: Splits 'Reward Anticipation' into 4 categories:
   - Prev Outcome Positive & Short ITI
   - Prev Outcome Positive & Long ITI
   - Prev Outcome Negative & Short ITI
   - Prev Outcome Negative & Long ITI
   * This is critical to diagnosing if the 'early peak' is signal bleed
     from a previous positive reward.

Assumptions:
- FEAT directories listed in 'feat_paths-unsmoothed.txt'
- EVs in standard path structure.
"""
from __future__ import annotations

import re
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# --------------------------- Config ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
MASKS_DIR  = ROOT_DIR / "masks"
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-final"
EV_BASE    = FSL_DERIV / "EVFiles"

VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
FEAT_LIST  = SCRIPT_DIR / "feat_paths-unsmoothed.txt"

# Params
TR_HARDCODE  = 1.615
DT           = 0.1    # Interpolation grid (seconds)
TMIN, TMAX   = -4.0, 16.0

# --------------------------- Helpers ---------------------------
def load_ev_3col(path: Path):
    if not path.exists(): return np.array([]), np.array([]), np.array([])
    try: arr = np.loadtxt(path, ndmin=2)
    except: return np.array([]), np.array([]), np.array([])
    if arr.size == 0: return np.array([]), np.array([]), np.array([])
    if arr.ndim == 1: arr = arr.reshape(1, -1)
    # Ensure 3 cols
    if arr.shape[1] == 2: return arr[:,0], arr[:,1], np.ones_like(arr[:,0])
    return arr[:,0], arr[:,1], arr[:,2]

def regress_confounds(ts: np.ndarray, feat_path: Path) -> np.ndarray:
    """Load confoundevs.txt and regress out of timeseries."""
    conf_file = feat_path / "confoundevs.txt"
    if not conf_file.exists():
        return ts
    try:
        X = np.loadtxt(conf_file)
        if X.ndim == 1: X = X.reshape(-1, 1)
        
        # Check dimensions
        if X.shape[0] != ts.shape[0]:
            print(f"[WARN] Confound dim mismatch: {X.shape} vs {ts.shape} in {feat_path.name}")
            return ts
            
        # Add intercept
        X = np.column_stack([X, np.ones(X.shape[0])])
        
        # Regress
        beta = np.linalg.lstsq(X, ts, rcond=None)[0]
        y_pred = X @ beta
        resid = ts - y_pred
        return resid + np.mean(ts) # Add mean back
    except Exception as e:
        print(f"[WARN] Confound regression failed for {feat_path.name}: {e}")
        return ts

def get_interpolated_windows(ts: np.ndarray, onsets: np.ndarray, tr: float) -> np.ndarray:
    if onsets.size == 0: return np.empty((0, 0))
    
    # Create high-res time axis
    n_pts = int((TMAX - TMIN) / DT) + 1
    t_win = np.linspace(TMIN, TMAX, n_pts)
    
    # Original time axis
    t_raw = np.arange(ts.size) * tr
    
    windows = []
    for onset in onsets:
        # Interpolate
        t_sample = onset + t_win
        val = np.interp(t_sample, t_raw, ts, left=np.nan, right=np.nan)
        windows.append(val)
    return np.vstack(windows)

def mean_sem(wins: np.ndarray):
    if wins.size == 0: return np.array([]), np.array([]), 0
    # Nan handling
    n = np.sum(~np.isnan(wins), axis=0)
    mu = np.nanmean(wins, axis=0)
    sem = np.nanstd(wins, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))
    return mu, sem, wins.shape[0]

# --------------------------- Main Logic ---------------------------
@dataclass
class RunResult:
    sub: str; ses: str; echo: str
    time_axis: np.ndarray
    # Data: {Condition: (Mean, SEM, N)}
    ant_std: Dict
    ant_bin: Dict # The binning diagnostic
    fb_std: Dict

def process_feat(feat: Path) -> Optional[RunResult]:
    # Parse path
    try:
        m = re.search(r"sub-(\d+).*ses-(\d+).*run-(\d+)", feat.as_posix())
        sub, ses, run = m.groups()
        echo = "multi-echo" if "multi-echo" in feat.name else "single-echo"
    except: return None

    # Load Data
    img_p = feat / "filtered_func_data.nii.gz"
    if not img_p.exists(): return None
    
    img = nib.load(img_p)
    data = img.get_fdata()
    mask = nib.load(VS_MNI).get_fdata() > 0
    ts = np.nanmean(data[mask], axis=0)
    
    # Confound Regression
    ts = regress_confounds(ts, feat)
    
    # Normalize
    ts_psc = ((ts / np.mean(ts)) - 1) * 100
    tr = img.header.get_zooms()[3]
    if abs(tr - TR_HARDCODE) > 0.01: tr = TR_HARDCODE # Force match

    # Load EVs
    ev_d = EV_BASE / f"sub-{sub}" / f"ses-{ses}" / "mid" / f"run-{run}"
    
    ant_r_on, _, _ = load_ev_3col(ev_d / "_anticipation_reward.txt")
    ant_n_on, _, _ = load_ev_3col(ev_d / "_anticipation_neutral.txt")
    
    # Feedback
    fb_pos_rew, fb_pos_rew_dur, _ = load_ev_3col(ev_d / "_feedback_positive_reward.txt")
    fb_neg_rew, fb_neg_rew_dur, _ = load_ev_3col(ev_d / "_feedback_negative_reward.txt")
    fb_pos_neu, fb_pos_neu_dur, _ = load_ev_3col(ev_d / "_feedback_positive_neutral.txt")
    fb_neg_neu, fb_neg_neu_dur, _ = load_ev_3col(ev_d / "_feedback_negative_neutral.txt")

    # --- BINNING LOGIC ---
    # 1. Define "Positive" vs "Negative" Feedback from *Previous* trial
    # Pos = Win Money (PosRew) OR Win Neutral (PosNeu - usually 'Hit' but no money)
    # Neg = Miss Money (NegRew) OR Miss Neutral (NegNeu)
    fb_all_pos = np.concatenate([fb_pos_rew, fb_pos_neu])
    fb_all_pos_dur = np.concatenate([fb_pos_rew_dur, fb_pos_neu_dur])
    
    fb_all_neg = np.concatenate([fb_neg_rew, fb_neg_neu])
    fb_all_neg_dur = np.concatenate([fb_neg_rew_dur, fb_neg_neu_dur])
    
    # Combine into (OffsetTime, ValenceCode) list. Code: 1=Pos, 0=Neg
    offsets = [] # (time, valence)
    for t, d in zip(fb_all_pos, fb_all_pos_dur): offsets.append((t+d, 1))
    for t, d in zip(fb_all_neg, fb_all_neg_dur): offsets.append((t+d, 0))
    
    # Sort by time
    offsets.sort(key=lambda x: x[0])
    fb_end_times = np.array([x[0] for x in offsets])
    fb_valences  = np.array([x[1] for x in offsets])
    
    # 2. For each Reward Anticipation, find prev trial
    bins = {"PrevPos_Short": [], "PrevPos_Long": [], "PrevNeg_Short": [], "PrevNeg_Long": []}
    
    if ant_r_on.size > 0 and fb_end_times.size > 0:
        # Calculate ITIs
        itis = []
        valences = []
        valid_mask = []
        
        for onset in ant_r_on:
            # Find closest feedback that ended BEFORE this anticipation started
            # searchsorted returns insertion index. idx-1 is the element before.
            idx = np.searchsorted(fb_end_times, onset) - 1
            if idx >= 0:
                gap = onset - fb_end_times[idx]
                if gap < 30: # Sanity check (ignore run breaks)
                    itis.append(gap)
                    valences.append(fb_valences[idx])
                    valid_mask.append(True)
                    continue
            valid_mask.append(False)
            
        itis = np.array(itis)
        valences = np.array(valences)
        valid_mask = np.array(valid_mask)
        
        # Median Split
        if len(itis) > 2:
            med = np.median(itis)
            # Assign bins
            # Using original 'ant_r_on' indices
            curr_valid_idx = 0
            for i, onset in enumerate(ant_r_on):
                if valid_mask[i]:
                    iti = itis[curr_valid_idx]
                    val = valences[curr_valid_idx]
                    curr_valid_idx += 1
                    
                    if val == 1: # Pos
                        k = "PrevPos_Short" if iti < med else "PrevPos_Long"
                    else: # Neg
                        k = "PrevNeg_Short" if iti < med else "PrevNeg_Long"
                    bins[k].append(onset)

    # Convert lists to arrays
    for k in bins: bins[k] = np.array(bins[k])

    # --- EXTRACTION ---
    time_ax = np.linspace(TMIN, TMAX, int((TMAX-TMIN)/DT)+1)
    
    ant_std = {
        "Reward": mean_sem(get_interpolated_windows(ts_psc, ant_r_on, tr)),
        "Neutral": mean_sem(get_interpolated_windows(ts_psc, ant_n_on, tr))
    }
    
    ant_bin = {k: mean_sem(get_interpolated_windows(ts_psc, v, tr)) for k,v in bins.items()}
    
    # Feedback (Pooled)
    fb_pos_on = np.sort(np.concatenate([fb_pos_rew, fb_pos_neu]))
    fb_neg_on = np.sort(np.concatenate([fb_neg_rew, fb_neg_neu]))
    
    fb_std = {
        "Positive": mean_sem(get_interpolated_windows(ts_psc, fb_pos_on, tr)),
        "Negative": mean_sem(get_interpolated_windows(ts_psc, fb_neg_on, tr))
    }

    return RunResult(sub, ses, echo, time_ax, ant_std, ant_bin, fb_std)

def aggregate_and_plot(results: List[RunResult]):
    subs = set(r.sub for r in results)
    for sub in subs:
        sub_res = [r for r in results if r.sub == sub]
        for echo in ["multi-echo", "single-echo"]:
            echo_res = [r for r in sub_res if r.echo == echo]
            if not echo_res: continue
            
            # Helper to pool
            def pool(cat, key):
                # cat: 'ant_std', 'ant_bin', etc
                # key: 'Reward', 'PrevPos_Short', etc
                all_means = []
                for r in echo_res:
                    m, s, n = getattr(r, cat)[key]
                    if n > 0: all_means.append(m)
                if not all_means: return np.array([]), np.array([]), 0
                # Weighted average (simplified to equal weight for visual check)
                stack = np.vstack(all_means)
                mu = np.mean(stack, axis=0)
                sem = np.std(stack, axis=0, ddof=1) / np.sqrt(len(all_means))
                return mu, sem, len(all_means)

            t = echo_res[0].time_axis
            out = OUT_TC_DIR / f"sub-{sub}" / echo
            out.mkdir(parents=True, exist_ok=True)

            # 1. Standard Anticipation
            plt.figure(figsize=(8,5))
            for k, col in zip(["Reward", "Neutral"], ["blue", "orange"]):
                m, s, n = pool("ant_std", k)
                if n: 
                    plt.plot(t, m, label=k, color=col)
                    plt.fill_between(t, m-s, m+s, alpha=0.2, color=col)
            plt.title(f"Anticipation (Standard) - {sub} {echo}")
            plt.axvline(0, color='k', ls=':')
            plt.legend()
            plt.savefig(out / "ant_standard.png"); plt.close()

            # 2. Diagnostic Binning (The critical one)
            plt.figure(figsize=(10,6))
            colors = {"PrevPos_Short": "red", "PrevPos_Long": "blue",
                      "PrevNeg_Short": "orange", "PrevNeg_Long": "green"}
            styles = {"Short": "-", "Long": "--"}
            
            for k, col in colors.items():
                m, s, n = pool("ant_bin", k)
                if n:
                    dur_label = k.split('_')[1] # Short/Long
                    st = styles[dur_label]
                    plt.plot(t, m, label=k, color=col, ls=st)
                    plt.fill_between(t, m-s, m+s, alpha=0.1, color=col)
            
            plt.title(f"Anticipation Binned by PREV Outcome & ITI - {sub}")
            plt.axvline(0, color='k', ls=':')
            plt.legend()
            plt.savefig(out / "ant_diagnostic_binning.png"); plt.close()

            # 3. Feedback
            plt.figure(figsize=(8,5))
            for k, col in zip(["Positive", "Negative"], ["green", "red"]):
                m, s, n = pool("fb_std", k)
                if n:
                    plt.plot(t, m, label=k, color=col)
                    plt.fill_between(t, m-s, m+s, alpha=0.2, color=col)
            plt.title(f"Feedback - {sub} {echo}")
            plt.axvline(0, color='k', ls=':')
            plt.legend()
            plt.savefig(out / "feedback.png"); plt.close()

def main():
    if not FEAT_LIST.exists(): return
    paths = [Path(x.strip()) for x in FEAT_LIST.read_text().splitlines() 
             if x.strip() and not x.startswith("#")]
    
    results = []
    for p in paths:
        r = process_feat(p)
        if r: results.append(r)
        
    aggregate_and_plot(results)
    print("Done.")

if __name__ == "__main__":
    main()