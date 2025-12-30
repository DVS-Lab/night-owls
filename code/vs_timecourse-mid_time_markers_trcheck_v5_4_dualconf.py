#!/usr/bin/env python3
"""
VS time-course extraction with Pre-ITI Binning Diagnostic
---------------------------------------------------------
This version adds a specific diagnostic analysis:
Splitting 'Reward Anticipation' trials into 'Short Pre-ITI' vs 'Long Pre-ITI'
to test if the 'early peak' (at ~2-3s) is actually hemodynamic carryover
from the previous trial's feedback.

Key Features:
- Standard Reward vs Neutral Anticipation plots.
- NEW: Reward Anticipation split by Median Pre-ITI duration.
- Discrete 4th-TR summaries.

Assumptions:
- Input FEAT directories are listed in ./feat_paths-unsmoothed.txt
- EV files live at <rootdir>/derivatives/fsl/EVFiles/...
- VS mask at <rootdir>/masks/space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz
"""
from __future__ import annotations

import re
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

# --------------------------- Fixed project structure ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent
FSL_DERIV  = ROOT_DIR / "derivatives" / "fsl"
MASKS_DIR  = ROOT_DIR / "masks"
OUT_TC_DIR = ROOT_DIR / "derivatives" / "extractions" / "timecourses-mid-iti-diagnostic"
SUMMARY_DIR= ROOT_DIR / "derivatives" / "extractions"
EV_BASE    = FSL_DERIV / "EVFiles"

VS_MNI     = MASKS_DIR / "space-MNI152NLin6Asym_desc-NAcc_mask.nii.gz"
FEAT_LIST_PATH = SCRIPT_DIR / "feat_paths-unsmoothed.txt"

# --------------------------- Analysis parameters ------------------------------
TASK         = "mid"
TR_HARDCODE  = 1.615
K_AFTER      = 3      # 4th TR after onset
TMIN         = -4.0
TMAX         = 16.0
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
    
    # Standard Conditions
    ant_psc: Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    ant_z:   Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_psc:  Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    fb_z:    Dict[str, Tuple[np.ndarray, np.ndarray, int]]
    
    # NEW: ITI Binned Conditions (Reward Ant only)
    ant_psc_iti: Dict[str, Tuple[np.ndarray, np.ndarray, int]] # Keys: "Short", "Long"

    # Discrete Summaries
    tp_means_psc: Dict[str, float | None]
    tp_means_z:   Dict[str, float | None]
    tp_counts:    Dict[str, int]

    # Markers
    cue_dur_s: float = 0.75
    isi_median_s: float = float('nan')
    fb_delay_median_s: float = float('nan')

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

def load_ev_3col(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load onset, duration, amplitude. Returns empty arrays if missing."""
    if not path.exists():
        return np.array([]), np.array([]), np.array([])
    try:
        arr = np.loadtxt(path, ndmin=2)
    except Exception:
        return np.array([]), np.array([]), np.array([])
    if arr.size == 0:
        return np.array([]), np.array([]), np.array([])
    # Handle case where file exists but might be malformed/empty
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 2: 
        return arr[:, 0], np.zeros_like(arr[:,0]), np.ones_like(arr[:,0])
    if arr.shape[1] == 2:
        return arr[:, 0], arr[:, 1], np.ones_like(arr[:,0])
    return arr[:, 0], arr[:, 1], arr[:, 2]

def warn_tr(nifti_img: nib.Nifti1Image, feat: Path) -> float:
    hdr_tr = float(nifti_img.header.get_zooms()[3])
    if abs(hdr_tr - TR_HARDCODE) > 1e-3:
        print(f"[WARN] TR mismatch in {feat.name}: header={hdr_tr:.3f}s vs hardcoded={TR_HARDCODE:.3f}s. Using hardcoded.")
    return TR_HARDCODE

def load_vs_timeseries(feat: Path) -> Tuple[np.ndarray, float]:
    img_path = feat / "filtered_func_data.nii.gz"
    if not img_path.exists():
        raise FileNotFoundError(f"Missing filtered_func_data.nii.gz in {feat}")
    img = nib.load(str(img_path))
    mask_img = nib.load(str(VS_MNI))
    
    data = img.get_fdata()
    mask = mask_img.get_fdata() > 0
    if data.shape[:3] != mask.shape:
        raise ValueError(f"Mask dims {mask.shape} != data dims {data.shape[:3]}")
        
    ts = np.nanmean(data[mask, :], axis=0)
    tr = warn_tr(img, feat)
    return ts.astype(float), tr

def to_psc(ts: np.ndarray) -> np.ndarray:
    mean = np.nanmean(ts)
    if mean == 0: return np.zeros_like(ts)
    return ((ts / mean) - 1.0) * 100.0

def to_z(ts: np.ndarray) -> np.ndarray:
    mu = np.nanmean(ts)
    sd = np.nanstd(ts, ddof=1)
    if sd == 0 or not np.isfinite(sd): return np.zeros_like(ts)
    return (ts - mu) / sd

def build_time_axis(tmin: float, tmax: float, tr: float) -> np.ndarray:
    n = int(np.floor((tmax - tmin) / tr)) + 1
    return (tmin + np.arange(n) * tr).astype(float)

def sample_windows(ts: np.ndarray, onsets_s: np.ndarray, tr: float, tmin: float, tmax: float) -> np.ndarray:
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
    if onsets.size == 0: return np.empty((0,), dtype=int)
    idx = np.floor((onsets + K_AFTER * tr) / tr).astype(int)
    return idx[(idx >= 0) & (idx < n_vols)]

# --------------------------- ITI Calculation Logic ---------------------------

def calculate_pre_iti_split(target_onsets: np.ndarray, 
                           all_fb_onsets: np.ndarray, 
                           all_fb_durs: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    For each target onset (Anticipation), find the time elapsed since the 
    end of the IMMEDIATELY PRECEDING feedback event.
    Returns:
        onsets_short: onsets where pre-ITI < median
        onsets_long:  onsets where pre-ITI >= median
    """
    if target_onsets.size == 0:
        return np.array([]), np.array([])
    
    # Calculate feedback offsets (end times)
    fb_offsets = all_fb_onsets + all_fb_durs
    # Sort to be safe
    sort_idx = np.argsort(fb_offsets)
    fb_offsets = fb_offsets[sort_idx]
    
    itis = []
    valid_mask = [] # To track trials where we actually found a previous event
    
    for t_onset in target_onsets:
        # Find index of the last feedback ending before this onset
        # searchsorted returns index where t_onset would be inserted. 
        # index-1 is the event before.
        idx = np.searchsorted(fb_offsets, t_onset) - 1
        
        if idx >= 0:
            prev_offset = fb_offsets[idx]
            gap = t_onset - prev_offset
            # Sanity check: gap should be positive and reasonable (e.g. < 20s)
            # If gap is huge (e.g. >30s), it might be a run break or start, treat carefully
            itis.append(gap)
            valid_mask.append(True)
        else:
            # First trial of the run, or no previous feedback found
            itis.append(np.nan)
            valid_mask.append(False)
            
    itis = np.array(itis)
    valid_mask = np.array(valid_mask, dtype=bool)
    
    # Calculate split threshold based on valid ITIs
    valid_itis = itis[valid_mask]
    if valid_itis.size < 2:
        return np.array([]), np.array([])
        
    threshold = np.nanmedian(valid_itis)
    
    # Identify indices
    is_short = (itis < threshold) & valid_mask
    is_long  = (itis >= threshold) & valid_mask
    
    return target_onsets[is_short], target_onsets[is_long]

# --------------------------- Plotting ---------------------------

def plot_two_conditions(time_axis: np.ndarray,
                        condA: Tuple[np.ndarray, np.ndarray, int], labelA: str,
                        condB: Tuple[np.ndarray, np.ndarray, int], labelB: str,
                        title: str, ylabel: str, out_png: Path,
                        vlines: Optional[List[Tuple[float, str]]] = None):
    mA, sA, nA = condA
    mB, sB, nB = condB
    
    # If empty
    if len(mA) == 0 and len(mB) == 0:
        return

    fig, ax = plt.subplots(figsize=(8, 4.8))
    
    if len(mA) > 0:
        ax.plot(time_axis, mA, label=f"{labelA} (n={nA})", color='tab:blue')
        ax.fill_between(time_axis, mA - sA, mA + sA, alpha=0.2, color='tab:blue')
        
    if len(mB) > 0:
        ax.plot(time_axis, mB, label=f"{labelB} (n={nB})", color='tab:orange')
        ax.fill_between(time_axis, mB - sB, mB + sB, alpha=0.2, color='tab:orange')

    if vlines:
        for x, lab in vlines:
            ax.axvline(x, ls=":", lw=1, color='gray')
            if lab:
                y_top = ax.get_ylim()[1]
                ax.text(x, y_top, f" {lab}", rotation=90, va="bottom", ha="left", fontsize=8)

    ax.set_xlabel("Time from event (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper right")
    
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)

# --------------------------- Main Processing ---------------------------

def process_one_feat(feat: Path) -> RunResult | None:
    try:
        sub, ses, run = parse_sub_ses_run_from_feat(feat)
        echo = echo_from_feat(feat)
        ts_raw, tr = load_vs_timeseries(feat)
    except Exception as e:
        print(f"[WARN] Skipping {feat.name}: {e}")
        return None

    ts_psc = to_psc(ts_raw)
    ts_z   = to_z(ts_raw)
    T = ts_raw.size

    # Load EVs (Onsets AND Durations now)
    ev_dir = get_ev_dir(sub, ses, run)
    
    # Anticipation
    ant_R_ons, ant_R_dur, _ = load_ev_3col(ev_dir / "_anticipation_reward.txt")
    ant_N_ons, ant_N_dur, _ = load_ev_3col(ev_dir / "_anticipation_neutral.txt")
    
    # Feedback (Detailed for discrete summaries)
    fb_PR_ons, fb_PR_dur, _ = load_ev_3col(ev_dir / "_feedback_positive_reward.txt")
    fb_NR_ons, fb_NR_dur, _ = load_ev_3col(ev_dir / "_feedback_negative_reward.txt")
    fb_PN_ons, fb_PN_dur, _ = load_ev_3col(ev_dir / "_feedback_positive_neutral.txt")
    fb_NN_ons, fb_NN_dur, _ = load_ev_3col(ev_dir / "_feedback_negative_neutral.txt")

    # Aggregate Feedback for ITI calculation
    all_fb_onsets = np.concatenate([fb_PR_ons, fb_NR_ons, fb_PN_ons, fb_NN_ons])
    all_fb_durs   = np.concatenate([fb_PR_dur, fb_NR_dur, fb_PN_dur, fb_NN_dur])
    
    # --- BINNING LOGIC ---
    ant_R_short, ant_R_long = calculate_pre_iti_split(ant_R_ons, all_fb_onsets, all_fb_durs)

    # Markers for plotting
    isi_median = np.nanmedian(np.concatenate([ant_R_dur, ant_N_dur])) if (ant_R_dur.size or ant_N_dur.size) else np.nan
    
    # Calc median feedback delay
    all_ant_ons = np.sort(np.concatenate([ant_R_ons, ant_N_ons]))
    all_fb_ons_sorted = np.sort(all_fb_onsets)
    delays = []
    for t in all_ant_ons:
        idx = np.searchsorted(all_fb_ons_sorted, t)
        if idx < all_fb_ons_sorted.size:
            delays.append(all_fb_ons_sorted[idx] - t)
    fb_delay_median = np.nanmedian(delays) if delays else np.nan

    # Build Time Axis
    t_axis = build_time_axis(TMIN, TMAX, tr)

    # --- EXTRACTIONS ---
    # Standard
    ant_psc = {
        "Reward": mean_and_sem(sample_windows(ts_psc, ant_R_ons, tr, TMIN, TMAX)),
        "Neutral": mean_and_sem(sample_windows(ts_psc, ant_N_ons, tr, TMIN, TMAX))
    }
    ant_z = {
        "Reward": mean_and_sem(sample_windows(ts_z, ant_R_ons, tr, TMIN, TMAX)),
        "Neutral": mean_and_sem(sample_windows(ts_z, ant_N_ons, tr, TMIN, TMAX))
    }
    
    # Binned (Reward Ant Only)
    ant_psc_iti = {
        "Short": mean_and_sem(sample_windows(ts_psc, ant_R_short, tr, TMIN, TMAX)),
        "Long":  mean_and_sem(sample_windows(ts_psc, ant_R_long, tr, TMIN, TMAX))
    }
    
    # Feedback (Pooled Valence)
    fb_pos_ons = np.sort(np.concatenate([fb_PR_ons, fb_PN_ons]))
    fb_neg_ons = np.sort(np.concatenate([fb_NR_ons, fb_NN_ons]))
    
    fb_psc = {
        "Positive": mean_and_sem(sample_windows(ts_psc, fb_pos_ons, tr, TMIN, TMAX)),
        "Negative": mean_and_sem(sample_windows(ts_psc, fb_neg_ons, tr, TMIN, TMAX))
    }
    fb_z = {
        "Positive": mean_and_sem(sample_windows(ts_z, fb_pos_ons, tr, TMIN, TMAX)),
        "Negative": mean_and_sem(sample_windows(ts_z, fb_neg_ons, tr, TMIN, TMAX))
    }

    # Discrete Summaries
    def val_at_4th(ts, onsets):
        idx = fourth_tr_indices(onsets, tr, T)
        return ts[idx] if idx.size else np.array([])
    
    def m(arr): return float(np.nanmean(arr)) if arr.size else None

    tp_psc = {
        "ANT_REWARD": m(val_at_4th(ts_psc, ant_R_ons)),
        "ANT_NEUTRAL": m(val_at_4th(ts_psc, ant_N_ons)),
        "FB_POS_REWARD": m(val_at_4th(ts_psc, fb_PR_ons)),
        "FB_NEG_REWARD": m(val_at_4th(ts_psc, fb_NR_ons)),
        "FB_POS_NEUTRAL": m(val_at_4th(ts_psc, fb_PN_ons)),
        "FB_NEG_NEUTRAL": m(val_at_4th(ts_psc, fb_NN_ons)),
    }
    tp_z = {
        "ANT_REWARD": m(val_at_4th(ts_z, ant_R_ons)),
        "ANT_NEUTRAL": m(val_at_4th(ts_z, ant_N_ons)),
        "FB_POS_REWARD": m(val_at_4th(ts_z, fb_PR_ons)),
        "FB_NEG_REWARD": m(val_at_4th(ts_z, fb_NR_ons)),
        "FB_POS_NEUTRAL": m(val_at_4th(ts_z, fb_PN_ons)),
        "FB_NEG_NEUTRAL": m(val_at_4th(ts_z, fb_NN_ons)),
    }
    tp_counts = {
        "ANT_REWARD": int(ant_R_ons.size),
        "ANT_NEUTRAL": int(ant_N_ons.size),
        "FB_POS_REWARD": int(fb_PR_ons.size),
        "FB_NEG_REWARD": int(fb_NR_ons.size),
        "FB_POS_NEUTRAL": int(fb_PN_ons.size),
        "FB_NEG_NEUTRAL": int(fb_NN_ons.size),
    }

    return RunResult(
        feat_path=feat, sub=sub, ses=ses, run=run, echo=echo, time_axis=t_axis,
        ant_psc=ant_psc, ant_z=ant_z, fb_psc=fb_psc, fb_z=fb_z,
        ant_psc_iti=ant_psc_iti,
        tp_means_psc=tp_psc, tp_means_z=tp_z, tp_counts=tp_counts,
        isi_median_s=isi_median, fb_delay_median_s=fb_delay_median
    )

def save_run_plots(res: RunResult) -> None:
    run_out = OUT_TC_DIR / "runs" / f"sub-{res.sub}" / f"ses-{res.ses}" / f"run-{res.run}" / res.echo
    
    # Timing Markers
    vlines_ant = [(-res.cue_dur_s, 'Cue On'), (0.0, 'Ant On')]
    if not np.isnan(res.isi_median_s): vlines_ant.append((res.isi_median_s, 'Tgt On'))
    if not np.isnan(res.fb_delay_median_s): vlines_ant.append((res.fb_delay_median_s, 'FB On'))
    vlines_ant.append((6.0, '+6s'))

    # Standard Plots
    plot_two_conditions(res.time_axis, res.ant_psc["Reward"], "Reward", res.ant_psc["Neutral"], "Neutral",
                        f"ANT PSC - {res.sub}", "PSC", run_out / "anticipation_psc.png", vlines_ant)
    
    # NEW DIAGNOSTIC PLOT: ITI BINNING
    plot_two_conditions(res.time_axis, res.ant_psc_iti["Long"], "Long Pre-ITI", res.ant_psc_iti["Short"], "Short Pre-ITI",
                        f"Reward Ant by ITI - {res.sub}", "PSC", run_out / "anticipation_by_iti_psc.png", vlines_ant)

def aggregate_subject(results: List[RunResult], subject: str) -> None:
    pool: Dict[str, Dict] = {}
    time_axis = None
    
    # Aggregate data
    for res in results:
        if res.sub != subject: continue
        if time_axis is None: time_axis = res.time_axis
        
        echo = res.echo
        pool.setdefault(echo, {
            "ANT_PSC": {"Reward": [], "Neutral": []},
            "ANT_PSC_ITI": {"Short": [], "Long": []},
            "FB_PSC": {"Positive": [], "Negative": []}
        })
        
        pool[echo]["ANT_PSC"]["Reward"].append((res.ant_psc["Reward"][0], res.ant_psc["Reward"][2]))
        pool[echo]["ANT_PSC"]["Neutral"].append((res.ant_psc["Neutral"][0], res.ant_psc["Neutral"][2]))
        
        pool[echo]["ANT_PSC_ITI"]["Short"].append((res.ant_psc_iti["Short"][0], res.ant_psc_iti["Short"][2]))
        pool[echo]["ANT_PSC_ITI"]["Long"].append((res.ant_psc_iti["Long"][0], res.ant_psc_iti["Long"][2]))
        
        pool[echo]["FB_PSC"]["Positive"].append((res.fb_psc["Positive"][0], res.fb_psc["Positive"][2]))
        pool[echo]["FB_PSC"]["Negative"].append((res.fb_psc["Negative"][0], res.fb_psc["Negative"][2]))

    if time_axis is None: return
    
    subj_out = OUT_TC_DIR / "subjects" / f"sub-{subject}"
    
    def w_mean(items):
        if not items: return np.array([]), np.array([]), 0
        means = np.vstack([m for m, n in items if m.size])
        if means.size == 0: return np.array([]), np.array([]), 0
        ns = [n for m, n in items if m.size]
        avg = np.nanmean(means, axis=0)
        sem = np.nanstd(means, axis=0, ddof=1) / np.sqrt(len(means))
        return avg, sem, sum(ns)

    for echo, data in pool.items():
        # Standard Ant
        plot_two_conditions(time_axis, w_mean(data["ANT_PSC"]["Reward"]), "Reward", w_mean(data["ANT_PSC"]["Neutral"]), "Neutral",
                            f"Subject {subject} ANT PSC [{echo}]", "PSC", subj_out / echo / "anticipation_psc.png", [(0,'On'),(6,'+6')])
        
        # ITI Diagnostic
        plot_two_conditions(time_axis, w_mean(data["ANT_PSC_ITI"]["Long"]), "Long Pre-ITI", w_mean(data["ANT_PSC_ITI"]["Short"]), "Short Pre-ITI",
                            f"Subject {subject} ITI DIAGNOSTIC [{echo}]", "PSC", subj_out / echo / "anticipation_by_iti_psc.png", [(0,'On'),(6,'+6')])

def main():
    if not FEAT_LIST_PATH.exists():
        print("No FEAT list found.")
        return
    feat_paths = [Path(line.strip()) for line in FEAT_LIST_PATH.read_text().splitlines() if line.strip() and not line.startswith('#')]
    
    results = []
    for feat in feat_paths:
        res = process_one_feat(feat)
        if res:
            save_run_plots(res)
            results.append(res)
            
    subjects = sorted({r.sub for r in results})
    for sub in subjects:
        aggregate_subject(results, sub)
        
    print(f"Done. Check {OUT_TC_DIR} for 'anticipation_by_iti_psc.png' plots.")

if __name__ == "__main__":
    main()