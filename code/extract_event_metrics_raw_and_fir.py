#!/usr/bin/env python3
"""
Extract summary metrics from:
  (1) Raw Wu-normalized ROI time courses (run_curves_long_wu.csv[.gz])
  (2) FIR beta time courses (fir_copes_roi_means_long.csv)

Overwrites (in --outdir):
  - metrics_raw_long.csv
  - metrics_raw_contrasts.csv
  - metrics_fir_long.csv
  - metrics_fir_contrasts.csv
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ----------------------------
# Core time-series utilities
# ----------------------------

def _ensure_boundary_points(
    t: np.ndarray, y: np.ndarray, t0: float, t1: float
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Return (t_seg, y_seg) clipped to [t0, t1] with boundary points inserted via linear interpolation."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if t.size == 0:
        return None, None

    order = np.argsort(t)
    t = t[order]
    y = y[order]

    tmin = float(t[0])
    tmax = float(t[-1])
    if t1 < tmin or t0 > tmax:
        return None, None

    t0c = max(t0, tmin)
    t1c = min(t1, tmax)
    if t1c <= t0c:
        return None, None

    mask = (t >= t0c) & (t <= t1c)
    tt = t[mask]
    yy = y[mask]

    # Insert boundary points if needed
    if tt.size == 0 or abs(tt[0] - t0c) > 1e-9:
        y0 = float(np.interp(t0c, t, y))
        tt = np.insert(tt, 0, t0c)
        yy = np.insert(yy, 0, y0)

    if tt.size == 0 or abs(tt[-1] - t1c) > 1e-9:
        y1 = float(np.interp(t1c, t, y))
        tt = np.append(tt, t1c)
        yy = np.append(yy, y1)

    return tt, yy


def _window_metrics(t: np.ndarray, y: np.ndarray, t0: float, t1: float) -> Dict[str, float]:
    """
    Metrics on y(t) within [t0, t1]:
      mean, peak, tpeak, auc_signed, auc_pos, auc_neg_abs
    """
    tt, yy = _ensure_boundary_points(t, y, t0, t1)
    if tt is None or yy is None:
        return dict(mean=np.nan, peak=np.nan, tpeak=np.nan,
                    auc_signed=np.nan, auc_pos=np.nan, auc_neg_abs=np.nan)

    dur = float(tt[-1] - tt[0])
    if dur <= 0:
        return dict(mean=np.nan, peak=np.nan, tpeak=np.nan,
                    auc_signed=np.nan, auc_pos=np.nan, auc_neg_abs=np.nan)

    auc_signed = float(np.trapezoid(yy, tt))
    auc_pos = float(np.trapezoid(np.maximum(yy, 0.0), tt))
    auc_neg_abs = float(np.trapezoid(np.abs(np.minimum(yy, 0.0)), tt))

    peak_idx = int(np.nanargmax(yy))
    peak = float(yy[peak_idx])
    tpeak = float(tt[peak_idx])

    mean = auc_signed / dur
    return dict(mean=mean, peak=peak, tpeak=tpeak,
                auc_signed=auc_signed, auc_pos=auc_pos, auc_neg_abs=auc_neg_abs)


def _value_at_time(t: np.ndarray, y: np.ndarray, t_query: float) -> float:
    """Linear interpolation at t_query; returns NaN if t_query is out of range."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if t.size == 0:
        return float("nan")
    order = np.argsort(t)
    t = t[order]
    y = y[order]
    if t_query < float(t[0]) or t_query > float(t[-1]):
        return float("nan")
    return float(np.interp(t_query, t, y))


# ----------------------------
# RAW (Wu PSC time courses)
# ----------------------------

def _raw_safe_end_seconds(g: pd.DataFrame) -> float:
    """
    Safe end time for raw curves:
      min(6.0, next_onset_med_s + 2.0) when available
      otherwise default by plot_group.
    """
    if "next_onset_med_s" in g.columns and pd.notna(g["next_onset_med_s"].iloc[0]):
        try:
            nxt = float(g["next_onset_med_s"].iloc[0])
            if math.isfinite(nxt) and nxt > 0:
                return float(min(6.0, nxt + 2.0))
        except Exception:
            pass

    pg = str(g["plot_group"].iloc[0])
    if pg == "mid_anticipation_primary":
        return 5.0
    return 6.0


def _compute_raw_long(df: pd.DataFrame) -> pd.DataFrame:
    required = ["sub", "ses", "run", "task", "echo", "cnfds", "plot_group", "condition", "time_s"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"RAW curves missing required column: {c}")

    if "mean" in df.columns:
        ycol = "mean"
    elif "psc" in df.columns:
        ycol = "psc"
    else:
        raise ValueError("RAW curves must contain a time-series column named 'mean' or 'psc'.")

    key_cols = ["sub", "ses", "run", "task", "echo", "cnfds", "plot_group", "condition"]

    out_rows: List[Dict[str, object]] = []
    for keys, g in df.groupby(key_cols, sort=False):
        t = g["time_s"].to_numpy(dtype=float)
        y = g[ycol].to_numpy(dtype=float)

        safe_end = _raw_safe_end_seconds(g)
        tmax = float(np.nanmax(t)) if t.size else float("nan")
        if math.isfinite(tmax):
            safe_end = float(min(safe_end, tmax))

        safe = _window_metrics(t, y, 0.0, safe_end)
        canon = _window_metrics(t, y, 4.0, 8.0)
        v6 = _value_at_time(t, y, 6.0)

        row = dict(zip(key_cols, keys))
        row.update(
            safe_end_s=safe_end,
            safe_mean=safe["mean"],
            safe_peak=safe["peak"],
            safe_tpeak=safe["tpeak"],
            safe_auc_signed=safe["auc_signed"],
            safe_auc_pos=safe["auc_pos"],
            safe_auc_neg_abs=safe["auc_neg_abs"],
            canon_mean=canon["mean"],
            canon_peak=canon["peak"],
            canon_tpeak=canon["tpeak"],
            canon_auc_signed=canon["auc_signed"],
            value_at_6s=v6,
        )
        out_rows.append(row)

    return pd.DataFrame(out_rows)


def _raw_contrast_defs(plot_group: str) -> List[Tuple[str, str]]:
    """List of (A, B) contrasts computed as A - B."""
    if plot_group == "mid_anticipation_primary":
        return [("Reward", "Neutral")]
    if plot_group.startswith("mid_feedback"):
        return [("Positive", "Negative")]
    if plot_group == "sr_decision_partner_extra":
        return [("Computer", "Stranger")]
    if plot_group in ("sr_outcome_primary", "sr_decision_by_future_outcome_extra"):
        return [("Reward", "Neutral"), ("Punish", "Neutral"), ("Reward", "Punish")]
    return []


def _compute_raw_contrasts(df: pd.DataFrame) -> pd.DataFrame:
    ycol = "mean" if "mean" in df.columns else "psc"
    base_cols = ["sub", "ses", "run", "task", "echo", "cnfds", "plot_group"]
    cond_col = "condition"

    out_rows: List[Dict[str, object]] = []

    for keys, g in df.groupby(base_cols, sort=False):
        plot_group = str(g["plot_group"].iloc[0])
        cdefs = _raw_contrast_defs(plot_group)
        if not cdefs:
            continue

        safe_end = _raw_safe_end_seconds(g)
        t_all = g["time_s"].to_numpy(dtype=float)
        tmax = float(np.nanmax(t_all)) if t_all.size else float("nan")
        if math.isfinite(tmax):
            safe_end = float(min(safe_end, tmax))

        cond_map: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for cond, gg in g.groupby(cond_col, sort=False):
            tt = gg["time_s"].to_numpy(dtype=float)
            yy = gg[ycol].to_numpy(dtype=float)
            cond_map[str(cond)] = (tt, yy)

        for a, b in cdefs:
            if a not in cond_map or b not in cond_map:
                continue
            ta, ya = cond_map[a]
            tb, yb = cond_map[b]

            t_union = np.unique(np.concatenate([ta, tb]))
            ya_u = np.interp(t_union, np.sort(ta), ya[np.argsort(ta)])
            yb_u = np.interp(t_union, np.sort(tb), yb[np.argsort(tb)])
            yd = ya_u - yb_u

            safe = _window_metrics(t_union, yd, 0.0, safe_end)
            canon = _window_metrics(t_union, yd, 4.0, 8.0)
            v6 = _value_at_time(t_union, yd, 6.0)

            row = dict(zip(base_cols, keys))
            row["contrast"] = f"{a}-{b}"
            row.update(
                safe_end_s=safe_end,
                safe_mean=safe["mean"],
                safe_peak=safe["peak"],
                safe_tpeak=safe["tpeak"],
                safe_auc_signed=safe["auc_signed"],
                safe_auc_pos=safe["auc_pos"],
                safe_auc_neg_abs=safe["auc_neg_abs"],
                canon_mean=canon["mean"],
                canon_peak=canon["peak"],
                canon_tpeak=canon["tpeak"],
                canon_auc_signed=canon["auc_signed"],
                value_at_6s=v6,
            )
            out_rows.append(row)

    return pd.DataFrame(out_rows)


# ----------------------------
# FIR (beta time courses)
# ----------------------------

def _normalize_fir_df(df: pd.DataFrame, *, tr_s: float) -> pd.DataFrame:
    """
    Accept either schema:
      A) time_s + amplitude + cnfds + roi
      B) bin + value + confounds (no cnfds/roi/time_s/amplitude)

    Returns a df containing, at minimum:
      feat, sub, ses, task, run, echo, cnfds, space, roi, time_s, condition, amplitude
    """
    dfx = df.copy()

    # feat column
    if "feat_dir" in dfx.columns and "feat" not in dfx.columns:
        dfx = dfx.rename(columns={"feat_dir": "feat"})
    if "feat" not in dfx.columns:
        dfx["feat"] = ""

    # cnfds column (your FIR file uses "confounds")
    if "cnfds" not in dfx.columns:
        if "confounds" in dfx.columns:
            dfx["cnfds"] = dfx["confounds"].astype(str)
        else:
            raise ValueError("FIR file missing required column: cnfds (or confounds)")

    # roi column (your FIR file is single-ROI)
    if "roi" not in dfx.columns:
        dfx["roi"] = "NAcc"

    # amplitude column (your FIR file uses "value")
    if "amplitude" not in dfx.columns:
        if "value" in dfx.columns:
            dfx["amplitude"] = pd.to_numeric(dfx["value"], errors="coerce")
        else:
            raise ValueError("FIR file missing required column: amplitude (or value)")

    # time_s column (your FIR file uses "bin")
    if "time_s" not in dfx.columns:
        if "bin" in dfx.columns:
            b = pd.to_numeric(dfx["bin"], errors="coerce")
            dfx["time_s"] = (b.astype(float) - 1.0) * float(tr_s)
        else:
            raise ValueError("FIR file missing required column: time_s (or bin)")

    return dfx


def _fir_plot_group_and_condition(task: str, condition: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Map FIR condition names to:
      - plot_group
      - condition_out (collapsed label)
    """
    task = str(task)
    condition = str(condition)

    if task == "mid":
        if condition in ("ant_rew", "ant_neu"):
            return "mid_anticipation_primary", ("Reward" if condition == "ant_rew" else "Neutral")
        if condition.startswith("out_pos"):
            return "mid_feedback_extra", "Positive"
        if condition.startswith("out_neg"):
            return "mid_feedback_extra", "Negative"
        return None, None

    if task == "sharedreward":
        if condition in ("computer", "stranger"):
            return "sr_decision_partner_extra", ("Computer" if condition == "computer" else "Stranger")
        if condition.startswith("rew_"):
            return "sr_outcome_primary", "Reward"
        if condition.startswith("neu_"):
            return "sr_outcome_primary", "Neutral"
        if condition.startswith("pun_"):
            return "sr_outcome_primary", "Punish"
        return None, None

    return None, None


def _fir_safe_end_seconds(plot_group: str) -> float:
    if plot_group == "mid_anticipation_primary":
        return 5.0
    return 6.0


def _compute_fir_long(df: pd.DataFrame, *, tr_s: float) -> pd.DataFrame:
    df = _normalize_fir_df(df, tr_s=tr_s)

    required = ["feat", "sub", "ses", "task", "run", "echo", "cnfds", "space", "roi", "time_s", "condition", "amplitude"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"FIR file missing required column: {c}")

    pg_list: List[Optional[str]] = []
    co_list: List[Optional[str]] = []
    for task, cond in zip(df["task"].astype(str), df["condition"].astype(str)):
        pg, co = _fir_plot_group_and_condition(task, cond)
        pg_list.append(pg)
        co_list.append(co)

    dfx = df.copy()
    dfx["plot_group"] = pg_list
    dfx["condition_out"] = co_list
    dfx = dfx[dfx["plot_group"].notna() & dfx["condition_out"].notna()].copy()

    dfx["confounds"] = dfx["cnfds"].astype(str)

    key_cols = ["feat", "sub", "ses", "task", "run", "echo", "confounds",
                "space", "roi", "plot_group", "condition_out"]

    out_rows: List[Dict[str, object]] = []
    for keys, g in dfx.groupby(key_cols, sort=False):
        t = g["time_s"].to_numpy(dtype=float)
        y = g["amplitude"].to_numpy(dtype=float)

        plot_group = str(g["plot_group"].iloc[0])
        safe_end = _fir_safe_end_seconds(plot_group)
        tmax = float(np.nanmax(t)) if t.size else float("nan")
        if math.isfinite(tmax):
            safe_end = float(min(safe_end, tmax))

        safe = _window_metrics(t, y, 0.0, safe_end)
        canon = _window_metrics(t, y, 4.0, 8.0)
        v6 = _value_at_time(t, y, 6.0)

        row = dict(zip(key_cols, keys))
        row["condition"] = row.pop("condition_out")
        row.update(
            safe_end_s=safe_end,
            safe_mean=safe["mean"],
            safe_peak=safe["peak"],
            safe_tpeak=safe["tpeak"],
            safe_auc_signed=safe["auc_signed"],
            safe_auc_pos=safe["auc_pos"],
            safe_auc_neg_abs=safe["auc_neg_abs"],
            canon_mean=canon["mean"],
            canon_peak=canon["peak"],
            canon_tpeak=canon["tpeak"],
            canon_auc_signed=canon["auc_signed"],
            value_at_6s=v6,
        )
        out_rows.append(row)

    return pd.DataFrame(out_rows)


def _fir_contrast_defs(plot_group: str) -> List[Tuple[str, str]]:
    if plot_group == "mid_anticipation_primary":
        return [("Reward", "Neutral")]
    if plot_group.startswith("mid_feedback"):
        return [("Positive", "Negative")]
    if plot_group == "sr_decision_partner_extra":
        return [("Computer", "Stranger")]
    if plot_group == "sr_outcome_primary":
        return [("Reward", "Neutral"), ("Punish", "Neutral"), ("Reward", "Punish")]
    return []


def _compute_fir_contrasts(df: pd.DataFrame, *, tr_s: float) -> pd.DataFrame:
    df = _normalize_fir_df(df, tr_s=tr_s)

    required = ["feat", "sub", "ses", "task", "run", "echo", "cnfds", "space", "roi", "time_s", "condition", "amplitude"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"FIR file missing required column: {c}")

    pg_list: List[Optional[str]] = []
    co_list: List[Optional[str]] = []
    for task, cond in zip(df["task"].astype(str), df["condition"].astype(str)):
        pg, co = _fir_plot_group_and_condition(task, cond)
        pg_list.append(pg)
        co_list.append(co)

    dfx = df.copy()
    dfx["plot_group"] = pg_list
    dfx["condition_out"] = co_list
    dfx = dfx[dfx["plot_group"].notna() & dfx["condition_out"].notna()].copy()

    dfx["confounds"] = dfx["cnfds"].astype(str)

    base_cols = ["feat", "sub", "ses", "task", "run", "echo", "confounds",
                 "space", "roi", "plot_group"]

    out_rows: List[Dict[str, object]] = []
    for keys, g in dfx.groupby(base_cols, sort=False):
        plot_group = str(g["plot_group"].iloc[0])
        cdefs = _fir_contrast_defs(plot_group)
        if not cdefs:
            continue

        safe_end = _fir_safe_end_seconds(plot_group)
        t_all = g["time_s"].to_numpy(dtype=float)
        tmax = float(np.nanmax(t_all)) if t_all.size else float("nan")
        if math.isfinite(tmax):
            safe_end = float(min(safe_end, tmax))

        cond_map: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for cond, gg in g.groupby("condition_out", sort=False):
            tt = gg["time_s"].to_numpy(dtype=float)
            yy = gg["amplitude"].to_numpy(dtype=float)
            cond_map[str(cond)] = (tt, yy)

        for a, b in cdefs:
            if a not in cond_map or b not in cond_map:
                continue
            ta, ya = cond_map[a]
            tb, yb = cond_map[b]

            t_union = np.unique(np.concatenate([ta, tb]))
            ya_u = np.interp(t_union, np.sort(ta), ya[np.argsort(ta)])
            yb_u = np.interp(t_union, np.sort(tb), yb[np.argsort(tb)])
            yd = ya_u - yb_u

            safe = _window_metrics(t_union, yd, 0.0, safe_end)
            canon = _window_metrics(t_union, yd, 4.0, 8.0)
            v6 = _value_at_time(t_union, yd, 6.0)

            row = dict(zip(base_cols, keys))
            row["contrast"] = f"{a}-{b}"
            row.update(
                safe_end_s=safe_end,
                safe_mean=safe["mean"],
                safe_peak=safe["peak"],
                safe_tpeak=safe["tpeak"],
                safe_auc_signed=safe["auc_signed"],
                safe_auc_pos=safe["auc_pos"],
                safe_auc_neg_abs=safe["auc_neg_abs"],
                canon_mean=canon["mean"],
                canon_peak=canon["peak"],
                canon_tpeak=canon["tpeak"],
                canon_auc_signed=canon["auc_signed"],
                value_at_6s=v6,
            )
            out_rows.append(row)

    return pd.DataFrame(out_rows)


# ----------------------------
# I/O + CLI
# ----------------------------

def _read_csv_auto(path: Path) -> pd.DataFrame:
    if path.suffix == ".gz":
        return pd.read_csv(path, compression="gzip")
    return pd.read_csv(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-curves", default="run_curves_long_wu.csv.gz",
                    help="Raw Wu-normalized curves file")
    ap.add_argument("--fir-curves", default="fir_copes_roi_means_long.csv",
                    help="FIR copes ROI means long file")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--tr-s", type=float, default=1.615,
                    help="TR in seconds (used to convert FIR bin -> time_s when needed)")
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    raw_path = Path(args.raw_curves)
    if raw_path.exists():
        df_raw = _read_csv_auto(raw_path)
        raw_long = _compute_raw_long(df_raw)
        raw_con = _compute_raw_contrasts(df_raw)
        raw_long.to_csv(outdir / "metrics_raw_long.csv", index=False)
        raw_con.to_csv(outdir / "metrics_raw_contrasts.csv", index=False)
    else:
        print(f"[raw] not found: {raw_path} (skipping)")

    fir_path = Path(args.fir_curves)
    if fir_path.exists():
        df_fir = _read_csv_auto(fir_path)
        fir_long = _compute_fir_long(df_fir, tr_s=args.tr_s)
        fir_con = _compute_fir_contrasts(df_fir, tr_s=args.tr_s)
        fir_long.to_csv(outdir / "metrics_fir_long.csv", index=False)
        fir_con.to_csv(outdir / "metrics_fir_contrasts.csv", index=False)
    else:
        print(f"[fir] not found: {fir_path} (skipping)")


if __name__ == "__main__":
    main()
