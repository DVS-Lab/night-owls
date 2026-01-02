#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # critical on login nodes
import matplotlib.pyplot as plt
import pandas as pd

# --------------------------- Config ---------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR   = SCRIPT_DIR.parent

CSV_PATH   = ROOT_DIR / "derivatives" / "extractions" / "fir" / "fir_copes_roi_means_long.csv"
PLOTS_DIR  = ROOT_DIR / "derivatives" / "extractions" / "fir" / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

DPI = 150

# Optional: nicer ordering for SharedReward conditions (fallback = alphabetical)
PREFERRED_COND_ORDER = [
    "computer", "stranger",
    "pun_c", "neu_c", "rew_c",
    "pun_s", "neu_s", "rew_s",
]


def cond_sort_key(cond: str) -> tuple[int, str]:
    if cond in PREFERRED_COND_ORDER:
        return (0, f"{PREFERRED_COND_ORDER.index(cond):02d}")
    return (1, cond)


def plot_timecourse(df: pd.DataFrame, title: str, outpath: Path) -> None:
    d = df.copy()
    d["bin"] = pd.to_numeric(d["bin"], errors="coerce")
    d = d.dropna(subset=["bin"])
    d["bin"] = d["bin"].astype(int)

    wide = (
        d.pivot_table(index="bin", columns="condition", values="value", aggfunc="mean")
         .sort_index()
    )
    if wide.empty:
        return

    # order conditions nicely when possible
    conds = sorted(list(wide.columns), key=cond_sort_key)
    wide = wide[conds]

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for cond in wide.columns:
        ax.plot(wide.index.values, wide[cond].values, label=str(cond))

    ax.set_title(title)
    ax.set_xlabel("FIR bin")
    ax.set_ylabel("VS ROI mean (COPE)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Missing CSV: {CSV_PATH}")

    df = pd.read_csv(CSV_PATH)

    needed = {"sub","ses","task","run","echo","confounds","condition","bin","value"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"CSV is missing expected columns: {sorted(missing)}")

    # ---------------- Session-level (avg across runs) ----------------
    ses_root = PLOTS_DIR / "session_level"
    ses_df = (
        df.groupby(["task","sub","ses","echo","confounds","condition","bin"], as_index=False)["value"]
          .mean()
    )

    for keys, g in ses_df.groupby(["task","sub","ses","echo","confounds"], sort=True):
        task, sub, ses, echo, conf = keys
        title = f"{task} | sub-{sub} ses-{ses} (avg runs) | {echo} | cnfds-{conf}"
        out = ses_root / task / f"sub-{sub}" / f"ses-{ses}_{echo}_cnfds-{conf}.png"
        plot_timecourse(g, title, out)

    # ---------------- Subject-level (avg across sessions+runs) ----------------
    sub_root = PLOTS_DIR / "subject_level"
    sub_df = (
        df.groupby(["task","sub","echo","confounds","condition","bin"], as_index=False)["value"]
          .mean()
    )

    for keys, g in sub_df.groupby(["task","sub","echo","confounds"], sort=True):
        task, sub, echo, conf = keys
        title = f"{task} | sub-{sub} (avg sessions+runs) | {echo} | cnfds-{conf}"
        out = sub_root / task / f"sub-{sub}_{echo}_cnfds-{conf}.png"
        plot_timecourse(g, title, out)

    print(f"Plots written under: {PLOTS_DIR}")


if __name__ == "__main__":
    main()
