#!/usr/bin/env python

import os, re, json
import pandas as pd

# Hard-coded input and output
mriqc_path_shared = "/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/mriqc"
mriqc_path_personal = "/home/tun47039/mriqc_sub103/derivatives/mriqc"
out_file = "/gpfs/scratch/tun47039/night-owls/derivatives/data-outputs/mriqc/mriqc_metrics.csv"

# Create output directory if it doesn't exist
os.makedirs(os.path.dirname(out_file), exist_ok=True)

# Collect all *_bold.json files from both locations
j_files = []
for mriqc_path in [mriqc_path_shared, mriqc_path_personal]:
    if os.path.exists(mriqc_path):
        for root, dirs, files in os.walk(mriqc_path):
            for f in files:
                if f.endswith("bold.json"):
                    j_files.append(os.path.join(root, f))

rows = []
for j in j_files:
    with open(j) as f:
        data = json.load(f)

    fname = os.path.basename(j)

    # Extract BIDS entities from filename
    try:
        sub  = re.search(r"(sub-[^_]+)", fname).group(1)
        ses  = re.search(r"(ses-[^_]+)", fname).group(1)
        task = re.search(r"(task-[^_]+)", fname).group(1)
        run  = re.search(r"(run-[^_]+)", fname).group(1)
        echo = re.search(r"(echo-[^_]+)", fname).group(1)
    except AttributeError:
        print(f"⚠️ Skipping unmatched file: {fname}")
        continue

    # Image ID = through run (ignores echo/part/etc.)
    image_id = f"{sub}_{ses}_{task}_{run}"

    rows.append({
        "image": image_id,
        "sub": sub,
        "ses": ses,
        "task": task,
        "run": run,
        "echo": echo,
        "mean_fd": data.get("fd_mean", None),
        "tsnr": data.get("tsnr", None)
    })

# Convert to dataframe
df = pd.DataFrame(rows)
print(f"✅ Parsed {len(df)} JSON files into dataframe")

if df.empty:
    raise RuntimeError("No files were parsed. Check regex patterns and filenames.")

# Compute averages across echoes
avg_df = (
    df.groupby(["image", "sub", "ses", "task", "run"], as_index=False)[["mean_fd", "tsnr"]]
      .mean()
)
avg_df["echo"] = "avg"

# Combine original + averages
df_final = pd.concat([df, avg_df], ignore_index=True)

# Simple sort by sub → ses → task → run → echo
df_final = df_final.sort_values(by=["sub", "ses", "task", "run", "echo"])

# Reorder columns
columns_order = ["sub", "ses", "task", "run", "echo", "image", "mean_fd", "tsnr"]
df_final = df_final[columns_order]

# Round numeric columns to 3 decimal places
df_final[["mean_fd", "tsnr"]] = df_final[["mean_fd", "tsnr"]].round(3)

# Save to CSV
df_final.to_csv(out_file, index=False)
print(f"✅ Summary CSV saved: {out_file}")

# --- Create per-subject per-task QC files (only echo-2) ---
for (sub, task), group in df_final.groupby(["sub", "task"]):
    # Keep only echo-2 rows
    qc_group = group[group["echo"] == "echo-2"].copy()
    if qc_group.empty:
        print(f"⚠️ No echo-2 data for {sub} {task}, skipping QC CSV")
        continue

    # Create the row label: sub_ses-run
    qc_group["row_label"] = qc_group["sub"] + "_" + qc_group["ses"] + "-" + qc_group["run"]
    qc_df = qc_group[["row_label", "mean_fd", "tsnr"]]
    
    qc_file = os.path.join(
        os.path.dirname(out_file),
        f"{sub}_{task}_qc.csv"
    )
    qc_df.to_csv(qc_file, index=False)
    print(f"✅ QC CSV saved: {qc_file}")