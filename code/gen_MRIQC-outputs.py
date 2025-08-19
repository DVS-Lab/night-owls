#!/usr/bin/env python

import os, re, json
import pandas as pd

# Hard-coded input and output
mriqc_path = "/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/mriqc"  
out_file = "/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/data-outputs/mriqc_metrics.csv"

# Collect all *_bold.json files
j_files = []
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

# Sort so echoes are in order and 'avg' appears last
def echo_sort_key(e):
    if e == "avg":
        return float('inf')
    else:
        return int(e.split("-")[1])

# First sort by sub, ses, task, run
df_final = df_final.sort_values(by=["sub", "ses", "task", "run"])

# Then sort echoes within each run
df_final = df_final.iloc[df_final.groupby(["sub", "ses", "task", "run"]).echo.apply(
    lambda x: x.map(echo_sort_key).sort_values().index
).explode().astype(int)]

# Reorder columns
columns_order = ["sub", "ses", "task", "run", "echo", "image", "mean_fd", "tsnr"]
df_final = df_final[columns_order]

# Round numeric columns to 3 decimal places
df_final[["mean_fd", "tsnr"]] = df_final[["mean_fd", "tsnr"]].round(3)

# Save to CSV
df_final.to_csv(out_file, index=False)
print(f"✅ Wrote {len(df_final)} rows (including averages) to {out_file}")
