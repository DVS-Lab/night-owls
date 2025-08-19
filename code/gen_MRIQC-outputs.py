import os, re, json
import pandas as pd

# Hard-coded input and output
mriqc_path = "/ZPOOL/data/projects/night-owls/bids/derivatives/mriqc"  
out_file = "/ZPOOL/data/projects/night-owls/bids/derivatives/data-outputs/mriqc_metrics.csv"

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

    # Example: sub-104_ses-01_task-mid_run-2_echo-3_part-mag_bold.json
    fname = os.path.basename(j)

    # Regex to capture subject, session, task, run, echo
    match = re.match(
        r"(sub-[^_]+)_"
        r"(ses-[^_]+)_"
        r"(task-[^_]+)_"
        r"(run-[^_]+)_"
        r"(echo-[^_]+)_", 
        fname
    )
    if not match:
        continue

    sub, ses, task, run, echo = match.groups()
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

# Compute averages across echoes
avg_df = (
    df.groupby(["image", "sub", "ses", "task", "run"], as_index=False)[["mean_fd", "tsnr"]]
      .mean()
)
avg_df["echo"] = "avg"

# Combine original + averages
df_final = pd.concat([df, avg_df], ignore_index=True)

# Save to CSV
df_final.to_csv(out_file, index=False)
print(f"✅ Wrote {len(df_final)} rows (including averages) to {out_file}")
