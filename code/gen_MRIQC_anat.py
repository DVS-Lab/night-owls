#!/usr/bin/env python

import os, re, json
import pandas as pd

# Hard-coded input and output
mriqc_path_shared = "/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/mriqc"
mriqc_path_personal = "/home/tun47039/mriqc_sub103/derivatives/mriqc"
out_file = "/gpfs/scratch/tun47039/night-owls/derivatives/data-outputs/mriqc/mriqc_metrics_anat.csv"

# Create output directory if it doesn't exist
os.makedirs(os.path.dirname(out_file), exist_ok=True)

# Collect all *_T1w.json files from both locations
j_files = []
for mriqc_path in [mriqc_path_shared, mriqc_path_personal]:
    if os.path.exists(mriqc_path):
        for root, dirs, files in os.walk(mriqc_path):
            for f in files:
                if f.endswith("T1w.json"):
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
    except AttributeError:
        print(f"⚠️ Skipping unmatched file: {fname}")
        continue

    # Image ID = sub_ses only
    image_id = f"{sub}_{ses}"

    rows.append({
        "image": image_id,
        "sub": sub,
        "ses": ses,
        "cnr": data.get("cnr", None),
        "snr_total": data.get("snr_total", None),
        "efc": data.get("efc", None),
        "fber": data.get("fber", None),
        "cjv": data.get("cjv", None),
        "qi_1": data.get("qi_1", None)
    })

# Convert to dataframe
df = pd.DataFrame(rows)
print(f"✅ Parsed {len(df)} JSON files into dataframe")

if df.empty:
    raise RuntimeError("No files were parsed. Check regex patterns and filenames.")

# Simple sort by sub → ses
df = df.sort_values(by=["sub", "ses"])

# Reorder columns
columns_order = ["sub", "ses", "image", "cnr", "snr_total", "efc", "fber", "cjv", "qi_1"]
df = df[columns_order]

# Round numeric columns to 3 decimal places
df[["cnr", "snr_total", "efc", "fber", "cjv", "qi_1"]] = df[["cnr", "snr_total", "efc", "fber", "cjv", "qi_1"]].round(3)

# Save to CSV
df.to_csv(out_file, index=False)
print(f"✅ Summary CSV saved: {out_file}")

# --- Create per-subject QC files ---
for sub, group in df.groupby("sub"):
    # Create the row label: sub_ses
    qc_group = group.copy()
    qc_group["row_label"] = qc_group["sub"] + "_" + qc_group["ses"]
    qc_df = qc_group[["row_label", "cnr", "snr_total", "efc", "fber", "cjv", "qi_1"]]
    
    qc_file = os.path.join(
        os.path.dirname(out_file),
        f"{sub}_qc.csv"
    )
    qc_df.to_csv(qc_file, index=False)
    print(f"✅ QC CSV saved: {qc_file}")