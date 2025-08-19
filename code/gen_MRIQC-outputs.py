#!/usr/bin/env python

import os, re, json
import pandas as pd

# Hard-coded input and output
mriqc_path = "/ZPOOL/data/projects/rf1-sra-data/bids/derivatives/mriqc"  
out_file = "/ZPOOL/data/projects/rf1-sra-data/bids/derivatives/mriqc_metrics.csv"

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
