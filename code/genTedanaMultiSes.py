#!/usr/bin/env python

import os
import re
import numpy as np
import pandas as pd
from natsort import natsorted

# adjust these to your actual base paths
BASE      = "/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives"
TEDANA    = os.path.join(BASE, "tedana")
FMRIPREP  = os.path.join(BASE, "fmriprep")
FSL_OUT   = os.path.join(BASE, "fsl", "confounds_tedana") #This needs to account for /sub/ses

# find all tedana_metrics files under sub-XX/ses-YY
metric_files = natsorted(
  os.path.join(root, f)
    for root, _, files in os.walk(TEDANA)
    for f in files
    if f.endswith("tedana_metrics.tsv")
)

for mfile in metric_files:
    # root = /…/tedana/sub-101/ses-01/
    root = os.path.dirname(mfile)
    rel  = os.path.relpath(root, TEDANA).split(os.sep)
    sub  = rel[0]                  # e.g. "sub-101"
    ses  = rel[1]                  # e.g. "ses-01"

    # extract run and task from filename
    run   = re.search(r"run-(\d+)_desc-tedana", mfile).group(1)
    task  = re.search(r"_task-([^_]+)_",      mfile).group(1)

    # build path to fmriprep confounds
    fprep_file = os.path.join(
      FMRIPREP,
      sub,
      ses,
      "func",
      f"{sub}_{ses}_task-{task}_run-{run}_part-mag_desc-confounds_timeseries.tsv"
    )

    if not os.path.exists(fprep_file):
        print(f"fmriprep missing for {sub} {ses} run-{run} task-{task}")
        continue

    print(f"Processing confounds for {sub} {ses} run-{run}")

    # load data
    fprep_df   = pd.read_csv(fprep_file, sep="\t")
    mixing_file = os.path.join(
    root,
    f"{sub}_{ses}_task-{task}_run-{run}_desc-ICA_mixing.tsv"
    )
    mixing_df = pd.read_csv(mixing_file, sep="\t")    metrics_df = pd.read_csv(os.path.join(root, "tedana_metrics.tsv"), sep="\t")

    # pick out the rejected components
    bad_idxs       = metrics_df.loc[metrics_df["classification"] == "rejected", "Component"]
    bad_components = mixing_df.loc[:, bad_idxs]

    # select fmriprep columns
    aCompCor = [f"a_comp_cor_{i:02}" for i in range(6)]
    cosine   = [c for c in fprep_df if c.startswith("cosine")]
    nss      = [c for c in fprep_df if c.startswith("non_steady_state")]
    motion   = ["trans_x","trans_y","trans_z","rot_x","rot_y","rot_z"]
    fd       = ["framewise_displacement"]
    keep     = aCompCor + cosine + nss + motion + fd

    conf_df = fprep_df.loc[:, keep].fillna(0) \
      .join(bad_components, how="left")

    # write out, no header/index as FSL wants
    out_dir = os.path.join(FSL_OUT, sub, ses)
    os.makedirs(out_dir, exist_ok=True)
    out_fn  = os.path.join(
      out_dir,
      f"{sub}_{ses}_task-{task}_run-{run}_desc-TedanaPlusConfounds.tsv"
    )
    conf_df.to_csv(out_fn, sep="\t", index=False, header=False)
