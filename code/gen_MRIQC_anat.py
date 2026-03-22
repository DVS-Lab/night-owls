#!/usr/bin/env python

import os, re, json
import pandas as pd

# Hard-coded input and output
mriqc_path_shared = "/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/mriqc"
mriqc_path_personal = "/home/tun47039/mriqc_sub103/derivatives/mriqc"
out_file = "/gpfs/scratch/tun47039/night-owls/derivatives/data-outputs/mriqc/mriqc_metrics_anat.csv"

# Create output directory if it doesn't