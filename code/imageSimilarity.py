import itertools
import pandas as pd
from pathlib import Path
from pyrelimri.similarity import pairwise_similarity

# -----------------------------
# Parameters
# -----------------------------
subs = [101, 103, 104, 105]
tasks = ['mid', 'sharedreward']
spaces = ['mni', 't1w']
echoes = ['multi-echo', 'single-echo']
confounds = ['cnfds-tedana', 'cnfds-fmriprep']
runs = [1, 2]

data_root = Path("/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/fsl")

#find all runs matching combination
def get_files_for_combination(sub, task, space, echo, confound):
    files = []
    for run in runs:
        # Determine cope number
        cope_n = 7 if (task == 'mid') else 11

        # Find all sessions for this subject
        ses_dirs = sorted(data_root.glob(f"sub-{sub}/ses-*/L1_sub-{sub}_ses-*_task-{task}_model-1_type-act_run-{run}_space-{space}_{echo}_{confound}.feat/stats/cope{cope_n}.nii.gz"))
        files.extend(ses_dirs)
    return [str(f) for f in files]

#similarity loop
all_combinations = list(itertools.product(tasks, spaces, echoes, confounds, runs))
results = []

for comb in all_combinations:
    task, space, echo, confound, run = comb
    row = {'Combination': f"{task}_{space}_{echo}_{confound}_run{run}"}

    for sub in subs:
        files = get_files_for_combination(sub, task, space, echo, confound)
        # Compute pairwise Spearman similarity
        df = pairwise_similarity(files, similarity_type='spearman')
        # take the mean of all pairwise correlations
        row[f"sub-{sub}"] = df['similar_coef'].mean()
    results.append(row)

#export
df_out = pd.DataFrame(results)
df_out.to_csv("/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/data-outputs/pairwise_icc_results.csv", index=False)
