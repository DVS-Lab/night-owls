import json
import os

# Define BIDS root and func directory name
bidsdir = "/gpfs/scratch/tug87422/smithlab-shared/night-owls/bids"
func_dir = "func"

# Mapping from acq-* in fieldmap files to task-* in functional files
acq_to_task = {
    "rest": "rest",
    "sharedreward": "sharedreward",
    "mid": "mid"
}

# Loop over subjects
subs = [d for d in os.listdir(bidsdir) if os.path.isdir(os.path.join(bidsdir, d)) and d.startswith('sub')]

for subj in subs:
    subj_path = os.path.join(bidsdir, subj)

    # Loop over sessions
    sessions = [s for s in os.listdir(subj_path) if os.path.isdir(os.path.join(subj_path, s)) and s.startswith('ses')]
    
    for ses in sessions:
        print(f"Running subject: {subj}, session: {ses}")
        fmap_dir = os.path.join(subj_path, ses, 'fmap')
        if not os.path.isdir(fmap_dir):
            continue

        json_files = [f for f in os.listdir(fmap_dir) if f.endswith('fieldmap.json') or f.endswith('magnitude.json')]

        for json_file in json_files:
            json_path = os.path.join(fmap_dir, json_file)
            with open(json_path, 'r') as f:
                data = json.load(f)
                intended_for = []

                # Extract acq and run from fieldmap filename
                file_parts = json_file.split('_')
                acq = next((p.split('-')[1] for p in file_parts if p.startswith('acq-')), None)
                run = next((p.split('-')[1] for p in file_parts if p.startswith('run-')), None)

                if not acq or not run:
                    print(f"Could not extract acq/run from: {json_file}")
                    continue

                # Map acquisition label to task name
                task = acq_to_task.get(acq)
                if not task:
                    print(f"Unknown acquisition type '{acq}' in {json_file}, skipping.")
                    continue

                # Construct intended file paths for 4 echoes
                for echo in range(1, 5):
                    intended_for.append(
                        f"{subj}/{ses}/{func_dir}/{subj}_{ses}_task-{task}_run-{run}_echo-{echo}_part-mag_bold.nii.gz"
                    )

                data["IntendedFor"] = intended_for
                data["Units"] = "Hz"
                data.pop("EchoTime1", None)
                data.pop("EchoTime2", None)

            with open(json_path, 'w') as f:
                json.dump(data, f, indent=4, sort_keys=True)

            print(f"✔ Added IntendedFor to {json_file}")

