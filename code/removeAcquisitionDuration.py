import json
import os

# Define paths
bidsdir = "/ZPOOL/data/projects/night-owls/bids/"

# Find all subject directories in the BIDS directory
subs = [d for d in os.listdir(bidsdir) if os.path.isdir(os.path.join(bidsdir, d)) and d.startswith('sub')]

for subj in subs:
    print("Running subject:", subj)

    for root, _, files in os.walk(os.path.join(bidsdir, subj)):
        if root.endswith('func'):
            json_files = [f for f in files if f.endswith('bold.json')]

            for json_file in json_files:
                json_path = os.path.join(root, json_file)

                with open(json_path, 'r') as f:
                    data = json.load(f)

                # Remove AcquisitionDuration if present
                if "AcquisitionDuration" in data:
                    del data["AcquisitionDuration"]

                    with open(json_path, 'w') as f:
                        json.dump(data, f, indent=4, sort_keys=True)

                    print("Removed AcquisitionDuration from", json_path)
                else:
                    print("No AcquisitionDuration in", json_path)
