import json
import os

# Test directory
playdir = "/ZPOOL/data/projects/night-owls/"

# Only run for sub-test
subj = "sub-test"
subj_path = os.path.join(playdir, subj)

print("Running subject:", subj)

# Walk through ses-xx subdirs
for root, _, files in os.walk(subj_path):
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
