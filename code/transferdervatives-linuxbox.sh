#!/bin/bash

# Ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

read -p "Enter AccessNet ID: " destination_user

source_directory="/gpfs/scratch/tug87422/smithlab-shared/night-owls/derivatives/ses-01"
destination_server="@cla18994.tu.temple.edu:"
destination_path="/ZPOOL/data/projects/night-owls/derivatives/fmriprep"

for sub in `cat ${basedir}/code/sublist.txt`; do
	der_files="$source_directory/sub-$sub"
	html_files="$source_directory/sub-$sub.html"
	rsync -avh --no-compress --progress "$der_files" "$destination_user""$destination_server""$destination_path"
	rsync -avh --no-compress --progress "$html_files" "$destination_user""$destination_server""$destination_path"
done

exit   