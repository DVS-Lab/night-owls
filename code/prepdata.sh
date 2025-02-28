#/usr/bin/env bash

# Example code for heudiconv and pydeface. This will get your data ready for analyses.
# This code will convert DICOMS to BIDS (PART 1). Will also deface (PART 2) and run MRIQC (PART 3).

# usage: bash prepdata.sh sub ses
# example: bash prepdata.sh 104 01

# Notes:
# 1) containers live under /data/tools on local computer. should these relative paths and shared? YODA principles would suggest so.
# 2) aside from containers, only absolute path in whole workflow (transparent to folks who aren't allowed to access to raw data)
sourcedata=/ZPOOL/data/sourcedata/sourcedata/NOSC

sub=$1
ses=$2


# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
dsroot="$(dirname "$scriptdir")"

echo ${dsroot}

# make bids folder if it doesn't exist
if [ ! -d $dsroot/bids ]; then
	mkdir -p $dsroot/bids
fi

# overwrite existing
rm -rf $dsroot/bids/sub-${sub}/ses-${ses}


# PART 1: running heudiconv and fixing fieldmaps
apptainer run --cleanenv \
-B $dsroot:/out \
-B $sourcedata:/sourcedata \
/ZPOOL/data/tools/heudiconv-1.3.2.sif \
-d /sourcedata/Smith-NOSC-{subject}-SES{session}/*/scans/*/*/DICOM/files/*.dcm \
-o /out/bids/ \
-f /out/code/heuristics.py \
-s $sub \
-ss $ses \
-c dcm2niix \
-b --minmeta --overwrite


## PART 2: Defacing anatomicals and date shifting to ensure compatibility with data sharing. (do we really need to shift the dates? ask IRB?)

## note that you may need to install pydeface via pip or conda
bidsroot=$dsroot/bids
echo "defacing subject $sub $ses"
pydeface ${bidsroot}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_T1w.nii.gz
mv -f ${bidsroot}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_T1w_defaced.nii.gz ${bidsroot}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_T1w.nii.gz

## shift dates on scans to reduce likelihood of re-identification
python $scriptdir/shiftdates.py $dsroot/bids/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_scans.tsv
