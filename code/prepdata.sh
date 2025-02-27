#/usr/bin/env bash

# Example code for heudiconv and pydeface. This will get your data ready for analyses.
# This code will convert DICOMS to BIDS (PART 1). Will also deface (PART 2) and run MRIQC (PART 3).

# usage: bash prepdata.sh sub nruns
# example: bash prepdata.sh 104 3

# Notes:
# 1) containers live under /data/tools on local computer. should these relative paths and shared? YODA principles would suggest so.
# 2) other projects should use Jeff's python script for fixing the IntendedFor
# 3) aside from containers, only absolute path in whole workflow (transparent to folks who aren't allowed to access to raw data)
sourcedata=/ZPOOL/data/sourcedata/sourcedata/rf1-sra
sub=$1

except_subs=(1001)
for i in "${except_subs[@]}" ; do
    if [ "$i" -eq "$sub" ] ; then
        echo "Exception ${sub}"
	      exit 1
    fi
done


# ensure paths are correct irrespective from where user runs the script
codedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
dsroot="$(dirname "$scriptdir")"

echo ${dsroot}

# make bids folder if it doesn't exist
if [ ! -d $dsroot/bids ]; then
	mkdir -p $dsroot/bids
fi

# overwrite existing
rm -rf $dsroot/bids/sub-${sub}


# PART 1: running heudiconv and fixing fieldmaps
apptainer run --cleanenv \
-B $dsroot:/out \
-B $sourcedata:/sourcedata \
/ZPOOL/data/tools/heudiconv-1.3.2.sif \
-d /sourcedata/Smith-SRA-{subject}/*/scans/*/*/DICOM/files/*.dcm \
-o /out/bids/ \
-f /out/code/heuristics_rf1.py \
-s $sub \
-c dcm2niix \
-b --minmeta --overwrite


## PART 2: Defacing anatomicals and date shifting to ensure compatibility with data sharing.

## note that you may need to install pydeface via pip or conda
bidsroot=$dsroot/bids
echo "defacing subject $sub"
pydeface ${bidsroot}/sub-${sub}/anat/sub-${sub}_T1w.nii.gz
mv -f ${bidsroot}/sub-${sub}/anat/sub-${sub}_T1w_defaced.nii.gz ${bidsroot}/sub-${sub}/anat/sub-${sub}_T1w.nii.gz
pydeface ${bidsroot}/sub-${sub}/anat/sub-${sub}_FLAIR.nii.gz
mv -f ${bidsroot}/sub-${sub}/anat/sub-${sub}_FLAIR_defaced.nii.gz ${bidsroot}/sub-${sub}/anat/sub-${sub}_FLAIR.nii.gz

## shift dates on scans to reduce likelihood of re-identification
python $codedir/shiftdates.py $dsroot/bids/sub-${sub}/sub-${sub}_scans.tsv
