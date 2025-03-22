import os

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
        if template is None or not template:
                raise ValueError('Template must be a valid format string')
        return template, outtype, annotation_classes

def infotodict(seqinfo):
    t1w = create_key('sub-{subject}/anat/sub-{subject}_T1w')
    mag = create_key('sub-{subject}/fmap/sub-{subject}_acq-bold_magnitude')
    phase = create_key('sub-{subject}/fmap/sub-{subject}_acq-bold_phasediff')
    t2_flair = create_key('sub-{subject}/anat/sub-{subject}_FLAIR')
    sharedreward_mag = create_key('sub-{subject}/func/sub-{subject}_task-sharedreward_run-{item:d}_part-mag_bold')
    sharedreward_phase = create_key('sub-{subject}/func/sub-{subject}_task-sharedreward_run-{item:d}_part-phase_bold')
    sharedreward_sbref = create_key('sub-{subject}/func/sub-{subject}_task-sharedreward_run-{item:d}_sbref')
    mid_mag = create_key('sub-{subject}/func/sub-{subject}_task-mid_run-{item:d}_part-mag_bold')
    mid_phase = create_key('sub-{subject}/func/sub-{subject}_task-mid_run-{item:d}_part-phase_bold')
    mid_sbref = create_key('sub-{subject}/func/sub-{subject}_task-mid_run-{item:d}_sbref')
    dwi = create_key('sub-{subject}/dwi/sub-{subject}_dwi')
    dwi_pa = create_key('sub-{subject}/fmap/sub-{subject}_acq-dwi_dir-PA_epi')
    dwi_ap = create_key('sub-{subject}/fmap/sub-{subject}_acq-dwi_dir-AP_epi')

    info = {t1w: [],
            mag: [], phase: [],
            dwi: [], dwi_pa: [], dwi_ap: [],
            t2_flair: [],
            sharedreward_mag: [], sharedreward_phase: [], sharedreward_sbref: [],
            mid_mag: [], mid_phase: [], mid_sbref: []}

    list_of_ids = [s.series_id for s in seqinfo]

    for s in seqinfo:

        # anatomicals and standard fmaps
        if ('T1w-anat_mpg_07sag_iso' in s.protocol_name) and ('NORM' in s.image_type):
            info[t1w] = [s.series_id]
        if ('gre_field' in s.protocol_name) and ('NORM' in s.image_type):
            info[mag] = [s.series_id]
        if ('gre_field' in s.protocol_name) and ('P' in s.image_type):
            info[phase] = [s.series_id]
        if ('t2_tse_dark-fluid_tra_p3' in s.protocol_name) and (s.dim3 == 47):
            info[t2_flair] = [s.series_id]

        # diffusion images and se fmaps
        if ('cmrr_fieldmapse_ap' in s.protocol_name) and (s.dim4 == 2):
            info[dwi_ap] = [s.series_id]
        if ('cmrr_fieldmapse_pa' in s.protocol_name) and (s.dim4 == 2):
            info[dwi_pa] = [s.series_id]
        if ('cmrr_mb3hydi_ipat2_64ch' in s.protocol_name) and (s.dim4 == 145):
            info[dwi] = [s.series_id]


        # functionals: mag, phase, and sbref
        if (s.dim4 == 1020) and ('Shared' in s.protocol_name) and ('_Pha' not in s.series_description):
            info[sharedreward_mag].append(s.series_id)
        if ('Shared' in s.protocol_name) and ('TR1615_SBRef' in s.series_description) and ('_Pha' not in s.series_description):
            info[sharedreward_sbref].append(s.series_id)
        if (s.dim4 == 1020) and ('Shared' in s.protocol_name) and ('TR1615_Pha' in s.series_description):
            info[sharedreward_phase].append(s.series_id)


        if (s.dim4 == 960) and ('MID' in s.series_description) and ('_Pha' not in s.series_description):
            info[mid_mag].append(s.series_id)
        if ('MID' in s.series_description) and ('TR1615_SBRef' in s.series_description) and ('_Pha' not in s.series_description):
            info[mid_sbref].append(s.series_id)
        if (s.dim4 == 960) and ('MID' in s.series_description) and ('TR1615_Pha' in s.series_description):
            info[mid_phase].append(s.series_id)



    return info

POPULATE_INTENDED_FOR_OPTS = {
                'matching_parameters': ['ModalityAcquisitionLabel'],
                'criterion': 'Closest'
}
