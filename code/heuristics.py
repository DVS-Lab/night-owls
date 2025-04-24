import os

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
        if template is None or not template:
                raise ValueError('Template must be a valid format string')
        return template, outtype, annotation_classes

def infotodict(seqinfo):
    t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')
    nm = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-NM_T2star')
    sharedreward_mag = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-sharedreward_run-{item:d}_part-mag_bold')
    sharedreward_phase = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-sharedreward_run-{item:d}_part-phase_bold')
    sharedreward_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-sharedreward_run-{item:d}_sbref')
    mid_mag = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-mid_run-{item:d}_part-mag_bold')
    mid_phase = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-mid_run-{item:d}_part-phase_bold')
    mid_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-mid_run-{item:d}_sbref')
    rest_mag = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-{item:d}_part-mag_bold')
    rest_phase = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-{item:d}_part-phase_bold')
    rest_sbref = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-{item:d}_sbref')
 

    info = {t1w: [], nm: [], 
            sharedreward_mag: [], sharedreward_phase: [], sharedreward_sbref: [],
            rest_mag: [], rest_phase: [], rest_sbref: [],
            mid_mag: [], mid_phase: [], mid_sbref: []}
    
    list_of_ids = [s.series_id for s in seqinfo]

    for s in seqinfo:

        # anatomicals and neuromelanin
        if ('T1w-anat_mpg_07sag_iso' in s.protocol_name) and ('NORM' in s.image_type):
            info[t1w] = [s.series_id]
        if ('neuromelanin' in s.protocol_name):
            info[nm] = [s.series_id]
        
        # functionals
        if (s.dim4 > 1000) and ('Shared' in s.protocol_name) and ('NORM' in s.image_type):
            info[sharedreward_mag].append(s.series_id)
            idx = list_of_ids.index(s.series_id)
            info[sharedreward_sbref].append(list_of_ids[idx -1])
        if (s.dim4 > 1000) and ('Shared' in s.protocol_name) and ('NORM' not in s.image_type):
            info[sharedreward_phase].append(s.series_id)

        if (s.dim4 > 1000) and ('MID' in s.protocol_name) and ('NORM' in s.image_type):
            info[mid_mag].append(s.series_id)
            idx = list_of_ids.index(s.series_id)
            info[mid_sbref].append(list_of_ids[idx -1])
        if (s.dim4 > 1000) and ('MID' in s.protocol_name) and ('NORM' not in s.image_type):
            info[mid_phase].append(s.series_id)

        if (s.dim4 > 1200) and ('resting-state' in s.protocol_name) and ('NORM' in s.image_type):
            info[rest_mag].append(s.series_id)
            idx = list_of_ids.index(s.series_id)
            info[rest_sbref].append(list_of_ids[idx -1])
        if (s.dim4 > 1200) and ('resting-state' in s.protocol_name) and ('NORM' not in s.image_type):
            info[rest_phase].append(s.series_id)



    return info

POPULATE_INTENDED_FOR_OPTS = {
                'matching_parameters': ['ModalityAcquisitionLabel'],
                'criterion': 'Closest'
}
