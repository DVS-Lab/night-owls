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
	trust = create_key('sub-{subject}/func/sub-{subject}_task-trust_run-{item:d}_bold')
	trust_sbref = create_key('sub-{subject}/func/sub-{subject}_task-trust_run-{item:d}_sbref')
	sharedreward = create_key('sub-{subject}/func/sub-{subject}_task-sharedreward_run-{item:d}_bold')
	sharedreward_sbref = create_key('sub-{subject}/func/sub-{subject}_task-sharedreward_run-{item:d}_sbref')
	srSocial = create_key('sub-{subject}/func/sub-{subject}_task-socialdoors_run-{item:d}_bold')
	srSocial_sbref = create_key('sub-{subject}/func/sub-{subject}_task-socialdoors_run-{item:d}_sbref')
	srDoors = create_key('sub-{subject}/func/sub-{subject}_task-doors_run-{item:d}_bold')
	srDoors_sbref = create_key('sub-{subject}/func/sub-{subject}_task-doors_run-{item:d}_sbref')
	UGR = create_key('sub-{subject}/func/sub-{subject}_task-ugr_run-{item:d}_bold')
	UGR_sbref = create_key('sub-{subject}/func/sub-{subject}_task-ugr_run-{item:d}_sbref')
	dwi = create_key('sub-{subject}/dwi/sub-{subject}_dwi')
	dwi_pa = create_key('sub-{subject}/fmap/sub-{subject}_acq-dwi_dir-PA_epi')
	dwi_ap = create_key('sub-{subject}/fmap/sub-{subject}_acq-dwi_dir-AP_epi')


	info = {t1w: [],
	      mag: [],
	      phase: [],
	      dwi: [],
	      dwi_pa: [],
	      dwi_ap: [],
	      t2_flair: [],
	      trust: [],
	      trust_sbref: [],
	      sharedreward: [],
	      sharedreward_sbref: [],
	      srSocial: [],
	      srSocial_sbref: [],
	      srDoors: [],
	      srDoors_sbref: [],
	      UGR: [],
	      UGR_sbref: [],
	      }
	
	
	list_of_ids = [s.series_id for s in seqinfo]
	
	for s in seqinfo:
	     if ('T1w-anat_mpg_07sag_iso' in s.protocol_name) and ('NORM' in s.image_type):
	         info[t1w] = [s.series_id]
	     #if ('gre_field' in s.protocol_name) and ('NORM' in s.image_type):
	     #    info[mag] = [s.series_id]
	     #if ('gre_field' in s.protocol_name) and ('P' in s.image_type):
	     #    info[phase] = [s.series_id]
	     if ('t2_tse_dark-fluid_tra_p3' in s.protocol_name) and (s.dim3 == 47):
	     		info[t2_flair] = [s.series_id]
	 
	     if ('cmrr_fieldmapse_ap' in s.protocol_name) and (s.dim4 == 2):
	     		info[dwi_ap] = [s.series_id]
	     if ('cmrr_fieldmapse_pa' in s.protocol_name) and (s.dim4 == 2):
	     		info[dwi_pa] = [s.series_id]
	     if ('cmrr_mb3hydi_ipat2_64ch' in s.protocol_name) and (s.dim4 == 145):
	     		info[dwi] = [s.series_id]
	 
	     if (s.dim4 == 1120) and ('Trust' in s.protocol_name) and ('NORM' in s.image_type):
	         info[trust].append(s.series_id)
	         idx = list_of_ids.index(s.series_id)
	         info[trust_sbref].append(list_of_ids[idx -1])
	
	     if (s.dim4 == 1020) and ('Shared' in s.protocol_name) and ('NORM' in s.image_type):
	         info[sharedreward].append(s.series_id)
	         idx = list_of_ids.index(s.series_id)
	         info[sharedreward_sbref].append(list_of_ids[idx -1])
	
	     if (s.dim4 == 872) and ('SocialDoors_face' in s.protocol_name) and ('NORM' in s.image_type):
	         info[srSocial] = [s.series_id]
	         idx = list_of_ids.index(s.series_id)
	         info[srSocial_sbref].append(list_of_ids[idx -1])
	
	     if (s.dim4 == 872) and ('SocialDoors_doors' in s.protocol_name) and ('NORM' in s.image_type):
	         info[srDoors] = [s.series_id]
	         idx = list_of_ids.index(s.series_id)
	         info[srDoors_sbref].append(list_of_ids[idx -1])
	
	     if (s.dim4 == 960) and ('UGR' in s.protocol_name) and ('NORM' in s.image_type):
	         info[UGR].append(s.series_id)
	         idx = list_of_ids.index(s.series_id)
	         info[UGR_sbref].append(list_of_ids[idx -1])
	
	return info

POPULATE_INTENDED_FOR_OPTS = {
        'matching_parameters': ['ModalityAcquisitionLabel'],
        'criterion': 'Closest'
}

