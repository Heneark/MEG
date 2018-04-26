# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:30:25 2018

@author: oussama.abdoun
"""

#IMPORT PACKAGES, DEFINE ENVIRONMENT VARIABLES, CREATE PATHES
#==============================================================================
import matplotlib
matplotlib.use('Agg') #not to display figures and run all subjects at once
from header import *
#==============================================================================


#%% PARAMETERS
#==============================================================================
t0 = time.perf_counter()
task = 'SMEG' #'MIMOSA'
states = ['RS','FA','OM']
subjects = get_subjlist(task)
subjects = subjects[:subjects.index('109')]
reject = ['072']#072, 109+: no MRI
for sub in reject:
    if sub in subjects:
        subjects.remove(sub)

# # Last subject preprocessed: 109, verbose='WARNING'
# # Future subjects list:
#subjects = subjects[subjects.index('109')+1:]
subjects.sort()
#==============================================================================


#%% ANATOMICAL RECONSTRUCTION
# FREESURFER
#==============================================================================
# recon-all -i <sub_T1.nii> -s <sub> -all
# # INPUT: MRI T1 raw data
# # OUPUT in FreeSurfer SUBJECTS_DIR
#==============================================================================


# # Run anatomy functions in a terminal (does not work with IPython).
# # To process all subjects in a loop, uncomment "import matplotlib; matplotlib.use('Agg')" at the top of this script
from anat import BEM, src_space

#for sub in subjects:
#    if op.isdir(op.join(os.environ['SUBJECTS_DIR'], sub)) and not op.isdir(op.join(os.environ['SUBJECTS_DIR'], sub, 'bem')):
#        watershed = not op.isfile(op.join(os.environ['SUBJECTS_DIR'], sub, 'bem', 'brain.surf'))
#        BEM(subject=sub, watershed=watershed)
#        src_space(subject=sub)


#subjectlist = ''
#for sub in subjects:
#    if op.isdir(op.join(os.environ['SUBJECTS_DIR'], sub, 'bem')) and not op.isfile(op.join(os.environ['SUBJECTS_DIR'], sub, 'bem', sub+'-head-dense.fif')):
#        subjectlist += sub+'\n'
#with open('bash_subject_list.txt', 'w') as fid:
#    fid.write(subjectlist)

# TERMINAL COMMAND
#==============================================================================
# while read p; do mne make_scalp_surfaces -s $p -o; done < bash_subject_list.txt
# # creates SUBJECTS_DIR/<subject>/bem/<subject>-head-dense.fif
#==============================================================================


#%% HEAD POSITION: MATLAB
#==============================================================================
# Adjust_Head_Pos_3DirbyCoil.m
# # Select which block to coregister --> Analyses/<task>/<meg>/Coregistration/<subject>/<precision>_<blocki_blockj>.hc
#==============================================================================


#%% PREPROCESSING
from preproc import process, R_T_ECG_events, epoch, empty_room_covariance, check_ecg_epoch

custom_ecg = {'004': {'R_sign': 1, 'heart_rate': 78, 'tstart': {'RS01': .5, 'OM02': .15, 'FA04': .7}, 'force':True},
              '010': {'R_sign': -1, 'heart_rate': 77, 'T_sign': 1},
              '012': {'R_sign': -1, 'heart_rate': 77, 'T_sign': 1},
              '028': {'R_sign': -1, 'heart_rate': 55},
              '069': {'R_sign': -1, 'heart_rate': 94}}

for sub in subjects:
    if not op.isfile(op.join(Analysis_path, task, 'meg', 'Covariance', sub, 'empty_room-cov.fif')):
        empty_room_covariance(task, sub)
    for state in states:
        for blk in get_blocks(sub, task=task, state=state):
            custom_args=dict()
            if sub in custom_ecg.keys():
                custom_args = custom_ecg[sub].copy()
            if 'tstart' in custom_args.keys():
                custom_args['tstart'] = custom_args['tstart'][state+blk]
            
#            data = process(task, sub, state, blk, ica_rejection={'mag':7000e-15}, ECG_threshold=0.2, EOG_threshold=5, check_ica=False, custom_args=custom_args)
#            events, event_id = R_T_ECG_events(task, sub, state, blk, data['ECG_included'], custom_args)
#            check_ecg_epoch(task, sub, state, blk, data['ECG_included'], events=events, save=True)
#            for i, (name, raw) in enumerate(data.items()):
#                epochs = epoch(task, sub, state, blk, raw, name, events, event_id)

#reject = list(custom_ecg.keys())
#for sub in reject:
#    if sub in subjects:
#        subjects.remove(sub)
#for sub in subjects[subjects.index('071'):]:
##    empty_room_covariance(task, sub)
#    for state in states:
#        for blk in get_blocks(sub, task=task, state=state):
##            raw,raw_ECG = process(task='MIMOSA', sub, state, blk, ica_rejection={'mag':7000e-15}, ECG_threshold=0.2, EOG_threshold=5)
#            epochs = epoch(task, sub, state, blk, save_t_timing=True, names=['T_ECG_included','T_ECG_excluded'])


#%% COREGISTRATION (https://www.slideshare.net/mne-python/mnepython-coregistration)
#==============================================================================
#%gui wx
#mne.gui.coregistration()
# # on "Save", creates SUBJECTS_DIR/<subject>/bem/<subject>-fiducials.fif
# # Files corresponding to the coregistration: Analyses/<task>/meg/Coregistration/<subject>/<HCfilename>'-trans.fif'
#==============================================================================


#%% SOURCE RECONSTRUCTION
from source_reconstruction import baseline_covariance, src_rec, fs_average

names = ['R_ECG_included', 'R_ECG_excluded']#, 'T_ECG_included', 'T_ECG_excluded']#
precision = '0.5cm'

for sub in subjects:#['004', '010', '054', '071']:#
    coreg_list = glob.glob(op.join(Analysis_path, task, 'meg', 'Coregistration', sub, '*'+precision+'*-trans.fif'))
    for c,coreg in enumerate(coreg_list):
        coreg_list[c] = set(op.split(coreg)[-1].split(precision)[-1].strip('-trans.fif').split('_')[1:])
    
    for state in states:
        blk_list = set(get_blocks(sub, task=task, state=state))
        coreg_by_state = []
        for coreg in coreg_list:
            if blk_list & coreg:
                coreg_by_state.append(sorted(list(blk_list & coreg)))
        
        for group in coreg_by_state:
#            noise_cov,evoked = baseline_covariance(task, sub, state, block_group=group, rejection={'mag':3500e-15}, baseline=(-.4,-.25), names=names)
            stc_surf,stc_vol = src_rec(task, sub, state, evoked=None, noise_cov=None, block_group=group, names=names, compute_fwd=False, compute_inv=False, compute_stc=True, tmin=.35, tmax=.5)
#
for name in names:
    for state in states:
        fs_average(task, state, name=name, subjects=subjects, do_morphing=False, overwrite=False)


#%%
t1 = time.perf_counter()
T = t1 - t0
print(colored(time.strftime('Finished %c',time.localtime()),'blue'))
print(colored('Elapsed time: {d}d {h}h {m}min {s}s'.format(s=round(T%60), m=round((T - T%60)%(60*60)/60), h=round((T - T%(60*60))%(24*60*60)/(60*60)), d=round((T - T%(24*60*60))/(24*60*60))), 'green'))