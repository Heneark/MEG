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

reject = ['004', '010']
for sub in reject:
    if sub in subjects:
        subjects.remove(sub)

# # Last subject preprocessed: 109
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
# # Select which block to coregister --> Analyses/<task>/<meg>/HC_for_coreg/<subject>/<precision>_<blocki_blockj>.hc
#==============================================================================


#%% PREPROCESSING
from preproc import process, epoch, empty_room_covariance

#for sub in subjects:
#    empty_room_covariance(task, sub)
#    for state in states:
#        for blk in get_blocks(sub, task=task, state=state):
#            raw,raw_ECG = process(task='MIMOSA', sub, state, blk, ica_rejection={'mag':7000e-15}, ECG_threshold=0.2, EOG_threshold=5)
#            epochs = epoch(task, sub, state, blk, high_pass=.5, low_pass=None)


#%% COREGISTRATION (https://www.slideshare.net/mne-python/mnepython-coregistration)
#==============================================================================
#%gui wx
#mne.gui.coregistration()
# # on "Save", creates SUBJECTS_DIR/<subject>/bem/<subject>-fiducials.fif
# # Files corresponding to the coregistration: Analyses/<task>/meg/Coregistration/<subject>/<HCfilename>'-trans.fif'
#==============================================================================


#%% SOURCE RECONSTRUCTION
from source_reconstruction import baseline_covariance, src_rec

names = ['R_ECG_included', 'R_ECG_excluded', 'T_ECG_included', 'T_ECG_excluded']
precision = '0.5cm'

for sub in subjects:
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
#            noise_cov,evoked = baseline_covariance(task, sub, state, block_group=group, baseline=(-.4,-.25), names=names)
            stc_surf,stc_vol = src_rec(task, sub, state, block_group=group, names=names)


#%%
t1 = time.perf_counter()
T = t1 - t0
print(colored(time.strftime('Finished %c',time.localtime()),'blue'))
print(colored('Elapsed time: {d}d {h}h {m}min {s}s'.format(s=round(T%60), m=round((T - T%60)%(60*60)/60), h=round((T - T%(60*60))%(24*60*60)/(60*60)), d=round((T - T%(24*60*60))/(24*60*60))), 'green'))