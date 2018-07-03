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
# # Last subject preprocessed: 101
# # Future subjects list:
#subjects = subjects[subjects.index('101')+1:]

reject = ['069', '072', '074', '079', '098', '109']#074, 079, 098, 109: no MRI ; '069': NOISE
for sub in reject:
    if sub in subjects:
        subjects.remove(sub)

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
# while read p; do mne make_scalp_surfaces -s $p -o --force; done < bash_subject_list.txt
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
              '057': {'R_sign': -1},
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
#            try:
#            step='process'
#            raw = process(task, sub, state, blk, ica_rejection={'mag':7000e-15}, ECG_threshold=0.2, EOG_threshold=5, check_ica=True, overwrite_ica=True, custom_args=custom_args)
#            step='epoch'
#            events, event_id = R_T_ECG_events(task, sub, state, blk, raw, custom_args)
#            step='ECG check'
#            check_ecg_epoch(task, sub, state, blk, raw, events, save=True)
#            except:
#                with open('run.log', 'a') as fid:
#                    fid.write(sub+'\t'+state+'\t'+blk+'\t'+'preproc bug\tstep\n')
#                pass


#%% COREGISTRATION (https://www.slideshare.net/mne-python/mnepython-coregistration)
#==============================================================================
#import os, mne
#os.environ["ETS_TOOLKIT"] = "wx"
#%gui wx
#mne.gui.coregistration()
# # on "Save", creates SUBJECTS_DIR/<subject>/bem/<subject>-fiducials.fif
# # Files corresponding to the coregistration: Analyses/<task>/meg/Coregistration/<subject>/<HCfilename>'-trans.fif'
#==============================================================================


#%% SOURCE RECONSTRUCTION
from source_reconstruction import ERP, src_rec, fs_average, src_rec15

names = ['ECG_included', 'ECG_excluded']
keys = ['R', 'T']
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
        
#        for group in coreg_by_state:
#            for name in names:
##                try:
#                noise_cov, evoked = ERP(task, sub, state, block_group=group, name=name, cov_keys=['R'], rejection={'mag':3500e-15}, baseline={'R':(None,None), 'T':(None,None)})#, tmin=-.8, tmax=.8)
#                src_rec15(task, sub, state, block_group=group, evoked=evoked, name=name, window={'R':(.35,.425), 'T':(.1,.25)}, volume=0, baseline_cov=False, method='beamformer')
#                stc_surf, stc_vol = src_rec(task, sub, state, block_group=group, keys=keys, name=name, compute_fwd=False, baseline_cov=False, window={'T':(-.175,-.075)})
#                stc_surf, stc_vol = src_rec(task, sub, state, block_group=group, keys=keys, name=name, compute_fwd=False, baseline_cov=False, window={'T':(.1,.2)})
#                except:
#                    with open('run.log', 'a') as fid:
#                        fid.write(sub+'\t'+state+'\t'+str(group)+'\t'+name+'\t'+'source reconstruction bug\n')
#                    pass

#windows = window={'R':(.35,.425), 'T':(.1,.25)}
#for name in names:
#    for key in keys:
#        for state in states:
##            for window in windows:
#                fs_average(task, state, name=name, key=key, subjects=subjects, do_morphing=True, overwrite=True, baseline_cov=False, window=windows[key])


#%% FULL EVOKED DATASET
verb = mne.set_log_level(False, return_old_level=True)

evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked')
names = ['ECG_included', 'ECG_excluded']
keys = ['R', 'T']
sfreq = 200
data_file = op.join(evoked_path, 'DATASET-{}Hz.nc'.format(sfreq))
#for name in names:
#    for k in keys:
#        kname = k+'_'+name
#        for st,state in enumerate(states):
#            for su,sub in enumerate(tqdm(subjects)):
#                evo_list = []
#                for b,blk in enumerate(get_blocks(sub, task=task, state=state)):
#                    evo = mne.Evoked(op.join(evoked_path, sub, '{}_{}_{}-{}-ave.fif'.format(sub, state, blk, name)), condition = k)
#                    if sfreq:
#                        evo.resample(sfreq)
#                    if not st and not su and not b:
#                        times = evo.times
#                        channels = [ch.split('-')[0] for c,ch in enumerate(evo.ch_names) if c in mne.pick_types(evo.info, ref_meg=False)]
#                        evoked = xr.DataArray(np.zeros((len(states), len(subjects), times.size, len(channels))), dims=['state', 'subject', 'time', 'sensor'], coords={'state':states, 'subject':subjects, 'time':times, 'sensor':channels})
#                    
#                    evo.info['bads'] = list(set(evo.info['bads']) | set(get_chan_name(sub, 'bad', evo)))
#                    evo.interpolate_bads(verbose='ERROR')
#                    evo_list.append(evo)
#                
#                evo = mne.combine_evoked(evo_list, 'nave')
#                evoked[st,su] = evo.data[mne.pick_types(evo.info, ref_meg=False)].T
#        
#        if not op.isfile(data_file):
#            mode = 'w'
#        else:
#            mode = 'a'
#        evoked.to_netcdf(path=data_file, group=kname, mode=mode)
#        
#        del evoked

mne.set_log_level(verb)


#%% FULL DATASET
verb = mne.set_log_level(False, return_old_level=True)

stc_path = op.join(Analysis_path, task, 'meg', 'SourceEstimate')
#names = ['T_ECG_included', 'T_ECG_excluded']#'R_ECG_included', 'R_ECG_excluded', 
noise_cov = 'empty_room_cov'
stc_ext = '-lh.stc'
sfreq = 200
win = (-.175,-.075)#(.1,.2)

surface = 'ico4'
win = {'R':(.35,.425), 'T':(.1,.25)}
data_file = op.join(stc_path, 'fsaverage', 'DATASET-{}-surface_{}-{}Hz.nc'.format(noise_cov, surface, sfreq))#, *win[]
for name in names:
    for k in keys:
        kname = k+'_'+name
        for st,state in enumerate(states):
            print(state, kname)
            for su,sub in enumerate(tqdm(subjects)):
                stc_file = op.join(stc_path, 'fsaverage', sub, '{}_{}-{}-{}*{}_{}*{}'.format(sub, state, kname, noise_cov, *win[k], stc_ext))
                file = glob.glob(stc_file)[0]
                data = mne.read_source_estimate(file[:file.index(stc_ext)])
                if sfreq:
                    data.resample(sfreq)
                if not st and not su:#and not n 
                    times = data.times
                    vertices = np.concatenate([['lh_' + str(x) for x in data.lh_vertno],['rh_' + str(x) for x in data.rh_vertno]])
                    stc = np.zeros((len(states), len(subjects), *data.data.T.shape))
                    stc = xr.DataArray(stc, dims=['state', 'subject', 'time', 'src'], coords={'state':states, 'subject':subjects, 'time':times, 'src':vertices})
                stc[st,su] = data.data.T
                del data
        if not op.isfile(data_file):
            mode = 'w'
        else:
            mode = 'a'
        stc.to_netcdf(path=data_file, group=kname, mode=mode)
        del stc

mne.set_log_level(verb)


#%%
t1 = time.perf_counter()
T = t1 - t0
print(colored(time.strftime('Finished %c',time.localtime()),'blue'))
print(colored('Elapsed time: {d}d {h}h {m}min {s}s'.format(s=round(T%60), m=round((T - T%60)%(60*60)/60), h=round((T - T%(60*60))%(24*60*60)/(60*60)), d=round((T - T%(24*60*60))/(24*60*60))), 'green'))