#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 14:32:51 2018

@author: benjamin.ador
"""

#IMPORT PACKAGES, DEFINE ENVIRONMENT VARIABLES, CREATE PATHES
#==============================================================================
from header import *
import visbrain
import xarray as xr
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


#%% STATS
names = ['R_ECG_included', 'R_ECG_excluded', 'T_ECG_included', 'T_ECG_excluded']
precision = '0.5cm'
noise_cov = 'baseline_cov'
fsaverage = '-fsaverage'#o
stc_ext = '-lh.stc'

stc_path = op.join(Analysis_path, task, 'meg', 'SourceEstimate')

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
#        
#        for group in coreg_by_state:
##            noise_cov,evoked = baseline_covariance(task, sub, state, block_group=group, baseline=(-.4,-.25), names=names)
#            stc_surf,stc_vol = src_rec(task, sub, state, block_group=group, names=names)


#%%
t1 = time.perf_counter()
T = t1 - t0
print(colored(time.strftime('Finished %c',time.localtime()),'blue'))
print(colored('Elapsed time: {d}d {h}h {m}min {s}s'.format(s=round(T%60), m=round((T - T%60)%(60*60)/60), h=round((T - T%(60*60))%(24*60*60)/(60*60)), d=round((T - T%(24*60*60))/(24*60*60))), 'green'))