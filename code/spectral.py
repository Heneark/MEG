#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:10:04 2018

@author: benjamin.ador
"""

from header import *
from mne.time_frequency import psd_welch, psd_multitaper
from typing import List, Union


#%%

def psd_average(task: str, subjects: List[str]=None, states: List[str]=['RS', 'FA', 'OM'], n_blk={'RS': 1, 'FA': 2, 'OM': 2}, psd_func: Union[psd_welch, psd_multitaper]=psd_welch, verbose=False):
    """
    
    """
    verb = mne.set_log_level(verbose, return_old_level=True)
    if not subjects:
        subjects = sorted(get_subjlist(task))
    
    state_blk = []
    for state in states:
        state_blk.extend([state + str(b+1) for b in range(n_blk[state])])
    
    for st,state in enumerate(states):
        print(state)
        for su,sub in enumerate(tqdm(subjects)):
            blocks = get_blocks(sub, task=task, state=state)
            for b,blk in enumerate(blocks):
                raw = load_preproc(task, sub, state, blk)
                psds, freqs = psd_func(raw, n_fft=int(10*raw.info['sfreq']), n_jobs=4)
                if not st and not su and not b:
                    channels = [ch.split('-')[0] for c,ch in enumerate(raw.ch_names) if c in mne.pick_types(raw.info, ref_meg=False)]
                    PSD = xr.DataArray(np.zeros((len(state_blk), len(subjects), len(channels), freqs.size)), 
                                       dims=['state', 'subject', 'chan', 'freq'], 
                                       coords={'state':state_blk, 'subject':subjects, 'chan':channels, 'freq':freqs})
                PSD.loc[state+str(b+1), sub] = psds
    
    PSD.to_netcdf(path=op.join(Analysis_path, task, 'meg', 'Raw', 'PSD.nc'))
    mne.set_log_level(verb)
    return PSD


#%%

task = 'SMEG'
subjects = get_subjlist(task) + ['053']
subjects.sort()

warnings.filterwarnings("ignore",category=DeprecationWarning)
PSD = psd_average(task, subjects)
