#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:10:04 2018

@author: benjamin.ador
"""

from header import *

from fnmatch import filter
from mne.time_frequency import psd_welch, psd_multitaper
from typing import List, Union


#%% PSD TO DATAARRAY FUNCTION

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


#%% COMPUTE OR LOAD PSD

task = 'SMEG'
subjects = get_subjlist(task)
subjects.sort()

#warnings.filterwarnings("ignore",category=DeprecationWarning)
#PSD = psd_average(task, subjects)

PSD = xr.open_dataarray(op.join(Analysis_path, task, 'meg', 'Raw', 'PSD.nc'))
PSD.load()


#%% SELECTION PARAMETERS

state = 'RS1'

subjects = {'all': sorted(get_subjlist(task))}
no_blk2 = ['002', '004', '007', '016'] #FA2 and OM2 don't exist for these subjects (PSD=0)

subjects['expert'] = list(); subjects['novice'] = list()
for sub in subjects['all']:
    if expertise(sub) is 'N':
        subjects['novice'].append(sub)
    elif expertise(sub) is 'E':
        subjects['expert'].append(sub)
    else:
        warnings.warn('Expertise of subject {} is not specified, check masterfile.'.format(sub))
sub_key = 'all'

channels = {'all': PSD.chan.values}
channels['left'] = filter(channels['all'], 'ML*')
channels['right'] = filter(channels['all'], 'MR*')
channels['central'] = filter(channels['all'], 'MZ*')
chan_key = 'all' #select a subset of channels

fmin = .5 #PSD.freq.values[0]
fmax = PSD.freq.values[-1]


#%% SELECT DATA

data = PSD.loc[state, subjects[sub_key], channels[chan_key], fmin:fmax].mean('chan')
average = data.mean('subject')
sem = data.std('subject')/sqrt(data.subject.size)


#%% PLOT

plt.figure()
plt.semilogx(data.freq, average)
plt.fill_between(data.freq, average+sem, average-sem, alpha=.4)
plt.title('Average PSD over {sub} subjects for state {s}\nAverage of {ch} channels'.format(s=state, sub=sub_key, ch=chan_key))
plt.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', '{}-{}_subs.png'.format(state, sub_key)))


#%%

plt.figure()
keys = ['novice', 'expert']
for sub_key in keys:
    data = PSD.loc[state, subjects[sub_key], channels[chan_key], fmin:fmax].mean('chan')
    average = data.mean('subject')
    sem = data.std('subject')/sqrt(data.subject.size)
    plt.semilogx(data.freq, average)
    plt.fill_between(data.freq, average+sem, average-sem, alpha=.4)
plt.title('Average PSD for state {s}\nAverage of {ch} channels'.format(s=state, sub=sub_key, ch=chan_key))
plt.legend(keys)
plt.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', '{}-{}.png'.format(state, '+'.join(keys))))
