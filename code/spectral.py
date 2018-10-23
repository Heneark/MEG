#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:10:04 2018

@author: benjamin.ador
"""
import time
t0 = time.perf_counter()

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
                if psd_func is psd_welch:
                    psds, freqs = psd_func(raw, n_fft=int(10*raw.info['sfreq']), fmax=int(raw.info['sfreq']/4), n_jobs=4)
                if psd_func is psd_multitaper:
                    psds, freqs = psd_func(raw, fmax=int(raw.info['sfreq']/4), n_jobs=4)
                if not st and not su and not b:
                    channels = [ch.split('-')[0] for c,ch in enumerate(raw.ch_names) if c in mne.pick_types(raw.info, ref_meg=False)]
                    PSD = xr.DataArray(np.zeros((len(state_blk), len(subjects), len(channels), freqs.size)), 
                                       dims=['state', 'subject', 'chan', 'freq'], 
                                       coords={'state':state_blk, 'subject':subjects, 'chan':channels, 'freq':freqs})
                PSD.loc[state+str(b+1), sub] = psds
    
    PSD.to_netcdf(path=op.join(Analysis_path, task, 'meg', 'Raw', 'PSD_{}.nc'.format(psd_func.__name__)))
    mne.set_log_level(verb)
    return PSD


#%% COMPUTE OR LOAD PSD

task = 'SMEG'
subjects = get_subjlist(task)
subjects.sort()

#warnings.filterwarnings("ignore",category=DeprecationWarning)
#PSD = psd_average(task, subjects, psd_func=psd_multitaper)

PSD = xr.open_dataarray(op.join(Analysis_path, task, 'meg', 'Raw', 'PSD.nc'))
PSD.load()

PSD_norm = PSD/PSD.mean('freq')


#%%

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

channels = {'all': PSD.chan.values}
channels['left'] = filter(channels['all'], 'ML*')
channels['left frontal'] = filter(channels['all'], 'MLF*')
channels['right'] = filter(channels['all'], 'MR*')
channels['central'] = filter(channels['all'], 'MZ[!O]*') #midline: all but occipital
channels['frontal'] = filter(channels['all'], 'M*F*')
channels['occipital'] = filter(channels['all'], 'M*O*')


#%% SELECTION PARAMETERS

state = 'FA1'

sub_key = 'all'
sub = '004'

chan_key = 'all' #select a subset of channels

fmin = 1 #PSD.freq.values[0]
fmax = PSD.freq.values[-1]


#%% SINGLE SUBJECT

data = PSD.loc[state, sub, channels[chan_key], fmin:fmax]
average = data.mean('chan')
sem = data.std('chan')/sqrt(data.chan.size)
plt.figure()
plt.semilogx(data.freq, average)
plt.fill_between(data.freq, average+sem, average-sem, alpha=.4)
plt.title('Average PSD of subject {sub} for state {s}\nAverage of {ch} channels'.format(s=state, sub=sub, ch=chan_key))


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


#%% PLOT BY GROUP

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


#%% PLOT BY STATE

#plt.figure()
#keys = ['novice', 'expert']
#sub_key = 'all'
#for stt in PSD.state.values:
##for sub_key in keys:
#    data = PSD.loc[stt, subjects[sub_key], channels[chan_key], fmin:fmax].mean('chan')
#    average = data.mean('subject')
#    sem = data.std('subject')/sqrt(data.subject.size)
#    plt.semilogx(data.freq, average)
#    plt.fill_between(data.freq, average+sem, average-sem, alpha=.4)
#plt.title('Average PSD over {sub} subjects\nAverage of {ch} channels'.format(s=state, sub=sub_key, ch=chan_key))
#plt.legend(PSD.state.values)
#plt.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', '{}-{}.png'.format(state, '+'.join(keys))))


#%% ALPHA PEAK

fmin = 8
fmax = 12

lateral = ['*', 'L', 'R', 'Z']
ante_post = ['*', 'C', 'F', 'P', 'O', 'T']

alpha_tsv = 'group\tsubject\tstate\tlateral\tante_post\tpeak_freq\tpeak_val\tpeak_norm\n'
alpha_fname = op.join(Analysis_path, task, 'meg', 'Raw', 'alpha_peak_{}_{}.tsv'.format(fmin, fmax))

for sub in tqdm(PSD.subject.values):
    for stt in PSD.state.values:
        if sub in no_blk2 and fnmatch.fnmatch(stt, '*2'):
            continue #Old subjects only did one session of meditation states
        
        for lat in lateral:
            for a_p in ante_post:
                chs = filter(PSD.chan.values.tolist(), 'M'+lat+a_p+'*')
                if lat is 'Z':
                    chs = filter(PSD.chan.values.tolist(), 'M'+lat+'[!O]*') #Exclude Occipital channels from midline
                
                psd_sel = PSD.loc[stt, sub, chs, fmin:fmax].mean('chan')
                f_peak = psd_sel.freq[psd_sel.argmax()].values
                peak_val = psd_sel[psd_sel.argmax()].values
                peak_norm = PSD_norm.loc[stt, sub, chs, fmin:fmax].mean('chan')[psd_sel.argmax()].values
                
                alpha_tsv += '{}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\n'.format(expertise(sub), sub, stt, lat, a_p, f_peak, peak_val, peak_norm)
                
                if lat is 'Z':
                    break #Only one set of channels for midline

with open(alpha_fname, 'w') as fid:
    fid.write(alpha_tsv)
    fid.close()


#%% RUNNING TIME

t1 = time.perf_counter()
T = t1 - t0
print(colored(time.strftime('Finished %c',time.localtime()),'blue'))
print(colored('Elapsed time: {d}d {h}h {m}min {s}s'.format(s=round(T%60), m=round((T - T%60)%(60*60)/60), h=round((T - T%(60*60))%(24*60*60)/(60*60)), d=round((T - T%(24*60*60))/(24*60*60))), 'green'))