# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 19:00:23 2017

@author: oussama.abdoun
"""

import mne
import numpy as np
import pandas as pd
import os, sys
import locale
import matplotlib.pyplot as plt

sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from csv_io_alex import *


locale.setlocale(locale.LC_ALL, "en_US.utf8")

subj = '040'

pathBase = "/home/oussama.abdoun/Documents/donneesMEG/"

ICA_path = '/home/oussama.abdoun/Documents/donneesMEG/' + get_id(subj) + '/ICA/'
#==============================================================================
# POUR SUJET 037
pathBase = '/dycog/meditation/ERC/Raw data/MIMOSA/meg/'
ICA_path = '/dycog/meditation/ERC/Analyses/MIMOSA/meg/ICA/' + get_id(subj)
#==============================================================================


n_components = 0.975






def cleanEpochsEventID(epochs):
    tmp = [epochs.events[k,2] for k in xrange(epochs.events.shape[0])]
    tmp = np.unique(np.array(tmp))
    event_id = dict()
    for code in tmp:
        event_id[str(code)] = code
    epochs.event_id = event_id
    
    return epochs




listblk = ['{:02d}'.format(b) for b in xrange(1,9)]
#==============================================================================
# EN UTILISANT LE CSV
# listblk = []
# for key in get_blocks(subj):
#     listblk.append(key)
#==============================================================================


listevoked_std = []
listevoked_std_phase0 = []
listevoked_std_phase45 = []
listevoked_std_blue = []
listevoked_std_red = []
listevoked_dev = []

listevoked_targ_on = []
listevoked_targ_off = []


if subj == "028":
    listblkFA = ['03','04','05','11','12','13']
    listblkOM = ['07','08','09','15','16','17']
    
if subj == "030":
    listblkFA = ['03','04','05','15','16','17']
    listblkOM = ['07','08','09','11','12','13']

if subj == "032":
    listblkFA = ['03','04','05','15','16','17']
    listblkOM = ['08','09','11','12','13']

if subj == "037":
    listblkFA = ['03','04','05','11','12','13']
    listblkOM = ['07','08','09','15','16','17']


for blk in listblkFA:
    filename = os.path.join(pathBase, get_rawpath(subj)[0], get_rawpath(subj)[1] + blk + '.ds')
        
    raw = mne.io.read_raw_ctf(filename,preload=True)
    
    layout = mne.channels.find_layout(raw.info, ch_type='mag')
    picks_meg = mne.pick_types(raw.info, meg=True)
    
    ica = mne.preprocessing.read_ica(os.path.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components) + '_blk' + blk + '-ica.fif'))
    
    raw.notch_filter(np.arange(50,251,50))
    raw.filter(l_freq=0.1, h_freq=40, n_jobs=4)
    
    #raw.plot_psd(fmin=0.1,fmax=40)
    
    events = mne.find_events(raw,stim_channel='UPPT001')
    
    picks_meg = mne.pick_types(raw.info, meg=True)
    
#==============================================================================
#   GABORS
#==============================================================================
    # Extract on-beat gabors
    epochs_targ_on = mne.Epochs(raw, events, event_id=[202,204,206], tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
    evoked_targ_on = epochs_targ_on.average(picks=picks_meg)
    evoked_targ_on = ica.apply(evoked_targ_on)    
    
    # Extract off-beat gabors
    epochs_targ_off = mne.Epochs(raw, events, event_id=[201,203,205,207], tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
    evoked_targ_off = epochs_targ_off.average(picks=picks_meg)
    evoked_targ_off = ica.apply(evoked_targ_off)
    
    listevoked_targ_on.append(evoked_targ_on)
    listevoked_targ_off.append(evoked_targ_off)


#==============================================================================
#   RINGS
#==============================================================================
    # Extract standard epochs    
    epochs_std = mne.Epochs(raw, events, event_id=range(70)+range(100,170), tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
#    epochs_std.plot()
    epochs_std = cleanEpochsEventID(epochs_std)
    evoked_std = epochs_std.average(picks=picks_meg)
    evoked_std = ica.apply(evoked_std)
    listevoked_std.append(evoked_std)
    
    # Standards - phase 0
    tmp = range(70)+range(100,170)
    event_id = [k for k in tmp if int(str(k)[-1])<5]
    epochs_std_phase0 = mne.Epochs(raw, events, event_id=event_id, tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
    epochs_std_phase0 = cleanEpochsEventID(epochs_std_phase0)
    evoked_std_phase0 = epochs_std_phase0.average(picks=picks_meg)
    evoked_std_phase0 = ica.apply(evoked_std_phase0)
    listevoked_std_phase0.append(evoked_std_phase0)
    
    # Standards - phase 45
    event_id = [k for k in tmp if int(str(k)[-1])>=5]
    epochs_std_phase45 = mne.Epochs(raw, events, event_id=event_id, tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
    epochs_std_phase45 = cleanEpochsEventID(epochs_std_phase45)
    evoked_std_phase45 = epochs_std_phase45.average(picks=picks_meg)
    evoked_std_phase45 = ica.apply(evoked_std_phase45)
    listevoked_std_phase45.append(evoked_std_phase45)
    
    # Standards - RED family
    epochs_std_red = mne.Epochs(raw, events, event_id=range(70), tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
    epochs_std_red = cleanEpochsEventID(epochs_std_red)
    evoked_std_red = epochs_std_red.average(picks=picks_meg)
    evoked_std_red = ica.apply(evoked_std_red)
    listevoked_std_red.append(evoked_std_red)
    
    # Standard - BLUE family
    epochs_std_blue = mne.Epochs(raw, events, event_id=range(100,170), tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
    epochs_std_blue = cleanEpochsEventID(epochs_std_blue)
    evoked_std_blue = epochs_std_blue.average(picks=picks_meg)
    evoked_std_blue = ica.apply(evoked_std_blue)
    listevoked_std_blue.append(evoked_std_blue)
  
    # Extract  deviant epochs
    epochs_dev = mne.Epochs(raw, events, event_id=range(90,100)+range(190,200), tmin=-1.0,tmax=1.0, baseline=(-0.1,0), preload=True,proj=False,on_missing='ignore')
#    epochs_dev.plot()  
    epochs_dev = cleanEpochsEventID(epochs_dev)
    evoked_dev = epochs_dev.average(picks=picks_meg)
    evoked_dev = ica.apply(evoked_dev)
    listevoked_dev.append(evoked_dev)
    
    
    
    
    
    
    
    
#==============================================================================
# POOL BLOCKS
#==============================================================================
evoked_targ_on = mne.combine_evoked(listevoked_targ_on, weights='nave')
evoked_targ_off = mne.combine_evoked(listevoked_targ_off, weights='nave')

evoked_std = mne.combine_evoked(listevoked_std, weights='nave')
evoked_std_phase0 = mne.combine_evoked(listevoked_std_phase0, weights='nave')
evoked_std_phase45 = mne.combine_evoked(listevoked_std_phase45, weights='nave')
evoked_std_red = mne.combine_evoked(listevoked_std_red, weights='nave')
evoked_std_blue = mne.combine_evoked(listevoked_std_blue, weights='nave')
evoked_dev = mne.combine_evoked(listevoked_dev, weights='nave')
    
evoked_mmn = mne.combine_evoked([evoked_std,evoked_dev], weights=[-1,1])


#==============================================================================
# PLOT
#==============================================================================
tcrop = (-500,1000)
if 0:
    evoked_std.plot_joint(times=[.105,.190,.270,.480], ts_args=dict(xlim=tcrop,ylim=dict(mag=[-150, 150])), topomap_args=dict(vmin=-100,vmax=100))
    evoked_dev.plot_joint(times=[.105,.190,.270,.445], ts_args=dict(xlim=tcrop,ylim=dict(mag=[-150, 150])), topomap_args=dict(vmin=-100,vmax=100))
    evoked_mmn.plot_joint(times=[.125,.170,.245,.460], ts_args=dict(xlim=tcrop,ylim=dict(mag=[-150, 150])), topomap_args=dict(vmin=-50,vmax=50))
    
    evoked_targ_off.plot_joint(times=[.100,.180,.285,.490], ts_args=dict(xlim=tcrop,ylim=dict(mag=[-150, 150])), topomap_args=dict(vmin=-75,vmax=75))


#picks_misc = mne.pick_types(epochs_std.info, misc = True, meg = False, exclude='bads')


