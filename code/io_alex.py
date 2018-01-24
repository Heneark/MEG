# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 13:30:06 2016

@author: ousabd
"""


import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from csv_io_alex import *
import mne
import numpy as np
import os
import os.path as op
import sip
import pandas as pd
sip.setapi('QString', 2)
#os.path.basename(unicode(p))


# Define general structure of data filename
def filename(pathIn,block,tag=None):
    return op.join(pathIn, 'Epochs', get_id(subj), get_id(subj) + ('_' + tag if tag else '') + '_epochs' + '_blk' + block + '-epo.fif')  

h_freq = 60
fmin = 1
fmax = 60
fstep = 0.5
freqs = np.arange(fmin,fmax+fstep,fstep)


def erp_save(subj,pathBase,blocks,tag=None,h_freq=60):

    for block in blocks:
        epoch = mne.read_epochs(filename(pathBase,block,tag))
        epoch.drop_bad()
        # filter, if requested
        if h_freq:  epoch.savgol_filter(h_freq, copy=False)
        data = epoch
        
        # Save evoked
        fname = op.join(pathBase, 'Evoked', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_evoked' + '_blk' + block + '_filt1-'+str(h_freq if h_freq else 60)+'Hz-ave.fif')
        data.average().save(fname)


def erp_load(subj,pathBase,blocks,tag=None,h_freq=60):
    
    data = dict()
    
    for block in blocks:
        fname = op.join(pathBase, 'Evoked', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_evoked' + '_blk' + block + '_filt1-'+str(h_freq if h_freq else 60)+'Hz'+'-ave.fif')
        data[block] = mne.read_evokeds(fname, verbose=False)[0]
                
    return data
            
            
#def cov_save(subj,pathBase,blocks,tag=None):
#
#    for block in blocks:
#        epoch = mne.read_epochs(filename(pathBase,block,tag))
#        epoch.drop_bad()
#        data = epoch
#        
#        # Compute and save noise covariance matrix
#        fname = op.join(pathBase, 'Covariance', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_covariance' + '_blk' + block + '-cov.fif')
#        noise_cov_baseline = mne.compute_covariance(data)
#        mne.write_cov(fname, noise_cov_baseline)
            
def cov_load(subj,pathBase,blocks,tag=None):
    
    data = dict()
    
    for block in blocks:
        fname = op.join(pathBase, 'Covariance', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_epochs_covariance' + '_blk' + block + '-cov.fif')
        data[block] = mne.read_cov(fname, verbose=False)
                
    return data
                      


def psd_save(subj,pathBase,blocks,method,tag=None,fmin=1,fmax=60,t_fft=0.5,t_overlap=0.25):
#==============================================================================
# Pour chaque électrode, normaliser par la puissance spectrale totale (pour cette électrode)
# DSPn(el1,freq1) = DSP(el1,freq1) / DSP(el1)
#==============================================================================
    info = []
    for block in blocks:
        
        # Load epochs data
        data = mne.read_epochs(filename(pathBase,block,tag))
        # Drop bad epochs
        data.drop_bad()
        
        data.pick_types(meg=True , eeg=True, ref_meg=False, stim=False, misc=False)
        
        #if not info:
        info = data.info
            
        # Concatenate deviant and standard epochs
        #data = mne.concatenate_epochs([data['de'],data['st']])
        nave = data.get_data().shape[0]

        # Compute PSD, according to specified method
        if method=='multitaper':
            data,freqs = mne.time_frequency.psd_multitaper(data, fmin=fmin, fmax=fmax)
        if method=='welch':
            sfreq = data.info['sfreq']
            n_fft = sfreq*t_fft
            n_overlap = sfreq*t_overlap
            data,freqs = mne.time_frequency.psd_welch(data, fmin=fmin, fmax=fmax, n_fft=n_fft, n_overlap=n_overlap, picks=None)
                    
        # Average over epochs
        data = data.mean(axis=0)
        
        # Create a normalized set of data 'ndata'
        ndata = data
        
        for i in range(ndata.shape[0]):
            totalpsd = sum(ndata[i,:])
            ndata[i,:] = ndata[i,:]/totalpsd
        
        # Store in a AverageTFR object (add an empty time axis)
        data = data[:,:,np.newaxis]
        data = mne.time_frequency.AverageTFR(info=info, data=data, times=[0], freqs=freqs, nave=nave, method=method)
        ndata = ndata[:,:,np.newaxis]
        ndata = mne.time_frequency.AverageTFR(info=info, data=ndata, times=[0], freqs=freqs, nave=nave, method=method)
        # Save
        fname = op.join(pathBase, 'Psds', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_psd_' + method + '_blk' + block + '_fmin'+str(fmin)+'_fmax'+str(fmax)+'-tfr.h5')
        data.save(fname, overwrite=True)
        fname = op.join(pathBase, 'Psds', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_normalized_psd_' + method + '_blk' + block + '_fmin'+str(fmin)+'_fmax'+str(fmax)+'-tfr.h5')
        data.save(fname, overwrite=True)    
            
            
def psd_load(subj,pathBase,blocks,method, tag=None,fmin=1,fmax=60,norm=True):
        
    data = dict()
    for block in blocks:
        
        if norm == True:
            fname = op.join(pathBase, 'Psds', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_normalized_psd_' + method + '_blk' + block + '_fmin'+str(fmin)+'_fmax'+str(fmax)+'-tfr.h5')
            data[block] = mne.time_frequency.read_tfrs(fname)[0]
        if norm == False:
            fname = op.join(pathBase, 'Psds', get_id(subj), get_id(subj) +('_' + tag if tag else '')+ '_psd_' + method + '_blk' + block + '_fmin'+str(fmin)+'_fmax'+str(fmax)+'-tfr.h5')
            data[block] = mne.time_frequency.read_tfrs(fname)[0]
            
    return data




# SAVE EVOKED AND/OR PSDS FOR ALL/SOME SUBJECTS
#==============================================================================
#
## Revoir task et protocol (protocol>task and not the other way around)
#
## Save evoked ?
#do_evoked = True
## Save psds ?
#do_psds = False
#
## For which subjects ?
##subjects = get_subjlist()
#subjects = ['012']
#
## For which protocol ?
#protocol = 'WIM'
#
## The tag of the epochs to use (if none, tag='None')
#tag = 'testWIM'
#
#for subj in subjects:
#    
#    pathBase = op.join('/dycog/meditation/ERC/Analyses', get_task(subj,protocol = protocol), 'meg')
#    
#    blocks = []
#    
#    for key in get_blocks(subj, protocol=protocol):
#        blocks.append(key)
#        
#    if do_evoked:
#        evoked_path = op.join(pathBase, 'Evoked', get_id(subj))
#        if not op.exists(evoked_path):
#            os.makedirs(evoked_path)
#        
#        erp_save(subj,pathBase,blocks,tag)
#        
#        
#    if do_psds:
#        psds_path = op.join(pathBase, 'Psds', get_id(subj))
#        if not op.exists(psds_path):
#            os.makedirs(psds_path)
#        
#        psd_save(subj,pathBase,blocks,'welch',tag)
#    
#==============================================================================
    
    
