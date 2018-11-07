# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 18:12:03 2017

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

subj = '081'

pathBase = '/dycog/meditation/ERC/Raw data/'

#ICA_path = '/dycog/meditation/ERC/Analyses/MIMOSA/meg/ICA/' + get_id(subj)
#n_components = 0.975
#listblk = ['{:02d}'.format(b) for b in xrange(1,10)]

            


# SUBJECTS
listSubj = ['040','028','030','042','053','052','063','064','090','096','050','081','093']
listSubj = ['052','064','090','096','050','081','093']

# DATAFRAME
df = pd.DataFrame(columns=['id','block','state','chan','freq','value'])


for subj in listSubj:
    
    #==============================================================================
    # PARAMETERS
    #==============================================================================
    
    # BLOCKS
    listblk = [k for k,v in get_blocks(subj).iteritems()]
    listblkFA = [k for k,v in get_blocks(subj).iteritems() if v=='FA']
    listblkOM = [k for k,v in get_blocks(subj).iteritems() if v=='OM']
    
    # PSD
    welch_width = 200 #12 for ring-rhythm peak
    welch_step = 50 #3 for ring-rhythm peak
    fmin = float(1)/welch_width
    fmax = 2
    
    # CHANNELS
    idxchan_EGG = [305,306,307,308,309,310,311,312]
    idxchan_Resp = [318]
    idxchan_ECG = [319]
    idxchan_allphysio = [305,306,307,308,309,310,311,312,318,319]
    

    #==============================================================================
    #     COMPUTATION
    #==============================================================================
    #psd_all = np.zeros((len(idxchan_EGG),int(float(fmax)/fmin)))
    
    list_psd = []
    
    for blk in listblk:
    
        filename = os.path.join(pathBase, get_task(subj,block=blk), 'meg', subj, get_rawpath(subj)[1] + blk + '.ds')
            
        raw = mne.io.read_raw_ctf(filename,preload=False)
        
        layout = mne.channels.find_layout(raw.info, ch_type='mag')
        picks_meg = mne.pick_types(raw.info, meg=True)
        
    #    ica = mne.preprocessing.read_ica(os.path.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components) + '_blk' + blk + '-ica.fif'))
        #raw = ica.apply(raw)  
        
        info = raw.info
        sfreq = info['sfreq']    
        
    #    mne.viz.plot_raw_psd(raw, fmin=0.0,fmax=2, n_fft=round(welch_width*sfreq),n_overlap=round((welch_width-welch_step)*sfreq), picks=idxchan_EGG)
    #    mne.viz.plot_raw_psd(raw, fmin=0.01,fmax=2, n_fft=round(welch_width*sfreq),n_overlap=round((welch_width-welch_step)*sfreq), picks=idxchan_Resp)    
    #    mne.viz.plot_raw_psd(raw, fmin=0.0,fmax=2, n_fft=round(welch_width*sfreq),n_overlap=round((welch_width-welch_step)*sfreq), picks=idxchan_ECG)
    
         
        psd, freqs = mne.time_frequency.psd_welch(raw, fmin=fmin,fmax=fmax, n_fft=round(welch_width*sfreq),n_overlap=round((welch_width-welch_step)*sfreq), picks=idxchan_EGG, n_jobs=4)
        list_psd.append(psd)
            
        #raw.save(os.path.join(pathBase,'egg_block'+'{:02d}'.format(b)+'_raw.fif'),picks=[305,306,307,308,309,310,311,312,318,319])
        
    #    psd = psd[:,:,np.newaxis]
    #    info_psd = info
    #    info_psd['chs'] = [info['chs'][i] for i in idxchan_EGG]
        
    # Average psd over all blocks
    psd_avg = 0*list_psd[0]
    for psd in list_psd:
        psd_avg += psd
    psd_avg /= len(list_psd)
    
    #==============================================================================
    # PLOT PSD OF EGG
    #==============================================================================
    if 0:
        # Large range
        fig = plt.figure('PSD EGG: '+'{:.3f}'.format(fmin)+' - '+'{:.3f}'.format(fmax) + ' Hz')
        for ch in xrange(psd_avg.shape[0]):
            plt.plot(freqs, 10*np.log10(psd_avg[ch]))
        
        # Short range
        fig = plt.figure('PSD EGG: '+'{:.3f}'.format(fmin)+' - '+'{:.3f}'.format(0.1) + ' Hz')
        for ch in xrange(psd_avg.shape[0]):
            plt.plot(freqs[:20], 10*np.log10(psd_avg[ch][:20,]))
        
    
    #==============================================================================
    # PLOT PSD OF MEG CHANNELS
    #==============================================================================  
    if 0:    
        fig = plt.figure('PSD MEG: '+'{:.3f}'.format(fmin)+' - '+'{:.3f}'.format(fmax) + ' Hz')
        for ch in xrange(psd_avg.shape[0]):
            plt.plot(freqs, 10*np.log10(psd_avg[ch]))    
        
        
        fig = plt.figure('PSD MEG')
        plt.plot(freqs, 10*np.log10(np.mean(psd_avg,axis=0)))
        
        # Zero-mean 
        psd_avg_0mean = 10*np.log10(psd_avg) - np.tile(np.mean(10*np.log10(psd_avg), axis=1),(psd_avg.shape[1],1)).swapaxes(0,1)
        
        fig = plt.figure('PSD MEG (zero-mean)')
        plt.plot(freqs, np.mean(psd_avg_0mean,axis=0))
    
    
    #==============================================================================
    # SAVE RESULTS OF EGG
    #==============================================================================
    for blk in xrange(len(list_psd)):
        state = get_blocks(subj)[listblk[blk]]
        psd = list_psd[blk]
        for ch in xrange(psd.shape[0]):
            for f in xrange(len(freqs[:20])):
                df.loc[len(df)] = [subj,listblk[blk],state,int(ch+1),freqs[f],10*np.log10(psd[ch][f])]
    
    df.chan = df.chan.astype(int)
    df.to_csv("/dycog/meditation/ERC/Analyses/MIMOSA/egg/results/EGG_data.csv", sep=" ", float_format="%.3f", index=False)
    
    ## Find peak location & amplitude for each channel
    #for ch in xrange(psd_avg.shape[0]):
    #    location = freqs[np.argmax(psd_avg[ch])]
    #    amplitude = np.max(freqs*psd_avg[ch])/location
        
                         
    
    #fig = plt.figure('PSD EGG: '+'{:.3f}'.format(fmin)+' - '+'{:.3f}'.format(fmax) + ' Hz')
    #for ch in xrange(len(idxchan_EGG)):
    #    plt.plot(freqs, 10*np.log10(psd_avg[ch]))
    #    
    #fig = plt.figure('PSD EGG: '+'{:.3f}'.format(fmin)+' - '+'{:.3f}'.format(20*fmin) + ' Hz')
    #for ch in xrange(len(idxchan_EGG)):
    #    plt.plot(freqs[:20], 10*np.log10(psd_avg[ch][:20]))