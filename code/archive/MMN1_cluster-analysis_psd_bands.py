# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:39:31 2016

@author: ousabd
"""

import mne
import numpy as np
import scipy as sp
import math, sys, os   
import matplotlib.pyplot as plt 
import matplotlib as mpl

# Add self path to searchable paths
# ! no need actually, it's done by default
sys.path.append(os.path.realpath(__file__))
# Import function to calculate matrix of sensors connectivity for cluster analysis
from GSN200_128chan_v21_neighbourhood import neighbourhood
# Import function to load data
from MMN1_io import MMN1_psd_load





def cmapLimits(evoked, margin, rounding, sym):
    if isinstance(evoked, list):
        mini = min([e.average().data.min() for e in evoked])
        maxi = max([e.average().data.max() for e in evoked])
    else:
        mini,maxi = evoked.data.min(), evoked.data.max()
    # Add requested margin
    mini,maxi = mini-margin*(maxi-mini), maxi+margin*(maxi-mini)
    # Round, if requested
    if rounding:    mini,maxi = round(mini), round(maxi)
    # Symmetrize if within-subject contrast
    if sym:
        mini = -max([abs(mini),abs(maxi)])
        maxi = abs(mini)
    
    return mini,maxi



# Define base path of data
pathBase = '/media/truecrypt1/pro/postdoc-lyon/projects/MMN1/data_preprocessed'
pathFig = '/media/truecrypt1/pro/postdoc-lyon/projects/MMN1/results'



allsubjs = dict(
    SC=['SC'+str(num) for num in [7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]],
    SP=['SP'+str(num) for num in [7,11,12,13,14,15,16,17,18,19,20,21,22,25,26,27]]
    )
    



#==============================================================================
# PARAMETERS
#==============================================================================
# Data selection and contrast
# Indicate the contrast using the list "groups":
# 1 group => one-sample test ; 2 grous => two-sample test
# Each group is a dictionary with the following mandatory keys:
# - "subjects": list of subjects in the form <group><number> (e.g. ["SP9","SP12"]
# - "states": a list containing 1 or 2 of "RE", "FA" and "OP"
# - "name": used for plot titles and file names

groups = list()
groups.append(dict(subjs=['SC7'], states=['RE'], name='SC'))
#groups.append(dict(subjs=allsubjs['SC']+allsubjs['SP'], states=['FA','OP'], name='FA > RE'))
#groups.append(dict(subjs=allsubjs['SC'], states=['RE'], name='SC'))
#groups.append(dict(subjs=allsubjs['SP'], states=['RE'], name='SP'))
name_general = 'RE'
ngroups = len(groups)


# Define frequency bands
bands_names = ['delta','theta','alpha','beta','gamma']
bands_freqs = dict(delta=[1,4], theta=[4,8], alpha=[8,12], beta=[15,25], gamma=[30,50])
nbands = len(bands_names)

fmin = 1
fmax = 60

blocks = [1,2,3]

# Load electrodes layout
montage = mne.channels.read_montage("/media/truecrypt1/pro/postdoc-lyon/projects/MMN1/meta/eloc_egi124.loc")
montage.ch_names = [s.strip('.') for s in montage.ch_names]

# Load matrix of sensors connectivity for cluster analysis
flayout = "/media/truecrypt1/pro/postdoc-lyon/projects/MMN1/meta/eloc_egi124.loc"
thres = 3.5
connectivity = neighbourhood(flayout,thres)

# Make-up of topomaps
topomap_args = dict(vmin=0, vmax=-math.log10(0.001), cmap='gist_heat_r')
mask_params = dict(marker='o', markerfacecolor='k', markersize=4, markeredgecolor='k', markeredgewidth=1)



#==============================================================================
# PLOT PSD
#==============================================================================
# ???




#==============================================================================
# SELECT DATA
#==============================================================================
# Load data
subjs = []
for g in groups:
    subjs += g['subjs']
data = MMN1_psd_load(pathBase,method='welch',ref='mastoids',subjs=subjs,fmin=fmin,fmax=fmax)


# Meta
shape = list(data['SC7']['RE'][1].data.shape[:2])
freqs = data['SC7']['RE'][1].freqs
ch_names = data['SC7']['RE'][1].info['ch_names']
nchan = len(ch_names)
sfreq = 1/(freqs[1]-freqs[0])
info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types='eeg', montage=montage)


# Arrange data
epochs = dict()
idxlim = dict()
for bn in bands_names:
    epochs[bn] = [[] for _ in groups]
    idxlim[bn] = list()
    for k in [0,1]:
        idxlim[bn].append(int((bands_freqs[bn][k]-fmin)*sfreq))

for g,group in enumerate(groups):
    for subj in group['subjs']:
        if len(group['states'])==1:
            tmp = np.mean(np.concatenate([data[subj][group['states'][0]][b].data for b in blocks], axis=2), axis=2)
        if len(group['states'])==2:
            tmp = [np.mean(np.concatenate([data[subj][group['states'][k]][b].data for b in blocks], axis=2), axis=2) for k in [0,1]]
#            tmp2 = np.mean(np.concatenate([data[subj][group['states'][0]][b].data for b in [0,1,2]], axis=2), axis=2)
            tmp = tmp[1]/tmp[0]
        if len(group['states'])==3:
            tmp = [np.mean(np.concatenate([data[subj][group['states'][k]][b].data for b in blocks], axis=2), axis=2) for k in [0,1,2]]
            tmp = np.mean(np.concatenate([tmp[k] for k in [0,1,2]], axis=2), axis=2)
            
        tmp = 10*np.log10(tmp)
            
        tmp = tmp.reshape(shape) # drop the 3rd dimension
    
        # Collapse into frequency bands
        for bn in bands_names:
            epochs[bn][g].append(tmp[:,idxlim[bn][0]:idxlim[bn][1]].mean(axis=1))
    
    # Convert to numpy arrays: n_epochs x n_channels
    for bn in bands_names:
        epochs[bn][g] = np.array(epochs[bn][g])

# Calculate colormap limits for topographies
if ngroups==1:
    vmin = min([epochs[bn][0].mean(axis=0).min() for bn in bands_names])
    vmax = max([epochs[bn][0].mean(axis=0).max() for bn in bands_names])
if ngroups==2:
    vmin = min([(epochs[bn][1].mean(axis=0)-epochs[bn][0].mean(axis=0)).min() for bn in bands_names])
    vmax = max([(epochs[bn][1].mean(axis=0)-epochs[bn][0].mean(axis=0)).max() for bn in bands_names])

margin = 0.1
vmin,vmax = vmin-margin*(vmax-vmin), vmax+margin*(vmax-vmin) # add margin
vmin,vmax = round(vmin), round(vmax) # round
vmin = -max([abs(vmin),abs(vmax)]) # symmetrize
vmax = abs(vmin) # symmetrize
        


#==============================================================================
# PERMUTATION-BASED CLUSTER DEFINITION
#==============================================================================
# Parameters
alpha = 0.05
p_sample = 0.05
#if ngroups==1:
t_sample = abs(sp.stats.t.ppf(p_sample/2,len(epochs[bands_names[0]][0])-1))
#else:
#    t_sample = None
threshold_tfce = dict(start=t_sample,step=0.2)
n_permutations = 10000
tfce = 1

# Prepare data
# ============
X = epochs

# Compute
# =======
clusters = dict()
cluster_pv = dict()
for b,bn in enumerate(bands_names):

    if tfce==1:
        threshold = threshold_tfce
    else:
        threshold = t_sample
        
    if ngroups==1:
        T_obs, clusters[bn], cluster_pv[bn], H0 = mne.stats.permutation_cluster_1samp_test(X[bn][0], threshold=threshold, n_permutations=n_permutations, connectivity=connectivity, t_power=1, out_type='mask', n_jobs=4)
    else:
        T_obs, clusters[bn], cluster_pv[bn], H0 = mne.stats.permutation_cluster_test(X[bn], threshold=threshold, n_permutations=n_permutations, connectivity=connectivity, t_power=1, out_type='mask', n_jobs=4)



#==============================================================================
# PLOT
#==============================================================================
# Prepare figure
# ==============
fig = plt.figure()
if ngroups==1:
    title = name_general+', ' + groups[0]['name']+', ' + str(n_permutations)+' permut., ' + r'$\alpha=$'+str(alpha)
elif ngroups==2:
    title = name_general+', ' + groups[1]['name']+' > '+groups[0]['name']+', ' + str(n_permutations)+' permut, ' + r'$\alpha=$'+str(alpha)
if tfce==1:
    title += ', TFCE (H=2, E=0.5, tmin='+'{:.2f}'.format(threshold_tfce['start']) + ', step='+'{:.1f}'.format(threshold_tfce['step']) + ')'
else:
    title += ', p='+str(p_sample)
 
# - Figure title
fig.suptitle(title)

# - Window title
fig.canvas.set_window_title(title)

fig = plt.gcf()
title = fig.canvas.get_window_title()
title = 'PSD-bands_' + title.replace(' > ','-').replace(', ','_').replace('\\','').replace('=','').replace(' ','').replace('$','')
fig.savefig(os.path.join(pathFig,title+'.png'), dpi=200)

# Plot topographies
# =================
for b,bn in enumerate(bands_names):
    # - Select significant clusters
    if tfce==1:
        nscl = sum((cluster_pv[bn] <= alpha))
        print("Number of significant points ["+bn+"]= "+str(nscl))
    else:
        significant_clusters = np.array(clusters[bn])[(cluster_pv[bn] <= alpha)]
        nscl = len(significant_clusters)
        print("Number of significant clusters ["+bn+"]= "+str(nscl))
        
    # - Plot power map    
    if ngroups==1:
        datapower = X[bn][0].mean(axis=0)
        axPower = fig.add_subplot(2,nbands+1,b+1)
        axPower.set_title(bn)
        im = mne.viz.plot_topomap(data=datapower, pos=info, axes=axPower, vmin=vmin, vmax=vmax, cmap='RdBu_r')
    if ngroups==2:
        datapower = X[bn][0].mean(axis=0)
        axPower = fig.add_subplot(2+2,nbands+1,b+1)
        axPower.set_title(bn)
        im = mne.viz.plot_topomap(data=datapower, pos=info, axes=axPower, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        datapower = X[bn][1].mean(axis=0)
        axPower = fig.add_subplot(2+2,nbands+1,nbands+1+b+1)
        axPower.set_title(bn)
        im = mne.viz.plot_topomap(data=datapower, pos=info, axes=axPower, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        datapower = X[bn][1].mean(axis=0)-X[bn][0].mean(axis=0)
        axPower = fig.add_subplot(2+2,nbands+1,2*(nbands+1)+b+1)
        axPower.set_title(bn)
        im = mne.viz.plot_topomap(data=datapower, pos=info, axes=axPower, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        
    

    # - Plot statistics provided by TFCE, only if any channel is significant
    if nscl:
        axStats = fig.add_subplot(2+ngroups*(ngroups==2),nbands+1,(nbands+1)*(1+ngroups*(ngroups==2))+b+1)
        axStats.set_title(bn)
        datastats = cluster_pv[bn]*(cluster_pv[bn] <= alpha) # threshold p-values at alpha
        datastats = -np.log10(datastats) # scale logarithmically
        datastats[datastats==np.inf] = 0 # correct infinite values
        mne.viz.plot_topomap(data=datastats, pos=info, axes=axStats, vmin=0, vmax=-np.log10(0.001), cmap='gist_heat_r', mask=(cluster_pv[bn] <= alpha), mask_params=mask_params)
     
# Plot colormaps
# ==============
# - For power
axPowerCb = fig.add_subplot(2,nbands+1,nbands+1)
axPowerCb.set_axis_off()
axPowerCb,_ = mpl.colorbar.make_axes(axPowerCb, location='left', shrink=0.6, aspect=12)
cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbPower = mpl.colorbar.ColorbarBase(axPowerCb, cmap=cmap, norm=norm)
cbPower.set_label('dB')

# - For statistics
if nscl:
    # - Plot colormap for stats
    axStatsCb = fig.add_subplot(2,nbands+1,2*(nbands+1))
    axStatsCb.set_axis_off()
    axStatsCb,_ = mpl.colorbar.make_axes(axStatsCb, location='left', shrink=0.6, aspect=12)
    cmap = mpl.cm.gist_heat_r
    norm = mpl.colors.Normalize(vmin=0, vmax=-np.log10(0.001))
    cbStats = mpl.colorbar.ColorbarBase(axStatsCb, cmap=cmap, norm=norm)
    cbStats.set_label('-log10(p)')

    
    
    
    
freqs = np.arange(1,51)
n_cycles = np.linspace(1,15,num=len(freqs))
power,itc = mne.time_frequency.tfr_morlet(data, freqs=freqs, n_cycles=n_cycles, return_itc=True, decim=1, n_jobs=4, average=False)