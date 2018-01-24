#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 17:22:34 2017

@author: alex
"""

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from header import *
#==============================================================================


#CONVERT NP_ARRAY TO AverageTFR
def cvt_np2atfr(data):
    
    data = data[:,:,np.newaxis]
    data = mne.time_frequency.AverageTFR(info=info, data=data, times=[0], freqs=freqs, nave=200, method='welch')
    
    return data
    
#CONTRAST BETWEEN BLOCKS (log) AND SAVE TOPOMAPS
def contrast(subj, n1, n2, topo=None, vmin=-0.4, vmax=0.4):
    cont = (np.log(data[n1].data)-np.log(data[n2].data)).squeeze()
    cont = cvt_np2atfr(cont)

    if topo == 'plot':
        for i in range(len(freq_bands)):
            im = cont.plot_topomap(fmin = freq_bands[i][0], fmax = freq_bands[i][1], title = freq_bands[i][2], vmin=vmin, vmax=vmax, size=7)
            im.savefig(op.join(contrast_path, get_id(subj), get_id(subj) + '_contrast_' + '_blk_' + n1 + '-' + n2 + '_' + dict_blocks[n1] + '-' + dict_blocks[n2] + '_' + freq_bands[i][2] + '.' + picformat),format=picformat)
            plt.close(im)
        return None
        
    else:
        return cont


def contrast_topomaps(subj, n1, n2, vmin=-0.4, vmax=0.4):
    
    
    contrast(subj, n1, n2, 'plot', vmin, vmax)

    fig = plt.figure(figsize = (18,6))
    plt.title(get_id(subj) + '_contrast_' + '_blk_' + n1 + '-' + n2 + '_' + dict_blocks[n1] + '-' + dict_blocks[n2])
    plt.axis('off')
    gs = gridspec.GridSpec(1, len(freq_bands))
    gs.update(wspace=0, hspace=0)
    
    for i in range(len(freq_bands)):
        fig.add_subplot(gs[i])
        img = mpimg.imread(op.join(contrast_path, get_id(subj), get_id(subj) + '_contrast_' + '_blk_' + n1 + '-' + n2 + '_' + dict_blocks[n1] + '-' + dict_blocks[n2] + '_' + freq_bands[i][2] + '.' + picformat))
        os.remove(op.join(contrast_path, get_id(subj), get_id(subj) + '_contrast_' + '_blk_' + n1 + '-' + n2 + '_' + dict_blocks[n1] + '-' + dict_blocks[n2] + '_' + freq_bands[i][2] + '.' + picformat))        
        plt.imshow(img)
        plt.axis('off')

    fig.savefig(op.join(contrast_path, get_id(subj), get_id(subj) + '_contrast_' + '_blk_' + n1 + '-' + n2 + '_' + dict_blocks[n1] + '-' + dict_blocks[n2] + '.svg'),format='svg')



# Run and save all contrast topomaps for all subjects for a given task
#==============================================================================
task = 'SMEG'

subjects = get_subjlist(task)

freq_bands = [[8, 12, 'alpha'], [12, 25, 'beta'], [25, 40, 'gamma']]
picformat = 'png'


data_path = op.join(Analysis_path, task, 'meg')
contrast_path = op.join(data_path, 'Contrast')

for subj in subjects:
    
    data = None
    
    dict_blocks = get_blocks(subj)
    
    blocks = []
    
    for key in get_blocks(subj):
        blocks.append(key)

    
    if not op.exists(op.join(contrast_path, get_id(subj))):
        os.makedirs(op.join(contrast_path, get_id(subj)))
    
    
    data = psd_load(subj, data_path, blocks, 'welch', fmin=1, fmax=60, norm=True)
            
    info = data[blocks[0]].info
    freqs = data[blocks[0]].freqs
    
    for i in range(len(blocks)):
        for j in range(i):
            n1 = blocks[i]
            n2 = blocks[j]
            contrast_topomaps(subj, n1, n2)
            plt.close('all')
#==============================================================================


