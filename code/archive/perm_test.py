# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:44:15 2017

@author: oussama.abdoun
"""

from mne.stats import permutation_cluster_test
from mne.viz import plot_sensors

import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from csv_io_alex import *
from io_alex import *

import os
import os.path as op

import numpy as np

from matplotlib import pyplot as plt

from termcolor import colored

task = 'SMEG'


pathBase = '/dycog/meditation/ERC/Analyses'
data_path = op.join(pathBase, task, 'meg')

#blocks = ...
data = psd_load(subj, data_path, blocks, 'welch', fmin=1, fmax=60, norm=True)
info = data['06'].info
#==============================================================================
# use getblocks with state argument 'FA', 'FA1', or 'FA2' etc
#==============================================================================

def perm_test(state1, state2, threshold = .0005, n_permutations = 10000):
    
    X1 = np.zeros((state_count(state1), 283, 30))
    X2 = np.zeros((state_count(state2), 283, 30))
    i = 0
    j = 0
    for subj in get_subjlist(task):
        
        
        blocks1 = []
        for key in get_blocks(subj, state=state1):
            blocks1.append(key)
        
        blocks2 = []
        for key in get_blocks(subj, state=state2):
            blocks2.append(key)
            
            
        data1 = psd_load(subj, data_path, blocks1, 'welch', fmin=1, fmax=60, norm=True)
        data2 = psd_load(subj, data_path, blocks2, 'welch', fmin=1, fmax=60, norm=True)
        

        for b in blocks1:
            X1[i,:,:] = (np.squeeze(data1[b].data))[:283]
            i = i + 1

        for b in blocks2:
            X2[j,:,:] = (np.squeeze(data2[b].data))[:283]
            j = j + 1
            
    T_obs, clusters, cluster_p_values, H0 = permutation_cluster_test([X1, X2], n_permutations=n_permutations, threshold=threshold, tail=0)
    
    perm_test_results = {'T_obs' : T_obs, 'clusters' : clusters, 'cluster_p_values' : cluster_p_values, 'H0' : H0}
    return perm_test_results
       

