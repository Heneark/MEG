# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 15:24:26 2017

@author: oussama.abdoun
"""

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from header import *
#==============================================================================

task = 'SMEG'
subjects = ['002','007','012','016','028','030','032','037','040','042']

time.sleep(7200)

for subj in subjects:

    
    subject_from = subj
    subject_to = 'fsaverage'
    
    sr_path = op.join(Analysis_path, task, 'meg', 'Source_Rec', get_id(subj))
    stc_path = op.join(sr_path, 'STC')
    fnames = glob.glob(op.join(stc_path, '*-rh.stc'))
    
    for fname in fnames:
        # Read input stc file
        stc_from = mne.read_source_estimate(fname)
        # Morph using one method (supplying the vertices in fsaverage's source
        # space makes it faster). Note that for any generic subject, you could do:
        #     vertices_to = mne.grade_to_vertices(subject_to, grade=5)
        # But fsaverage's source space was set up so we can just do this:
        vertices_to = [np.arange(10242), np.arange(10242)]
        stc_to = mne.morph_data(subject_from, subject_to, stc_from, n_jobs=6, grade=vertices_to)
        stc_to.save(fname[:-7] + '_fsaverage')
        
        
        # View source activations
        #plt.plot(stc_from.times, stc_from.data.mean(axis=0), 'r', label='from')
        #plt.plot(stc_to.times, stc_to.data.mean(axis=0), 'b', label='to')
        #plt.plot(stc_to_2.times, stc_to.data.mean(axis=0), 'g', label='to_2')
        #plt.xlabel('time (ms)')
        #plt.ylabel('Mean Source amplitude')
        #plt.legend()
        #plt.show()
