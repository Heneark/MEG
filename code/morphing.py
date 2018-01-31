# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 15:24:26 2017

@author: oussama.abdoun
"""

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/MEG/code/")
from header import *
#==============================================================================


def average(task, subjects=None):
    """Output (Source_Rec directory): *'-fsaverage-rh.stc' + *'-fsaverage-rh.stc' per source estimate file processed.
    Parameters:
        subjects: list of the subjects, default to all subjects available for the task."""
        
    if not subjects:
        subjects = get_subjlist(task)
    
    for subj in subjects:
    
        
        subject_from = subj
        subject_to = 'fsaverage'
        
        stc_path = op.join(Analysis_path, task, 'meg', 'Source_Rec', get_id(subj), 'STC')
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
