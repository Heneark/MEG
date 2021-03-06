# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 17:42:45 2017

@author: oussama.abdoun
"""

#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# # # # # # # # # # # # # # 
# # # # # # # # # # # # # # OUT OF DATE
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# 
#==============================================================================

#==============================================================================
#==============================================================================
# # Compute FORWARD SOLUTION
#==============================================================================
#==============================================================================

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from header import *
#==============================================================================

task = 'SMEG'
states = ['RS','OM','FA']
subjects = ['002', '007', '012', '016', '028', '030', '032', '037', '040', '042']

#==============================================================================
tag = 'Cardiac'
#==============================================================================

#==============================================================================
spacing = 'oct6'
vol = 0
#==============================================================================


#==============================================================================
compute_fwd = False
compute_inv = True
compute_stc = False
#==============================================================================

for subj in subjects:
    
    cov_path = op.join(Analysis_path, 'Covariance', get_id(subj))
    sr_path = op.join(Analysis_path, task, 'meg', 'Source_Rec', get_id(subj))
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', get_id(subj))
    
    src = mne.read_source_spaces(op.join(sr_path, subj +'_source_space_' + spacing + '-src.fif'))
    bem_sol = mne.read_bem_solution(op.join(sr_path, subj + '_bem_solution-sol.fif'))
    
    trans_paths = glob.glob(op.join(sr_path, '*-trans.fif'))
    trans_names = []
    
    for path in trans_paths:
        trans_names.append(op.split(path)[1].strip('-trans.fif').split('_')[2:]) 

    
    for state in states:

        blocks = []
        
        for key in get_blocks(subj, task = task, state = state):
            blocks.append(key)
        
        for i in range(len(trans_names)):
            b = list(set(trans_names[i]) & set(blocks))
            b.sort()
            
            if b:
                
                evoked_list = []
                a = ''
                for blk in b:
                   evoked_list.append(mne.read_evokeds(op.join(evoked_path, get_id(subj) + '_' + tag + '_' + (state if state else '') + '_evoked' + '_blk' + blk + '-ave.fif'), baseline = (-.400,-.300))[0])
                   a = a + '_' + blk
                evoked = mne.combine_evoked(evoked_list, 'nave')
                
                if compute_fwd == True:                
                    trans_path = trans_paths[i]
                    fwd = mne.make_forward_solution(evoked.info, trans=trans_path, src=src, bem=bem_sol, meg=True, eeg=False, mindist=5.0, n_jobs=2)
                    print(fwd)
                    mne.write_forward_solution(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') +'_fwd_solution' + a + '-fwd.fif'), fwd, overwrite = True)
    
                if compute_inv == True:
                    fwd = mne.read_forward_solution(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') +'_fwd_solution' + a + '-fwd.fif'))
                    fwd = mne.pick_types_forward(fwd, meg=True, eeg=False)
                    noise_cov = mne.read_cov(op.join(cov_path, get_id(subj) + '_' + task + '_raw_covariance-cov.fif'))
                    inv = make_inverse_operator(evoked.info, fwd, noise_cov, verbose='DEBUG')
                    write_inverse_operator(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') +'_inv_mne' + a + '-inv.fif'), inv)
                    
                if compute_stc == True:
                    method = "dSPM"
                    snr = 3.
                    lambda2 = 1. / snr ** 2
                    stc = apply_inverse(evoked, inv, lambda2, method=method, pick_ori=None)
                    stc.save(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') + '_source_estimate' + a))

                    
            
        
