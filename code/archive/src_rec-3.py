# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 17:42:45 2017

@author: oussama.abdoun
"""

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


# Evoked parameters
#==============================================================================
l_freq = 1
h_freq = 40
tag = 'Cardiac'
#remove_ECG_artifacts = True
REM = [True, False]
# To remove all of them, set to 'None'
#number_of_eog_components_to_remove = None
NUM = [None]
#==============================================================================
#==============================================================================

# Source space parameters
#==============================================================================
spacing = 'oct6'
vol = 0
#==============================================================================

# What will it be today sir ?
#==============================================================================
compute_fwd = False
compute_inv = True
compute_stc = True
#==============================================================================

# Shall we compute the covariance matrix based on empty room noise or on epochs ?
#==============================================================================
#cov = 'raw'
cov = 'epochs'

covariance = ['raw', 'epochs']
#==============================================================================

for subj in subjects:
    
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', get_id(subj))
    sr_path = op.join(Analysis_path, task, 'meg', 'Source_Rec', get_id(subj))
    stc_path = op.join(sr_path, 'STC')
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', get_id(subj))
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', get_id(subj))
    src = mne.read_source_spaces(op.join(sr_path, subj +'_source_space_' + spacing + '-src.fif'))
    bem_sol = mne.read_bem_solution(op.join(sr_path, subj + '_bem_solution-sol.fif'))
    
    if not op.exists(stc_path):
        os.makedirs(stc_path)
    
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
                for remove_ECG_artifacts in REM:
                    for number_of_eog_components_to_remove in NUM:
                        for cov in covariance:                        
                        
                            evoked_list = []
                            a = '_blk'
                            for blk in b:
                                
                                ev = mne.read_evokeds(op.join(evoked_path, get_id(subj) + '_' + tag + '_' + (state if state else '') + '_evoked_{}-{}_blk'.format(l_freq,h_freq) + blk + '-ave.fif'), baseline = (-.400,-.300))[0]
                                
                                #==============================================================================
                                #load ICA file, and apply ICA
                                ica_fname = sorted(glob.glob(op.join(ICA_path, get_id(subj) + '_ICA_ncomp0.99*_blk' + blk + '-ica.fif')), key = op.getmtime)[-1]
                                ica = mne.preprocessing.read_ica(ica_fname)
                                
                                n_comp = ica.get_components().shape[1]
                                exclude_inds = []
                                
                                t = '_ECGincluded'
                                
                                if remove_ECG_artifacts == True:
                                    t = '_ECGexcluded'
                                    exclude_inds.extend(ica.labels_['ecg'])
                                    
                                if number_of_eog_components_to_remove:
                                    s = '_{}EOGcomp'.format(number_of_eog_components_to_remove)
                                    exclude_inds.extend(ica.labels_['eog'][:number_of_eog_components_to_remove])
                                else:
                                    s = '_{}EOGcomp'.format(len(ica.labels_['eog']))
                                    exclude_inds.extend(ica.labels_['eog'])
                                                    
                                include = list(set(range(n_comp)) - set(exclude_inds))
                                
                                ica.apply(ev, include = include)
                                #==============================================================================     
                                
                                evoked_list.append(ev)
                                a = a + '_' + blk
                            evoked = mne.combine_evoked(evoked_list, 'nave')
                             
                            if compute_fwd == True:                
                                trans_path = trans_paths[i]
                                fwd = mne.make_forward_solution(evoked.info, trans=trans_path, src=src, bem=bem_sol, meg=True, eeg=False, mindist=5.0, n_jobs=2)
                                print(fwd)
                                mne.write_forward_solution(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') +'_fwd_solution' + a[4:] + '-fwd.fif'), fwd, overwrite = True)
                
                            if compute_inv == True:
                                fwd = mne.read_forward_solution(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') +'_fwd_solution' + a[4:] + '-fwd.fif'))

                                fwd = mne.pick_types_forward(fwd, meg=True, eeg=False)
                                if cov == 'raw':
                                    c = '_raw-cov'
                                    noise_cov = mne.read_cov(op.join(cov_path, get_id(subj) + '_' + task + '_raw_covariance-cov.fif'))
                                if cov == 'epochs':
                                    c = '_epo-cov'
                                    noise_cov = mne.read_cov(op.join(cov_path, get_id(subj) + ('_' + tag if tag else '') + ('_' + state if state else '') + '_epochs_covariance' + a[4:] + '-cov.fif'))
                                inv = make_inverse_operator(evoked.info, fwd, noise_cov, verbose='DEBUG')
                                write_inverse_operator(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') +'_inv_mne' + c + a[4:] + '-inv.fif'), inv)
                                
                            if compute_stc == True:
                                if cov == 'raw':
                                    c = '_raw-cov'
                                if cov == 'epochs':
                                    c = '_epo-cov'
                                inv = read_inverse_operator(op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') +'_inv_mne' + c + a[4:] + '-inv.fif'))
                                method = "dSPM"
                                snr = 3.
                                lambda2 = 1. / snr ** 2
                                stc = apply_inverse(evoked, inv, lambda2, method=method, pick_ori=None)
                                stc.save(op.join(stc_path, subj + '_' + task + ('_vol' if vol else '') + '_' + tag + '_{}-{}'.format(l_freq,h_freq) + t + s + c + a))

                    
            
        
