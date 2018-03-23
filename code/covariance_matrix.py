# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:43:30 2017

@author: oussama.abdoun
"""
#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/MEG/code/")
from header import *
#==============================================================================


#mne.set_log_level(verbose = 'WARNING')


def covariance(task, states, subjects=None, tag='Cardiac', ECG_channel = ['EEG062-2800', 'EEG062'], do_empty_room_cov=False, do_epochs_cov=False):
    """Output (Covariance directory): *'-cov.fif'.
    Parameters:
        states: list of the states.
        subjects: list of the subjects, default to all subjects available for the task.
        tag: evoked series to process."""
    
    if not subjects:
        subjects = get_subjlist(task)
    
    data_path_noise = op.join(Raw_data_path, 'MEG_noise')
    
    l_freq = 0.1
    h_freq = 40
    baseline = (-.400, -.300)
    
    for subj in subjects:
        
        print((' - subj {} :'.format(subj)))   
                 
        data_path = op.join(Raw_data_path, task, 'meg')
        cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', get_id(subj))
        if not op.exists(cov_path):
            os.makedirs(cov_path)
        # EMPTY ROOM NOISE COVARIANCE
        if do_empty_room_cov:
            
            try:
                raw_empty_room = mne.io.read_raw_ctf(op.join(data_path_noise, * get_rawpath(subj, task = task, noise=1)), preload = True)
            except:
                print(colored('Did not find any noise file for subject {}'.format(subj), 'red'))
        
            else:
                picks_all = mne.pick_types(raw_empty_room.info,meg=True,eeg=True,stim=False,eog=True,exclude='bads')
                raw_empty_room = raw_empty_room.notch_filter(np.arange(50,251,50),picks=picks_all,n_jobs=6)
                raw_empty_room = raw_empty_room.filter(l_freq=l_freq,h_freq=h_freq, filter_length=20479, picks=picks_all, n_jobs=6)
                noise_cov = mne.compute_raw_covariance(raw_empty_room, tmin=0, tmax=None, n_jobs = 6)
                
                raw_empty_room.close()
                
                raw_cov_fname = op.join(cov_path, get_id(subj) + '_' + task + '_raw_covariance-cov.fif')
                mne.write_cov(raw_cov_fname, noise_cov)
                
        
        # NEEDS REWRITING
#==============================================================================
#       # EPOCHS BASELINE NOISE COVARIANCE
#==============================================================================
        if do_epochs_cov:
            
            
            sr_path = op.join(Analysis_path, task, 'meg', 'Source_Rec', get_id(subj))
            
            
            trans_paths = glob.glob(op.join(sr_path, '*-trans.fif'))
            trans_names = []
            
            for path in trans_paths:
                trans_names.append(op.split(path)[1].strip('-trans.fif').split('_')[2:])
        
            
            for state in states:
                
                print(('State : {}'.format(state))) 
                
                blocks = []
                
                for key in get_blocks(subj, task = task, state = state):
                    blocks.append(key)
                
                for i in range(len(trans_names)):
                    b = list(set(trans_names[i]) & set(blocks))
                    b.sort()
                    
                    if b:
                        epochs_list = []
                        a = ''
                        for blk in b:    
                
                            raw_fname = op.join(data_path, op.join(* get_rawpath(subj, task = task)) + blk + '.ds')
                            #load file
                            raw = mne.io.read_raw_ctf(raw_fname, preload=True)
                            picks_all = mne.pick_types(raw.info,meg=True,eeg=True,stim=False,eog=True,exclude='bads')
                            #Notch filter 50Hz + harmonics & band -------------------------------------------------------------------------------------
                            #Power line
                            filtered_raw = raw.notch_filter(np.arange(50,251,50),picks=picks_all,n_jobs=6)
                            #Bandpass
                            filtered_raw = filtered_raw.filter(l_freq=l_freq,h_freq=h_freq,picks=picks_all, n_jobs=6, fir_window='hann')
                            picks = mne.pick_types(filtered_raw.info,meg=True,exclude='bads')
                            try:
                                epochs = create_ecg_epochs(filtered_raw, tmin=-.5, tmax=.8, ch_name=ECG_channel[0], picks=picks, baseline=baseline)
                            except:
                                epochs = create_ecg_epochs(filtered_raw, tmin=-.5, tmax=.8, ch_name=ECG_channel[1], picks=picks, baseline=baseline)
                                
                            epochs_list.append(epochs)
                            a = a + '_' + blk
                        try:
                            epochs = mne.concatenate_epochs(epochs_list)
                        except:
                            epochs = epochs_list[0]
                    
            
                        epochs_noise_cov = mne.compute_covariance(epochs, n_jobs=6, method = 'empirical')
                        epochs_cov_fname = op.join(cov_path, get_id(subj) + ('_' + tag if tag else '') + ('_' + state if state else '') + '_epochs_covariance' + a + '-cov.fif')
                        mne.write_cov(epochs_cov_fname, epochs_noise_cov)
                        
            