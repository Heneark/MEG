# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:22:13 2017

@author: alex
"""
#IMPORT PACKAGES, DEFINE ENVIRONMENT VARIABLES, CREATE PATHES
#==============================================================================
from header import *
#==============================================================================


# MANUAL EXECUTION
#==============================================================================
# subject='094'; task='SMEG'; state='OM'; block='15'; n_components=.975; method='fastica'; notch=np.arange(50,301,50); high_pass=0.1; low_pass=None; rejection={'mag':2.5e-12}; ica_rejection={'mag':4000e-15}; ECG_channel=['EEG062-2800', 'EEG062']; EOG_channel='EOGV'; stim_channel='UPPT001'

# ica = run_ica(task, subject, state, block, save=False, ica_rejection={'mag':4000e-15}, ECG_threshold=0.25, EOG_threshold=3.5)

# raw, raw_ECG = process(task, subject, state, block, check_ica=True, save_ica=False, overwrite_ica=False, notch=np.arange(50,301,50), high_pass=1, low_pass=40)

# epochs,evoked = epoch(task, subject, state, block, save=False, rejection={'mag':2.5e-12}, tmin=-.5, tmax=.8, baseline=(-.4,-.3), overwrite_ica=False, ica_rejection={'mag':4000e-15}, notch=np.arange(50,301,50), high_pass=0.5, low_pass=None, ECG_threshold=0.25, EOG_threshold=3.5)
#==============================================================================


def run_ica(task, subject, state, block, raw=None, save=True, fit_ica=False, n_components=0.975, method='fastica', ica_rejection={'mag':4000e-15}, ECG_channel=['EEG062-2800', 'EEG062'], ECG_threshold=0.25, EOG_channel='EOGV', EOG_threshold=3.5, stim_channel='UPPT001'):
    """
    Fit ICA on raw MEG data and return ICA object.
    If save, save ICA, save ECG and EOG artifact scores plots (deleting previously existing ones), and write log (default to True).
    If fit_ica, fit ICA even if there is already an ICA file (default to False).
    Output:
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-ica.fif'
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_ecg.svg'
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_eog.svg'
    Log:
        'Analyses/<task>/meg/ICA/ICA_log.tsv'
    Parameters (see mne.preprocessing.ICA):
        raw: raw data to fit ICA on. If None (default), will be loaded according to previous parameters.
        n_components: number of components used for ICA decomposition
        method: the ICA method to use
        ica_rejection: epoch rejection threshold (default to 4000 fT for magnetometers)
        ECG_threshold: ECG artifact detection threshold (mne default to 0.25; may be too low)
        ECG_threshold: ECG artifact detection threshold (mne default to 3.0, increased to 3.5)
        ECG_channel: channel or list of possible channels corresponding to the ECG
        EOG_channel: channel corresponding to the EOG
        stim_channel: channel corresponding to start and end experiment events
    """
    # Load data
    data_path = op.join(Raw_data_path, task, 'meg')
    raw_fname = op.join(data_path, op.join(* get_rawpath(subject, task=task)) + block + '.ds')
    if not raw:
        raw = mne.io.read_raw_ctf(raw_fname, preload=True)
    
    # ICA path
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    if save and not op.isdir(ICA_path):
        os.makedirs(ICA_path)
    ICA_file = op.join(ICA_path, '{}_{}-{}_components-ica.fif'.format(state, block, n_components))
    
    # ICA log
    ICA_log = op.join(Analysis_path, task, 'meg', 'ICA', 'ICA_log.tsv')
    if save and not op.isfile(ICA_log):
        with open(ICA_log, 'w') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('date','time','subject','state','block','start','end','n_components','n_selected_comps','ncomp_ECG','pulse','ncomp_EOG','rejection','dropped_epochs'))
    
    # Crop recording
    raw.set_channel_types({stim_channel: 'stim'})
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[1][0]] > 300 else None
    raw.crop(tmin=start, tmax=end)
    
    # Channels fix
    if type(ECG_channel) is list:
        ECG_list = ECG_channel
        for chan in ECG_list:
            if chan in raw.ch_names:
                ECG_channel = chan
    
    raw.pick_types(meg=True, exclude='bads', include=[ECG_channel, EOG_channel], ref_meg=False)
    
    # Filter for ICA
    raw.filter(l_freq=1, h_freq=40, fir_design='firwin', n_jobs=4)
    
    # Fit ICA
    if fit_ica or not op.isfile(ICA_file):
        ica = ICA(n_components=n_components, method=method) #create ICA object
        ica.drop_inds_ = []
        ica.fit(raw, reject=ica_rejection, decim=6, picks=mne.pick_types(raw.info, meg=True)) #decimate: 200Hz is more than enough for ICA, saves time; picks: fit only on MEG
        dropped = len(ica.drop_inds_)
    else:
        ica = read_ica(ICA_file)
        ica_rejection = '{previous threshold}'
        dropped = 'previously rejected'
    
    ica.scores_ = dict()
    ica.labels_ = dict()
    ica.exclude = []
    
    # Detect ECG and EOG artifacts
    ica.scores_['ecg'] = ica.find_bads_ecg(raw, ch_name=ECG_channel, threshold=ECG_threshold)[1]
    pulse = mne.preprocessing.find_ecg_events(raw, l_freq=8, h_freq=16, ch_name=ECG_channel)[2]
    ica.exclude.extend(ica.labels_['ecg'])
    
    ica.scores_['eog'] = ica.find_bads_eog(raw, ch_name=EOG_channel, threshold=EOG_threshold)[1]
    ica.exclude.extend(ica.labels_['eog'])
    
    # Plot scores
    ica.plot_scores(ica.scores_['ecg'], exclude=ica.labels_['ecg'], title="ECG artifacts")
    if save:
        plt.savefig(op.join(ICA_path, '{}_{}-{}_components-scores_ecg.png'.format(state, block, n_components)), dpi=360, transparent=True)
        plt.close()
    
    ica.plot_scores(ica.scores_['eog'], exclude=ica.labels_['eog'], title="EOG artifacts")
    if save:
        plt.savefig(op.join(ICA_path, '{}_{}-{}_components-scores_eog.png'.format(state, block, n_components)), dpi=360, transparent=True)
        plt.close()
    
    # Save ICA
    if save:
        ica.save(ICA_file)
        # Write ICA log
        with open(ICA_log, 'a') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(time.strftime('%Y_%m_%d\t%H:%M:%S',time.localtime()),subject,state,block,int(round(start)),int(round(end)) if end else int(round(raw.times[-1])),n_components,ica.n_components_,len(ica.labels_['ecg']),int(round(pulse)),len(ica.labels_['eog']),ica_rejection,dropped))
    
    return ica


def process(task, subject, state, block, n_components=.975, ica=None, check_ica=True, save_ica=True, overwrite_ica=False, fit_ica=False, ica_rejection={'mag':4000e-15}, notch=np.arange(50,301,50), high_pass=0.5, low_pass=None, ECG_channel=['EEG062-2800', 'EEG062'], ECG_threshold=0.25, EOG_channel='EOGV', EOG_threshold=3.5, stim_channel='UPPT001'):
    """
    Run preprocessing and return preprocessed raw data.
    If check_ica, plot overlay and properties of ECG and EOG components (default to True).
    If save_ica, save ICA and ICA plots (deleting previously existing ones), and write ICA log (default to True).
    If overwrite_ica, artifact detection will be perfomed again even if an ICA file exists (default to False).
    Output:
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_ecg.svg'
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_eog.svg'
    Parameters (see mne.filter):
        ica: ICA object. If None (default), will be loaded according to previous parameters. ICA parameters:
            ica_rejection: epoch rejection threshold (default to 4000 fT for magnetometers)
            ECG_threshold: ECG artifact detection threshold (mne default to 0.25)
            ECG_threshold: ECG artifact detection threshold (mne default to 3.0, increased to 3.5)
        low_pass: frequency (in Hz) for low-pass filtering
        high_pass: frequency (in Hz) for high-pass filtering
        notch: frequency (in Hz) or list of frequencies to notch filter
        ECG_channel: channel or list of possible channels corresponding to the ECG
        EOG_channel: channel corresponding to the EOG
        stim_channel: channel corresponding to start and end experiment events
    """
    # Load data
    data_path = op.join(Raw_data_path, task, 'meg')
    raw_fname = op.join(data_path, op.join(* get_rawpath(subject, task=task)) + block + '.ds')
    raw = mne.io.read_raw_ctf(raw_fname, preload=True)
    
    # Load ICA
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    ICA_file = op.join(ICA_path, '{}_{}-{}_components-ica.fif'.format(state, block, n_components))
    if not ica:
        if overwrite_ica or not op.isfile(ICA_file):
            ica = run_ica(task, subject, state, block, raw=raw.copy(), n_components=n_components, ECG_channel=ECG_channel, EOG_channel=EOG_channel, stim_channel=stim_channel, save=save_ica, fit_ica=fit_ica, ica_rejection=ica_rejection, ECG_threshold=ECG_threshold, EOG_threshold=EOG_threshold)
        else:
            ica = read_ica(ICA_file)
    
    # Crop recording
    raw.set_channel_types({stim_channel: 'stim'})
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[-1][0]] > 300 else None
    raw.crop(tmin=start, tmax=end)
    
    # Channels fix
    if type(ECG_channel) is list:
        ECG_list = ECG_channel
        for chan in ECG_list:
            if chan in raw.ch_names:
                ECG_channel = chan
    
    raw.pick_types(meg=True, exclude='bads', include=[ECG_channel, EOG_channel], ref_meg=False)
    
    # Filter
    raw.notch_filter(notch, fir_design='firwin', n_jobs=4)
    raw.filter(l_freq=high_pass, h_freq=low_pass, fir_design='firwin', n_jobs=4)
    
    # Visual check
    if check_ica:
        # ECG components
        check_ecg = create_ecg_epochs(raw, ch_name=ECG_channel, reject=ica_rejection)
        ica.plot_overlay(check_ecg.average(), exclude=ica.labels_['ecg'])
        if save_ica:
            for pic in glob.glob(op.join(ICA_path, state+'*'+block+'*overlay_ecg*')):
                os.remove(pic)
            plt.savefig(op.join(ICA_path, '{}_{}-{}_components-overlay_ecg.png'.format(state, block, n_components)), dpi=360, transparent=True)
            plt.close()
        
        for comp in ica.labels_['ecg']:
            ica.plot_properties(check_ecg, picks=comp)
            if save_ica:
                for pic in glob.glob(op.join(ICA_path, state+'*'+block+'*properties_ecg*')):
                    os.remove(pic)
                plt.savefig(op.join(ICA_path, '{}_{}-{}_components-properties_ecg{}.png'.format(state, block, n_components, comp)), dpi=360, transparent=True)
                plt.close()
        
        # EOG components
        check_eog = create_eog_epochs(raw, ch_name=EOG_channel, reject=ica_rejection)
        ica.plot_overlay(check_eog.average(), exclude=ica.labels_['eog'])
        if save_ica:
            for pic in glob.glob(op.join(ICA_path, state+'*'+block+'*overlay_eog*')):
                os.remove(pic)
            plt.savefig(op.join(ICA_path, '{}_{}-{}_components-overlay_eog.png'.format(state, block, n_components)), dpi=360, transparent=True)
            plt.close()
        
        for comp in ica.labels_['eog']:
            ica.plot_properties(check_eog, picks=comp)
            if save_ica:
                for pic in glob.glob(op.join(ICA_path, state+'*'+block+'*properties_eog*')):
                    os.remove(pic)
                plt.savefig(op.join(ICA_path, '{}_{}-{}_components-properties_eog{}.png'.format(state, block, n_components, comp)), dpi=360, transparent=True)
                plt.close()
    
    # Apply ICA
    raw_ECG = raw.copy()
    ica.apply(raw)
    ica.apply(raw_ECG, exclude=ica.labels_['eog'])
    
    return raw, raw_ECG


def epoch(task, subject, state, block, raw=None, save=True, rejection={'mag':2.5e-12}, name=['ECG_included','ECG_excluded'], tmin=-.5, tmax=.8, baseline=(-.4,-.3), ECG_channel=['EEG062-2800', 'EEG062'], EOG_channel='EOGV', overwrite_ica=False, fit_ica=False, ica_rejection={'mag':4000e-15}, notch=np.arange(50,301,50), high_pass=0.5, low_pass=None, ECG_threshold=0.25, EOG_threshold=3.5):
    """
    Epoch preprocessed raw data and average to evoked response, and return them.
    Output:
        'Analyses/<task>/meg/Epochs/<subject>/<epoch_name>-<state>_<block>-epo.fif'
        'Analyses/<task>/meg/Evoked/<subject>/<epoch_name>-<state>_<block>-ave.fif'
    Parameters (see mne.Epochs):
        raw: preprocessed raw data to epoch. If None (default), will be loaded according to previous parameters.
            Raw parameters:
                low_pass: frequency (in Hz) for low-pass filtering
                high_pass: frequency (in Hz) for high-pass filtering
                notch: frequency (in Hz) or list of frequencies to notch filter
            ICA parameters:
                ica_rejection: epoch rejection threshold (default to 4000 fT for magnetometers)
                ECG_threshold: ECG artifact detection threshold (mne default to 0.25)
                ECG_threshold: ECG artifact detection threshold (mne default to 3.0, increased to 3.5)
        rejection: epoch rejection threshold (default to 2500 fT for magnetometers)
        name: list of events to epoch (only 'ECG_included' and 'ECG_excluded' are supported yet)
        tmin, tmax, baseline: epoching parameters
    """
    # Load data
    raw_path = op.join(Analysis_path, task, 'meg', 'Preprocessed_raw', subject)
    raw_file = op.join(raw_path, '{}_{}-raw.fif'.format(state, block))
    if not raw:
        raw, raw_ECG = process(task, subject, state, block, overwrite_ica=overwrite_ica, fit_ica=fit_ica, ica_rejection=ica_rejection, notch=notch, high_pass=high_pass, low_pass=low_pass, ECG_threshold=ECG_threshold, EOG_threshold=EOG_threshold)
    
    # Channels fix
    if type(ECG_channel) is list:
        ECG_list = ECG_channel
        for chan in ECG_list:
            if chan in raw.ch_names:
                ECG_channel = chan
    
    epochs = dict()
    evoked = dict()
    
    for epo in name:
        
        # Save paths
        epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs', subject)
        if not op.exists(epochs_path):
            os.makedirs(epochs_path)
        epochs_file = op.join(epochs_path, '{}-{}_{}-epo.fif'.format(epo, state, block))
        
        evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
        if not op.exists(evoked_path):
            os.makedirs(evoked_path)
        evoked_file = op.join(evoked_path, '{}-{}_{}-ave.fif'.format(epo, state, block))
        
        # Epoch
        if epo == 'ECG_included':
            epochs[epo] = create_ecg_epochs(raw_ECG, reject=rejection, tmin=tmin, tmax=tmax, baseline=baseline, ch_name=ECG_channel)
        if epo == 'ECG_excluded':
            epochs[epo] = create_ecg_epochs(raw, reject=rejection, tmin=tmin, tmax=tmax, baseline=baseline, ch_name=ECG_channel)
        
        # Rejection
        epochs[epo].plot_drop_log()
        plt.savefig(op.join(epochs_path, '{}-{}_{}-drop_log.png'.format(epo, state, block)), dpi=360, transparent=True)
        plt.close()
        epochs[epo].drop_bad()
        
        # Save epochs
        if save:
            epochs[epo].save(epochs_file)
        
        # Save evoked
        evoked[epo] = epochs[epo].average()
        if save:
            evoked[epo].save(evoked_file)
        
        drop_log = op.join(Analysis_path, task, 'meg', 'Epochs', 'drop_log.txt')
        with open(drop_log, 'a') as fid:
            fid.write('{} {} epochs dropped\t{}\n'.format(epochs_file.split('/')[-2:], len(np.array(epochs[epo].drop_log)[np.where(epochs[epo].drop_log)]), rejection))
        
    return epochs,evoked


def process_alex(tasks, states, subjects=None, run_ICA=True, n_components=0.975, ECG_channel=['EEG062-2800', 'EEG062'], EOG_channel='EOGV', save_evokeds = True, ecg_epochs = True):
    """If run_ICA, outputs *'-ica.fif' along with three SVG for checking in the ICA directory.
    If save_evokeds, outputs *'-ave.fif' in the Evoked directory.
    Epoching parameters to be adjusted within the function code.
    Parameters:
        tasks: single task or list of the tasks to process.
        states: list of the states to process.
        subjects: list of the subjects to process, default to all subjects available for the task."""
    
    picformat = 'svg'
    
    if save_evokeds:
        if not ecg_epochs:
            # If resting state, 'event_id = None'
            event_id1 = None
            event_id2 = None
        tmin = -1.0
        tmax = 1.0
        baseline = (-0.4, -0.3)
        epochs_reject = dict(mag=4e-12)
        # Give a name tag to your epochs file (or None)
        epochs_name = 'Cardiac'
        # Pass-band filter between h_freq and l_freq
        h_freq = 40
        l_freq = 1
        resample_freq = 200
    
    if type(tasks) is not list:
        tasks = [tasks]
    for task in tasks:
    
        if not op.exists(op.join(Analysis_path, task, 'meg', 'ICA')):
            os.makedirs(op.join(Analysis_path, task, 'meg', 'ICA'))
        
        if run_ICA:
            ICA_csv_filename = op.join(Analysis_path, task, 'meg', 'ICA', 'ICA_info' + time.strftime('_%y-%m-%d_%Hh%M',time.localtime())  + '.csv')
            
            with open(ICA_csv_filename, 'a') as fid:
                fid.write("{}\t{}\t{}\t{}\t{}\t{}\n".format('subject','state','block','ncomponents','ncomp_ECG','ncomp_EOG'))
        
        if not subjects:
            subjects = get_subjlist(task)
        
        
        for subj in subjects:
            
            for state in states:
                
                data_path = op.join(Raw_data_path, task, 'meg')
                raw = None
                epochs = None
                ica = None
                #==============================================================================
                # Paths to PREPROCESSED DATA
                #==============================================================================
            
                ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', get_id(subj))
                evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', get_id(subj))
                epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs', get_id(subj))
                
       
                #create the corresponding directories if they don't exist
                
                if not op.exists(ICA_path):
                    os.makedirs(ICA_path)
                
                if save_evokeds == True:
                    if not op.exists(evoked_path):
                        os.makedirs(evoked_path)
                    
                    if not op.exists(epochs_path):
                        os.makedirs(epochs_path)
                
                blocks = []
                
        #            IF STATE, get_blocks(..., state=state)
                for key in get_blocks(subj, task = task, state=state):
                        blocks.append(key)
                
                
                for blk in blocks:
                    
                    #file name
                    raw_fname = op.join(data_path, op.join(* get_rawpath(subj, task = task)) + blk + '.ds')
                    
                    #load file
                    raw = mne.io.read_raw_ctf(raw_fname, preload=True)
                    #print(raw.info)
                    
                    
                    #load bad channels
                    raw.load_bad_channels(op.join(raw_fname, 'BadChannels'))
                    
                    #define picks
                    picks_all = mne.pick_types(raw.info,meg=True,eeg=True,stim=False,eog=True,exclude='bads')
                    
                    picks_meg = mne.pick_types(raw.info, ref_meg=False, meg=True, eeg=False, eog=False, stim=False, exclude='bads')
                    
                    if save_evokeds:
                        #Notch filter 50Hz + harmonics & band -------------------------------------------------------------------------------------
                        #Power line
                        filtered_raw = raw.notch_filter(np.arange(50,251,50),picks=picks_all,n_jobs=6)
                        #Bandpass
                        filtered_raw = filtered_raw.filter(l_freq=l_freq,h_freq=h_freq,picks=picks_all, n_jobs=6, fir_window='hann')
                        #Resample
                        filtered_raw.resample(resample_freq, window = 'hann', n_jobs=6)
                    
                    
                    
                    
                    
                    #ICA--------------------------------------------------------------------------------------------------------------------
                    if run_ICA:
                        
                        rawforica = raw.copy()
                        
                        #Filter for ICA 
                        rawforica = rawforica.filter(l_freq=1, h_freq=40, l_trans_bandwidth=0.5, h_trans_bandwidth=0.5, filter_length = '10s', phase='zero-double',n_jobs=6)
                        
                        method = 'fastica'  
                        decim = 3  
                        
                        
                        ica = ICA(n_components=n_components, method=method)
                        print(ica)
                        
                        try:
                            reject = dict(mag=6e-12, grad=4000e-13)
                            ica.fit(rawforica, picks=picks_meg, decim=decim, reject=reject)
                            print(ica)
                        except:
                            reject = dict(mag=17e-12, grad=4000e-13)
                            ica.fit(rawforica, picks=picks_meg, decim=decim, reject=reject)
                            print(ica)
                        
                        
                        
                        title = 'Sources related to %s artifacts (red)'
                        
                        #ECG--------------------------------------                
                        
    #                    n_max_ecg = 3  
                        
                        try:
                            try:
                                ecg_epochs = create_ecg_epochs(rawforica, tmin=-.5, tmax=.5, ch_name=ECG_channel[0], picks=picks_meg)
                            except:
                                ecg_epochs = create_ecg_epochs(rawforica, tmin=-.5, tmax=.5, ch_name=ECG_channel[1], picks=picks_meg)
                            
                            ecg_inds, scores = ica.find_bads_ecg(ecg_epochs, method='ctps')
                            fig = ica.plot_scores(scores, exclude=ecg_inds, title=title % 'ecg', labels='ecg')
                            fig.savefig(op.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_ecg_scores.' + picformat),format=picformat) # save figure
                            plt.close(fig)
                            
                            show_picks = np.abs(scores).argsort()[::-1][:5]
                            
                            fig = ica.plot_sources(rawforica, show_picks, exclude=ecg_inds, title=title % 'ecg')
                            fig.savefig(op.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_ecg_sources.' + picformat),format=picformat) # save figure
                            plt.close(fig)
                            
                            fig = ica.plot_components(ecg_inds, title=title % 'ecg', colorbar=True)
                            fig.savefig(op.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_ecg_components.' + picformat),format=picformat) # save figure
                            plt.close(fig)
                            
    #                        ecg_inds = ecg_inds[:n_max_ecg]
                            
                            
                            ica.exclude.extend(ecg_inds)
                        except:
                            from termcolor import colored
                            print(colored('Error finding ECG related components', 'green'))
                            plt.close('all')
                        
                        #EOGV--------------------------------------
                       
    #                    n_max_eog = 1 
                        
                        #eog_epochs = create_eog_epochs(rawforica, ch_name='EOGV', reject=dict(mag=5e-12))  # get single EOG trials
                        
                        try:
                            eog_inds, scores = ica.find_bads_eog(rawforica, ch_name=EOG_channel)
                            
                            fig = ica.plot_scores(scores, exclude=eog_inds, title=title % 'eog', labels='eog')
                            fig.savefig(op.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_eog_scores.' + picformat),format=picformat) # save figure
                            plt.close(fig)
                            
                            show_picks = np.abs(scores).argsort()[::-1][:5]
                            
                            fig = ica.plot_sources(rawforica, show_picks, exclude=eog_inds, title=title % 'eog')
                            fig.savefig(op.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_eog_sources.' + picformat),format=picformat) # save figure
                            plt.close(fig)
                            
                            fig = ica.plot_components(eog_inds, title=title % 'eog', colorbar=True)
                            fig.savefig(op.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_eog_components.' + picformat),format=picformat) # save figure
                            plt.close(fig)
                            
    #                        eog_inds = eog_inds[:n_max_eog]
                            
                            #Double Check---------------------------
                            #eog_evoked = create_eog_epochs(rawforica, ch_name='EOGV', tmin=-.5, tmax=.5, picks=picks_meg).average()
                            #ica.plot_sources(eog_evoked, exclude=eog_inds)  # EOG sources + selection
                            #ica.plot_overlay(eog_evoked, exclude=eog_inds)  # EOG cleaning
                            #---------------------------------------
                            
                            #Exclude EOGV related components
                            ica.exclude.extend(eog_inds)
                        except:
                            from termcolor import colored
                            print(colored('Error finding EOG related components', 'green'))
                            plt.close('all')
                            
                            
                        #Save ICA to path
                        
                        with open(ICA_csv_filename, 'a') as fid:
                            fid.write("{}\t{}\t{}\t{:.3f}\t{:d}\t{:d}\n".format(subj,state,blk,n_components,len(ecg_inds),len(eog_inds)))
                            
                        ica.save(op.join(ICA_path, get_id(subj) + '_ICA_ncomp{}'.format(n_components) + '_blk' + blk + '-ica.fif'))
                        
                        rawforica.close()
                    
                    #EPOCHS 2sec---------------------------------------------------------------------------------------------------------------
                    
                    if save_evokeds:
                        
                        rawforerf = filtered_raw.copy()
                        picks = mne.pick_types(rawforerf.info, meg=True)
                        
                        if ecg_epochs:
                            
                            try:
                                epochs = create_ecg_epochs(rawforerf, tmin=-.5, tmax=.8, ch_name=ECG_channel[0], picks=picks, baseline=baseline)
                            except:
                                epochs = create_ecg_epochs(rawforerf, tmin=-.5, tmax=.8, ch_name=ECG_channel[1], picks=picks, baseline=baseline)
                            
                        else:
                            
                            if event_id1:                    
                                events = mne.find_events(rawforerf)
                            
                            #Otherwise, create fixed length events on which to generate epochs
                            else:
                                event_id = 0                
                                events = make_fixed_length_events(rawforerf, event_id, start = 1, stop = 360, duration = 2)
                    
                            #create epochs
                            try:
                                epochs = mne.Epochs(rawforerf, events=events, event_id=event_id1, tmin=tmin, tmax=tmax, reject = epochs_reject, picks = picks, verbose=True, preload=True, baseline = baseline)
                            except:
                                epochs = mne.Epochs(rawforerf, events=events, event_id=event_id2, tmin=tmin, tmax=tmax, reject = epochs_reject, picks = picks, verbose=True, preload=True, baseline = baseline)
                            #epochs.plot(block = True)
                        
                        rawforerf.close()
                        
                        #average to create evoked
                        #picks = mne.pick_types(epochs.info, meg = True)
                        #evoked = epochs.average(picks = picks) 
                        
                        
    #==============================================================================
    #                     #load ICA file, visualize and apply ICA
    #                     ica_fname = glob.glob(op.join(ICA_path, get_id(subj) + '_ICA_*_blk' + blk + '-ica.fif'))[-1]
    #                     ica = mne.preprocessing.read_ica(ica_fname)
    #                     exclude = ica.exclude                    
    #                     
    #                     
    #                     if remove_ECG_artifacts == False:
    #                         for key in ica.labels_:
    #                             if 'ecg' in key:
    #                                 for ind in ica.labels_[key]:
    #                                     if ind in exclude:
    #                                         exclude.remove(ind)
    #                     
    #                     if remove_EOG_artifacts == False:
    #                         for key in ica.labels_:
    #                             if 'eog' in key:
    #                                 for ind in ica.labels_[key]:
    #                                     if ind in exclude:
    #                                         exclude.remove(ind)
    #                                         
    #                     
    #                     ica.apply(epochs, exclude = exclude)
    #==============================================================================
                        
                        #Plot and manually select bad segments
                        #chan_overview = ['MLC15-2805', 'MLC51-2805', 'MLF14-2805', 'MLF65-2805', 'MLO15-2805', 'MLO51-2805', 'MLP15-2805', 'MLP51-2805', 'MLT14-2805', 'MLT41-2805', 
                        #                 'MRC15-2805', 'MRC51-2805', 'MRF15-2805', 'MRF51-2805', 'MRO15-2805', 'MRO51-2805', 'MRP15-2805', 'MRP51-2805', 'MRT15-2805', 'MRT51-2805',
                        #                 'MZC02-2805', 'MZF02-2805', 'MZO02-2805', 'MZP01-2805',]
                        
                        #picks_viz_epochs = mne.pick_channels(epochs.ch_names, chan_overview)
                        
                        #epochs.plot(picks=picks_viz_epochs)
                        
                        #Drop bad segments
                        epochs.drop_bad()
                        
                        #Save Evokeds to path
                        epochs.save(op.join(epochs_path, get_id(subj) + ('_' + epochs_name if epochs_name else '') + ('_' + state if state else '') + '_epochs' + '_blk' + blk + '-epo.fif'))
                        epochs.average().save(op.join(evoked_path, get_id(subj) + ('_' + epochs_name if epochs_name else '') + ('_' + state if state else '') + '_evoked_{}-{}_blk'.format(l_freq,h_freq) + blk + '-ave.fif'))
                        
                        filtered_raw.close()
                        raw.close
    

#==============================================================================
# #SPECTRAL ANALYSIS---------------------------------------------------------------------------------------------------------
# 
# for subj in subjects:
#     
#     epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs', get_id(subj))
# 
# 
#     blocks = []
#     
#     for key in get_blocks(subj):
#         blocks.append(key)
#     
#     for blk in blocks:
#         
#         epochs = mne.read_epochs(op.join(epochs_path, get_id(subj) + '_epochs' + '_blk' + blk + '-epo.fif'), preload = True)
#         epochs.resample(400., npad='auto', window='boxcar')
#         #plot PSD across all epochs
#         fig = epochs.plot_psd(fmin=2., fmax=40.)
#         fig.savefig(op.join(epochs_path, get_id(subj) + '_epochs' + '_blk' + blk + 'psd.' + picformat),format=picformat) # save figure
#         plt.close(fig)
#         #Spatial distribution of PSD
#         fig = epochs.plot_psd_topomap(vmin=0, vmax=0.03, normalize=True, agg_fun=np.mean)
#         fig.savefig(op.join(epochs_path, get_id(subj) + '_epochs' + '_blk' + blk + 'psd_topomap.' + picformat),format=picformat) # save figure
#         plt.close(fig)
# 
#==============================================================================
