# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:22:13 2017

@author: alex
"""
#IMPORT PACKAGES, DEFINE ENVIRONMENT VARIABLES, CREATE PATHES
#==============================================================================
from header import *
#==============================================================================


def run_ica(subject, task, state, block, raw=None, save=True, n_components=0.975, method='fastica', rejection={'mag':4000e-15}, ECG_channel=['EEG062-2800', 'EEG062'], EOG_channel='EOGV'):
    """
    Fit ICA on raw MEG data.
    Output:
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-ica.fif'
    Log:
        'Analyses/<task>/meg/ICA/ICA_log.tsv'
    Parameters (see mne.preprocessing.ICA):
        raw: raw data to fit ICA on. If None (default), will be loaded according to previous parameters.
        n_components: number of components used for ICA decomposition
        method: the ICA method to use
        rejection: epoch rejection threshold (default to 4000 fT for magnetometers)
        ECG_channel: channel(s) corresponding to the ECG
        EOG_channel: channel(s) corresponding to the EOG
    """
    # Load data
    data_path = op.join(Raw_data_path, task, 'meg')
    raw_fname = op.join(data_path, op.join(* get_rawpath(subject, task=task)) + block + '.ds')
    if not raw:
        raw = mne.io.read_raw_ctf(raw_fname, preload=True)
    
    picks_meg = mne.pick_types(raw.info, meg=True, ref_meg=False, exclude='bads')
    
    # ICA path
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    if not op.exists(ICA_path):
        os.makedirs(ICA_path)
    
    # ICA log
    ICA_log = op.join(Analysis_path, task, 'meg', 'ICA', 'ICA_log.tsv')
    if not op.isfile(ICA_log):
        with open(ICA_log, 'w') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('date','time','subject','state','block','n_components','n_selected_comps','ncomp_ECG','pulse','ncomp_EOG','rejection','dropped_epochs'))
    
    # Filter for ICA
    raw.filter(l_freq=1, h_freq=40, fir_design='firwin', picks=picks_meg, n_jobs=4)
    
    # Fit ICA
    ica = ICA(n_components=n_components, method=method) #create ICA object
    ica.fit(raw, reject=rejection, decim=6, picks=picks_meg) #decimate: 200Hz is more than enough for ICA, saves time
    
    # Detect ECG and EOG artifacts
    ica.scores_ = dict()
    try:
        ica.scores_['ecg'] = ica.find_bads_ecg(raw, ch_name=ECG_channel[0], threshold=0.3)[1] #default threshold at 0.25 too low
        pulse = mne.preprocessing.find_ecg_events(raw, l_freq=8, h_freq=16, ch_name=ECG_channel[0])[2]
    except:
        ica.scores_['ecg'] = ica.find_bads_ecg(raw, ch_name=ECG_channel[1], threshold=0.3)[1]
        pulse = mne.preprocessing.find_ecg_events(raw, l_freq=8, h_freq=16, ch_name=ECG_channel[1])[2]
    ica.exclude.extend(ica.labels_['ecg'])
    
    ica.scores_['eog'] = ica.find_bads_eog(raw, ch_name=EOG_channel)[1]
    ica.exclude.extend(ica.labels_['eog'])
    
    # Plot scores
    ica.plot_scores(ica.scores_['ecg'], exclude=ica.labels_['ecg'], title="ECG artifacts")
    plt.savefig(op.join(ICA_path, '{}_{}-{}_components-scores_ecg.svg'.format(state, block, n_components)))
    plt.close()
    
    ica.plot_scores(ica.scores_['eog'], exclude=ica.labels_['eog'], title="EOG artifacts")
    plt.savefig(op.join(ICA_path, '{}_{}-{}_components-scores_eog.svg'.format(state, block, n_components)))
    plt.close()
    
    # Save ICA
    if save:
        ica.save(op.join(ICA_path, '{}_{}-{}_components-ica.fif'.format(state, block, n_components)))
    
    # Write ICA log
    with open(ICA_log, 'a') as fid:
        fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(time.strftime('%Y_%m_%d\t%H:%M:%S',time.localtime()),subject,state,block,n_components,ica.n_components_,len(ica.labels_['ecg']),int(round(pulse)),len(ica.labels_['eog']),rejection,len(ica.drop_inds_)))

    return ica
    

#subject='010'; task='SMEG'; state='OM'; block='07'; n_components=.975; notch=np.arange(50,301,50); high_pass=0.1; low_pass=None; rejection={'mag':2.5e-12}; epoching={'name':'Cardiac','tmin':-.5,'tmax':.8,'baseline':(-.4,-.3)}; ECG_channel=['EEG062-2800', 'EEG062']; EOG_channel='EOGV'
def process0(subject, task, state, block, raw=None, n_components=.975, ica=None, check_ica=True, notch=np.arange(50,301,50), high_pass=0.1, low_pass=None, rejection={'mag':2.5e-12}, epoching={'name':'Cardiac','tmin':-.5,'tmax':.8,'baseline':(-.4,-.3)}, ECG_channel=['EEG062-2800', 'EEG062'], EOG_channel='EOGV'):
    """
    Run preprocessing to create epochs and evoked response.
    Output:
        'Analyses/<task>/meg/Epochs/<subject>/<epoch_name>-<state>_<block>-epo.fif'
        'Analyses/<task>/meg/Evoked/<subject>/<epoch_name>-<state>_<block>-ave.fif'
    Parameters (see mne.filter):
        raw: raw data to process. If None (default), will be loaded according to previous parameters.
        ica: ICA object. If None (default), will be loaded according to previous parameters.
        low_pass: frequency (in Hz) for low-pass filtering
        high_pass: frequency (in Hz) for high-pass filtering
        notch: frequency (in Hz) or list of frequencies to notch filter
        rejection: epoch rejection threshold (default to 2500 fT for magnetometers)
        epoching: dictionary of epoching parameters. Keys (see mne.Epochs): name (only 'Cardiac' is supported yet), tmin, tmax, baseline
    """
    # Load data
    data_path = op.join(Raw_data_path, task, 'meg')
    raw_fname = op.join(data_path, op.join(* get_rawpath(subject, task=task)) + block + '.ds')
    if not raw:
        raw = mne.io.read_raw_ctf(raw_fname, preload=True)
    
    picks = mne.pick_types(raw.info, meg=True, exclude='bads')
    
    #Load ICA
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    ICA_file = op.join(ICA_path, '{}_{}-{}_components-ica.fif'.format(state, block, n_components))
    if not ica:
        try:
            ica = read_ica(ICA_file)
        except IOError:
            ica = run_ica(subject, task, state, block, raw=raw.copy(), n_components=n_components, ECG_channel=ECG_channel, EOG_channel=EOG_channel)
    
    # Save paths
    epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs', subject)
    if not op.exists(epochs_path):
        os.makedirs(epochs_path)
    epochs_file = op.join(epochs_path, '{}-{}_{}-epo.fif'.format(epoching['name'], state, block))
    
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    if not op.exists(evoked_path):
        os.makedirs(evoked_path)
    evoked_file = op.join(evoked_path, '{}-{}_{}-ave.fif'.format(epoching['name'], state, block))
    
    # Filter
#    raw.notch_filter(notch, picks=picks, n_jobs=4) #python2 has mne 0.14, which does not support 'fir_design' yet
#    raw.filter(l_freq=high_pass, h_freq=low_pass, picks=picks, n_jobs=4, l_trans_bandwidth=0.1)
    raw.notch_filter(notch, fir_design='firwin', picks=picks, n_jobs=4)
    raw.filter(l_freq=high_pass, h_freq=low_pass, fir_design='firwin', picks=picks, n_jobs=4)
    
    if check_ica:
        # Visual check
        ica.plot_overlay(raw)
        plt.savefig(op.join(ICA_path, '{}_{}-{}_components-overlay.svg'.format(state, block, n_components)))
        plt.close()
        
        try:
            check_ecg = create_ecg_epochs(raw, ch_name=ECG_channel[0], picks=picks)
        except:
            check_ecg = create_ecg_epochs(raw, ch_name=ECG_channel[1], picks=picks)
        for comp in ica.labels_['ecg']:
            ica.plot_properties(check_ecg, picks=comp)
            plt.savefig(op.join(ICA_path, '{}_{}-{}_components-properties_ecg{}.svg'.format(state, block, n_components, comp)))
            plt.close()
        
        check_eog = create_eog_epochs(raw, ch_name=EOG_channel, picks=picks)
        for comp in ica.labels_['eog']:
            ica.plot_properties(check_eog, picks=comp)
            plt.savefig(op.join(ICA_path, '{}_{}-{}_components-properties_eog{}.svg'.format(state, block, n_components, comp)))
            plt.close()
    
    # Apply ICA
    ica.apply(raw)
    
    # Epoch
    if epoching['name'] == 'Cardiac':
        try:
            epochs = create_ecg_epochs(raw, reject=rejection, picks=picks, tmin=epoching['tmin'], tmax=epoching['tmax'], baseline=epoching['baseline'], ch_name=ECG_channel[0])
        except:
            epochs = create_ecg_epochs(raw, reject=rejection, picks=picks, tmin=epoching['tmin'], tmax=epoching['tmax'], baseline=epoching['baseline'], ch_name=ECG_channel[1])
    
    # Rejection
#    rejected = epochs.copy()
    epochs.drop_bad()
#    rejected.drop_bad(reject=None, flat=rejection)
    
    # Save epochs
    epochs.save(epochs_file)
    
    # Save evoked
    evoked = epochs.average()
    evoked.save(evoked_file)
    
    drop_log = op.join(Analysis_path, task, 'meg', 'Epochs', 'drop_log.txt')
    with open(drop_log, 'a') as fid:
        fid.write('{} {} epochs dropped (threshold = {})\n'.format(epochs_file.split('/')[-2:], len(np.array(epochs.drop_log)[np.where(epochs.drop_log)]), rejection))


def filtering_parameters(do=True, low_pass=40, high_pass=1, notch=np.arange(50,251,50), resample=200, window='hann'):
    """
    Creates a dictionary of the filtering parameters to be passed to process().
    If do is not True, skip other parameters: no filtering will be applied.
    Parameters (see mne.filter):
        low_pass: frequency (in Hz) for low-pass filtering
        high_pass: frequency (in Hz) for high-pass filtering
        notch: frequency (in Hz) or list of frequencies to notch filter
        resample: frequency (in Hz) to resample to
    """
    filt_param = dict()
    
    if not do:
        return filt_param
    
    filt_param['low_pass'] = low_pass
    filt_param['high_pass'] = high_pass
    filt_param['notch'] = notch
    filt_param['resample'] = resample
    filt_param['window'] = window
    
    return filt_param


def epoching_parameters(do=True, name='Cardiac', function=mne.preprocessing.create_ecg_epochs, tmin=-.5, tmax=.8, baseline=(-0.4, -0.3), reject={'mag':4e-12}):
    """
    Creates a dictionary of the epoching parameters to be passed to process().
    If do is not True, skip other parameters: no epoching will be done.
    Parameters (see mne.Epochs):
        name: name of the epoch (e.g. 'Cardiac' for ECG epochs)
        function: particular function to apply for epoching (e.g. create_ecg_epochs for ECG epochs
        tmin: start time (in s) before event
        tmax: end time (in s) after event
        baseline: time interval (in s) as a tuple (start, end) to apply baseline correction
        reject: rejection parameters
    """
    epoch_param = dict()
    
    if not do:
        return epoch_param
    
    epoch_param['name'] = name
    epoch_param['function'] = function
    epoch_param['tmin'] = tmin
    epoch_param['tmax'] = tmax
    epoch_param['baseline'] = baseline
    epoch_param['reject'] = reject
    
    return epoch_param


def ICA_parameters(do=True, n_components=0.975, method='fastica', ECG_channel=['EEG062-2800', 'EEG062'], EOG_channel='EOGV'):
    """
    Creates a dictionary of the ICA parameters to be passed to process().
    If do is not True, skip other parameters: ICA will not be applied.
    Parameters (see mne.preprocessing.ICA):
        n_components: number of components used for ICA decomposition
        method: the ICA method to use
        ECG_channel: channel(s) corresponding to the ECG
        EOG_channel: channel(s) corresponding to the EOG
    """
    ICA_param = dict()
    
    if not do:
        return ICA_param
    
    ICA_param['n_components'] = n_components
    ICA_param['method'] = method
    ICA_param['ECG_channel'] = ECG_channel
    ICA_param['EOG_channel'] = EOG_channel
    
    return ICA_param


#filter_param = filtering_parameters(do=True, notch=None)
#epoch_param = epoching_parameters(do=True)
#ICA_param = ICA_parameters(do=True)

def process2(task, state, subject, filter_param, epoch_param, ICA_param, picformat='svg'):
    """
    Use the _parameters functions from preprocessing.py to adjust preprocessing parameters.
    Output:
        (ICA directory) *'-ica.fif' along with three images (as picformat files) for checking
        (Evoked directory) *'-ave.fif'
        (Epochs directory) *'-epo.fif'
    """
    
    data_path = op.join(Raw_data_path, task, 'meg')
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs', subject)
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    
    #create the corresponding directories if they don't exist
    if ICA_param and not op.exists(ICA_path):
        os.makedirs(ICA_path)
    
    if epoch_param and not op.exists(epochs_path):
        os.makedirs(epochs_path)
    
    if epoch_param and not op.exists(evoked_path):
        os.makedirs(evoked_path)
    
    
    if ICA_param:
        ICA_csv_filename = op.join(Analysis_path, task, 'meg', 'ICA', 'ICA_info' + time.strftime('_%y-%m-%d_%Hh%M',time.localtime())  + '.csv')
        with open(ICA_csv_filename, 'a') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\n".format('subject','state','block','ncomponents','ncomp_ECG','ncomp_EOG'))
    
    
    blocks = []
    for key in get_blocks(subject, task=task, state=state):
            blocks.append(key)
    
    for blk in blocks:
        
        #file name
        raw_fname = op.join(data_path, op.join(* get_rawpath(subject, task = task)) + blk + '.ds')
        
        #load file
        raw = mne.io.read_raw_ctf(raw_fname, preload=True)
        #print(raw.info)
        
        #load bad channels
        raw.load_bad_channels(op.join(raw_fname, 'BadChannels'))
        
        #define picks
        picks_all = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True, exclude='bads')
        picks_meg = mne.pick_types(raw.info, ref_meg=False, meg=True, eeg=False, eog=False, stim=False, exclude='bads')
        
        
        filtered_raw = raw.copy()
        if filter_param:
            if filter_param['low_pass'] or filter_param['high_pass']:
                filtered_raw.filter(l_freq=filter_param['high_pass'], h_freq=filter_param['low_pass'], picks=picks_all, fir_design='firwin')
                rawforica = filtered_raw.copy()
            
            if filter_param['notch']:
                filtered_raw.notch_filter(filter_param['notch'], picks=picks_all, n_jobs=len(filter_param['notch']))
            
            if filter_param['resample']:
                filtered_raw.resample(filter_param['resample'], window=filter_param['window'])
        
        
        if ICA_param:
            if not rawforica:
                rawforica = raw.copy()
                rawforica.filter(l_freq=1, h_freq=40, l_trans_bandwidth=0.5, h_trans_bandwidth=0.5, filter_length = '10s', phase='zero-double',n_jobs=6, fir_design='firwin')
            
            ica = ICA(n_components=ICA_param['n_components'], method=ICA_param['method'])
            print(ica)
            
            try:
                ica.fit(rawforica, picks=picks_meg, decim=3, reject={mag:6e-12, grad:4000e-13})#mag:2500e-15, no grad
                print(ica)
            except:
                ica.fit(rawforica, picks=picks_meg, decim=3, reject={mag:17e-12, grad:4000e-13})
                print(ica)
            
            
            title = 'Sources related to %s artifacts (red)'
            
            #ECG
            try:
                try:
                    ecg_epochs = create_ecg_epochs(rawforica, tmin=-.5, tmax=.5, ch_name=ECG_channel[0], picks=picks_meg)
                except:
                    ecg_epochs = create_ecg_epochs(rawforica, tmin=-.5, tmax=.5, ch_name=ECG_channel[1], picks=picks_meg)
                
                ecg_inds, scores = ica.find_bads_ecg(ecg_epochs, method='ctps')
                fig = ica.plot_scores(scores, exclude=ecg_inds, title=title % 'ecg', labels='ecg')
                fig.savefig(op.join(ICA_path, subject + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_ecg_scores.' + picformat),format=picformat) # save figure
                plt.close(fig)
                
                show_picks = np.abs(scores).argsort()[::-1][:5]
                
                fig = ica.plot_sources(rawforica, show_picks, exclude=ecg_inds, title=title % 'ecg')
                fig.savefig(op.join(ICA_path, subject + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_ecg_sources.' + picformat),format=picformat) # save figure
                plt.close(fig)
                
                fig = ica.plot_components(ecg_inds, title=title % 'ecg', colorbar=True)
                fig.savefig(op.join(ICA_path, subject + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_ecg_components.' + picformat),format=picformat) # save figure
                plt.close(fig)
                
                
#                ecg_inds = ecg_inds[:3]
                
                ica.exclude.extend(ecg_inds)
                
            except:
                print(colored('Error finding ECG related components', 'green'))
                plt.close('all')
            
            #EOGV
            #eog_epochs = create_eog_epochs(rawforica, ch_name='EOGV', reject=dict(mag=5e-12))  # get single EOG trials
            
            try:
                eog_inds, scores = ica.find_bads_eog(rawforica, ch_name=EOG_channel)
                
                fig = ica.plot_scores(scores, exclude=eog_inds, title=title % 'eog', labels='eog')
                fig.savefig(op.join(ICA_path, subject + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_eog_scores.' + picformat),format=picformat) # save figure
                plt.close(fig)
                
                show_picks = np.abs(scores).argsort()[::-1][:5]
                
                fig = ica.plot_sources(rawforica, show_picks, exclude=eog_inds, title=title % 'eog')
                fig.savefig(op.join(ICA_path, subject + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_eog_sources.' + picformat),format=picformat) # save figure
                plt.close(fig)
                
                fig = ica.plot_components(eog_inds, title=title % 'eog', colorbar=True)
                fig.savefig(op.join(ICA_path, subject + '_ICA_ncomp{}'.format(n_components)+ '_blk' + blk + '_eog_components.' + picformat),format=picformat) # save figure
                plt.close(fig)
                
#                eog_inds = eog_inds[:1]
                
                #Double Check---------------------------
                #eog_evoked = create_eog_epochs(rawforica, ch_name='EOGV', tmin=-.5, tmax=.5, picks=picks_meg).average()
                #ica.plot_sources(eog_evoked, exclude=eog_inds)  # EOG sources + selection
                #ica.plot_overlay(eog_evoked, exclude=eog_inds)  # EOG cleaning
                #---------------------------------------
                
                #Exclude EOGV related components
                ica.exclude.extend(eog_inds)
            except:
                print(colored('Error finding EOG related components', 'green'))
                plt.close('all')
            
            #Save ICA to path
            
            with open(ICA_csv_filename, 'a') as fid:
                fid.write("{}\t{}\t{}\t{:.3f}\t{:d}\t{:d}\n".format(subject,state,blk,n_components,len(ecg_inds),len(eog_inds)))
            
            ica.save(op.join(ICA_path, subject + '_ICA_ncomp{}'.format(n_components) + '_blk' + blk + '-ica.fif'))
            
            rawforica.close()
        
        
        if epoch_param:
            
            rawforerf = filtered_raw.copy()
            picks = mne.pick_types(rawforerf.info, meg=True)
            
            if epoch_param['function']:
                
                try:
                    epochs = epoch_param['function'](rawforerf, tmin=-.5, tmax=.8, ch_name=ECG_channel[0], picks=picks, baseline=baseline)
                except:
                    epochs = epoch_param['function'](rawforerf, tmin=-.5, tmax=.8, ch_name=ECG_channel[1], picks=picks, baseline=baseline)
                
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
    


def process(tasks, states, subjects=None, run_ICA=True, n_components=0.975, ECG_channel=['EEG062-2800', 'EEG062'], EOG_channel='EOGV', save_evokeds = True, ecg_epochs = True):
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
