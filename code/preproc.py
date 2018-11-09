# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:22:13 2017

@author: alex
"""
#IMPORT PACKAGES, DEFINE ENVIRONMENT VARIABLES, CREATE PATHES
#==============================================================================
from header import *
from mne_custom import *
from scipy.signal import detrend, hilbert
#==============================================================================
warnings.filterwarnings("ignore",category=DeprecationWarning)

# MANUAL EXECUTION
#==============================================================================
# method='fastica'; notch=np.arange(50,301,50); high_pass=.5; low_pass=None
# ECG_threshold=0.2; ECG_max=3; EOG_threshold=5; EOG_min=1; EOG_max=2; rejection={'mag':3.5e-12}; ica_rejection={'mag':7e-12}

# subject='050'; state='RS'; block='01'; task='SMEG'; n_components=.975; name='ECG_included'

#raw = mne.io.read_raw_fif(op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}-raw.fif'.format(state, block)), preload=True)
#ica = read_ica(op.join(Analysis_path, task, 'meg', 'ICA', subject, '{}_{}-{}_components-ica.fif'.format(state, block, n_components)))
#ica.exclude = ica.labels_['eog']
#if name != 'ECG_included':
#    ica.exclude += ica.labels_['ecg']
#ica.apply(raw)
#ecg = raw.copy().pick_types(meg=False, ref_meg=False, ecg=True)
#epo = create_ecg_epochs(ecg, reject_by_annotation=False)
#epo.set_channel_types({ch:'eeg' for ch in ecg.ch_names})
#erp = epo.average()
#erp.get_peak(tmin=.15, tmax=.4, mode='neg')

# ica = run_ica(task, subject, state, block, save=False, ECG_threshold=ECG_threshold, EOG_threshold=EOG_threshold, ica_rejection=ica_rejection)
# # ica.labels_['ecg_scores'][np.where(ica.labels_['ecg_scores']>0.1)]

# raw = process(task, subject, state, block, check_ica=False, save_ica=False, ica=None, update_HPI=False)[name]
# # eog = raw[name].copy().pick_types(meg=False, ref_meg=False, eog=True)
# # ecg = raw[name].copy().pick_types(meg=False, ref_meg=False, ecg=True)

# epo = create_ecg_epochs(raw_ECG[name], reject={'mag':3000e-15}, tmin=-.5, tmax=.8, baseline=None, picks=mne.pick_types(raw_ECG.info, meg=True, ecg=True, eog=True, stim=True, exclude='bads'))
# epochs = epoch(task, subject, state, block, save=False, rejection=None, tmin=-.5, tmax=.8, baseline=None, overwrite_ica=False)
#==============================================================================


def load_raw(task, subject, block):
    """
    
    """
    raw = mne.io.read_raw_ctf(get_rawpath(subject, task, block), preload=True, clean_names=True)
    
    raw.info['subject_info'].update({'sub':subject})
    raw.set_channel_types({get_chan_name(subject, 'ecg_chan', raw):'ecg', get_chan_name(subject, 'eogV_chan', raw):'eog', 'UPPT001':'stim'})
    raw.pick_types(meg=True, ecg=True, eog=True, stim=True)
    
    # Crop recording
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0 #Crop recording if first event is less than 2 min after the beginning
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[-1][0]] > 300 else None #Crop recording if last event is more than 5 min after the beginning
    raw.crop(tmin=start, tmax=end)
    
    return raw


def process(task, subject, state, block, notch=np.arange(50,301,50), high_pass=0.5, low_pass=None, EOG_threshold=3, EOG_min=1, EOG_max=None, EOG_score=None, ECG_threshold=.25, ECG_max=3, update_HPI=True, HPI_kwargs=dict(), ICA_kwargs=dict(), custom_args=dict()):
    """
    Run preprocessing and return preprocessed raw data.
    If check_ica, plot overlay and properties of ECG and EOG components (default to True).
    If save_ica, save ICA and ICA plots (deleting previously existing component properties), and write ICA log (default to True).
    If overwrite_ica, artifact scoring will be perfomed again (default to False).
    Output:
        'Analyses/<task>/meg/ICA/fixed_length/<subject>/<state>_<block>-overlay_eog.png'
        'Analyses/<task>/meg/ICA/fixed_length/<subject>/<state>_<block>-properties_eog<i>.png'
    Parameters (see mne.filter):
        ica: ICA object. If None (default), will be loaded according to previous parameters. ICA parameters:
            ica_rejection: epoch rejection threshold (default to 4000 fT for magnetometers)
            ECG_threshold: ECG artifact detection threshold (mne default to 0.25)
            EOG_threshold: EOG artifact detection threshold (mne default to 3.0)
        low_pass: frequency (in Hz) for low-pass filtering (default to None)
        high_pass: frequency (in Hz) for high-pass filtering (default to 0.5)
        notch: frequency (in Hz) or list of frequencies to notch filter (set to np.array([]) for no notch filtering)
        update_HPI:
    """
    # Load data
    raw = mne.io.read_raw_ctf(get_rawpath(subject, task, block), preload=True, clean_names=True)
    
    raw.info['subject_info'].update({'sub':subject})
    raw.set_channel_types({get_chan_name(subject, 'ecg_chan', raw):'ecg', get_chan_name(subject, 'eogV_chan', raw):'eog', 'UPPT001':'stim'})
    raw.pick_types(meg=True, ecg=True, eog=True, stim=True)
    
    # Update head coordinates
    if update_HPI:
        raw = HPI_update(task, subject, block, raw.copy(), **HPI_kwargs)
    
    # Crop recording
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0 #Crop recording if first event is less than 2 min after the beginning
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[-1][0]] > 300 else None #Crop recording if last event is more than 5 min after the beginning
    raw.crop(tmin=start, tmax=end)
    
    ica = raw_ica(task, subject, state, block, raw.copy(), **ICA_kwargs)
    
    # Filter
    if notch.size:
        raw.notch_filter(notch, fir_design='firwin', n_jobs=4)
    raw.filter(l_freq=high_pass, h_freq=low_pass, fir_design='firwin', n_jobs=4)
    
    # Detect EOG and ECG artifacts
    ica.labels_['eog_scores'] = ica.find_bads_eog(raw.copy(), threshold=EOG_threshold)[1].tolist()
    if custom_args:
        ica.labels_['ecg_scores'] = custom_bads_ecg(ica, raw.copy(), custom_args, threshold=ECG_threshold)[1].tolist()
    else:
        ica.labels_['ecg_scores'] = ica.find_bads_ecg(raw.copy(), threshold=ECG_threshold)[1].tolist()
    
    # Fix number of artifactual components
    ica.labels_['ecg'] = ica.labels_['ecg'][:ECG_max]
    ica.labels_['eog'] = ica.labels_['eog'][:EOG_max]
    if EOG_score:
        ica.labels_['eog'] = sorted(np.where(np.abs(ica.labels_['eog_scores']) >= EOG_score)[0])
    if EOG_min and not ica.labels_['eog']:
        ica.labels_['eog'] = np.argsort(np.abs(ica.labels_['eog_scores']))[::-1].tolist()[:EOG_min]
    
    # AutoReject
    raw_clean = raw.copy()
    ica.apply(raw_clean, exclude=ica.labels_['eog']+ica.labels_['ecg'])
    raw_clean, threshold = auto_annotate(raw_clean)
    raw.annotations = raw_clean.annotations
    
    # Save artefact detection
    ica.labels_['ECG_threshold'] = ECG_threshold
    ica.labels_['EOG_threshold'] = EOG_score if EOG_score else EOG_threshold
    ica.labels_['raw_rejection'] = threshold
    ica.save(ica.labels_['filename'])
    
    # Save pre-processed data
    raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}_{}-raw.fif'.format(subject, state, block))
    os.makedirs(op.dirname(raw_file), exist_ok=True)
    raw.save(raw_file, overwrite=True)
    
    return raw_clean, raw, ica


# # /!\ Custom attributes (e.g., ica.scores_) are not kept upon .save(), which calls _write_ica() whose dict ica_misc is not editable on call.
# # => exploit the attribute labels_
# # # /!\ Numpy arrays are not supported --> convert to type list with .tolist()

def raw_ica(task, subject, state, block, raw=None, save=True, overwrite_fit=False, n_components=None, method='fastica', ica_rejection='auto', max_iter=300):
    """
    Fit ICA on raw MEG data and return ICA object.
    If save, save ICA, save ECG and EOG artifact scores plots, and write log (default to True).
    If fit_ica, fit ICA even if there is already an ICA file (default to False).
    
    Output
    ------
    'Analyses/<task>/meg/Raw/ICA/<subject>/<state>_<block>-ica.fif'
    
    'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-scores_eog.png'
    
    Log
    ---
    'Analyses/<task>/meg/Raw/ICA/ICA_log.tsv'
    
    Parameters (see mne.preprocessing.ICA)
    ----------
    raw : raw data to fit ICA on
        If None (default), will be loaded according to previous parameters.
    n_components : 
        number of components used for ICA decomposition
    method : 
        the ICA method to use
    ica_rejection : 
        epoch rejection threshold (default to 4000 fT for magnetometers)
    EOG_threshold : 
        EOG artifact detection threshold (mne default to 3.0)
    EOG_min : 
        minimum number of EOG components to exclude (default to 1)
    EOG_max : 
        maximum number of EOG components to exclude (default to 2)
    
    """
    # ICA path
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    if save and not op.isdir(ICA_path):
        os.makedirs(ICA_path)
    ICA_file = op.join(ICA_path, '{}_{}_{}-{}_raw-ica.fif'.format(subject, state, block, method))
    
    if op.isfile(ICA_file) and not overwrite_fit:
        ica = read_ica(ICA_file)
        return ica
    
    ica = ICA(n_components=n_components, method=method, max_iter=max_iter) #create ICA object
    ica.drop_inds_ = []
    ica.labels_ = dict()
    
    # Filter for ICA
    if not raw:
        raw = load_raw(task, subject, block)
    raw.filter(l_freq=1, h_freq=40, fir_design='firwin', n_jobs=4)
    
    if ica_rejection is 'auto':
        auto_epochs = mne.Epochs(raw.copy(), make_fixed_length_events(raw.copy(), 111, duration=2, first_samp=False), tmin=0, tmax=2, baseline=None, picks=mne.pick_types(raw.info), reject=None)
        ica_rejection = get_rejection_threshold(auto_epochs)
    
    # Fit ICA
    ica.fit(raw, reject=ica_rejection, decim=6, picks=mne.pick_types(raw.info, meg=True, ref_meg=False)) #decimate: 200Hz is more than enough for ICA, saves time; picks: fit only on MEG
    ica.labels_['rejection'] = ica_rejection
    ica.labels_['drop_inds_'] = ica.drop_inds_
    ica.labels_['filename'] = ICA_file
    
    # Save ICA
    if save:
        ica.save(ICA_file)
    
    return ica


def check_preproc(task, subject, state, block, raw=None, ica=None, report=None, save_report=True, custom_args=dict()):
    """
    
    """
    # Initialisation
    if not raw:
        raw = load_preproc(task, subject, state, block, exclude_eog=False, exclude_ecg=False)
    if not ica:
        ica = raw_ica(task, subject, state, block)
    if not report:
        report = Report(subject=subject, title='{} {} {} - Preprocessing report'.format(subject, state, block), image_format='svg')
    
    figs = dict()
    
    # EOG plots
    check_eog = create_eog_epochs(raw.copy(), reject=ica.labels_['rejection'])
    
    figs['{} EOG scores'.format(state+block)] = ica.plot_scores(ica.labels_['eog_scores'], exclude=ica.labels_['eog'], labels='eog', axhline=ica.labels_['EOG_threshold'] if ica.labels_['EOG_threshold'] <= 1 else None, figsize=(8,2), show=False)
    for comp in ica.labels_['eog']:
        figs['{} EOG {}'.format(state+block, comp)] = ica.plot_properties(check_eog, picks=comp, show=False)
    figs['{} EOG overlay'.format(state+block)] = ica.plot_overlay(check_eog.average(), exclude=ica.labels_['eog'], show=False)
    
    # ECG Plots
    if custom_args:
        check_ecg, pulse = custom_ecg_epochs(raw.copy(), custom_args, reject=ica.labels_['rejection'])
    else:
        check_ecg = create_ecg_epochs(raw.copy(), reject=ica.labels_['rejection'])
    
    figs['{} ECG scores'.format(state+block)] = ica.plot_scores(ica.labels_['ecg_scores'], exclude=ica.labels_['ecg'], labels='ecg', axhline=ica.labels_['ECG_threshold'], figsize=(8,2), show=False)
    for comp in ica.labels_['ecg']:
        figs['{} ECG {}'.format(state+block, comp)] = ica.plot_properties(check_ecg, picks=comp, show=False)
    figs['{} ECG overlay'.format(state+block)] = ica.plot_overlay(check_ecg.average(), exclude=ica.labels_['ecg'], show=False)
    
    # Save report
    t_start = time.perf_counter()
    report.add_htmls_to_section('{}s of data rejected (threshold = {:.0f} fT).'.format(rejected_duration(raw, ica.labels_['raw_rejection'])), 'Total rejected duration', state+block)
    report.add_figs_to_section(list(figs.values()), list(figs.keys()), state+block)
    logger.info("Editing Report took {:.1f}s.".format(time.perf_counter() - t_start))
    if save_report:
        report_file = op.join(Analysis_path, task, 'meg', 'Reports', subject, '{}_Preprocessing-report.html'.format(subject))
        os.makedirs(op.dirname(report_file), exist_ok=True)
        report.save(report_file, open_browser=False, overwrite=True)
    
    # ICA log
    ICA_log = op.join(Analysis_path, task, 'meg', 'ICA', 'ICA_log.tsv')
    if not op.isfile(ICA_log):
        with open(ICA_log, 'w') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('subject','state','block','start','end','n_selected_comps','ncomp_EOG','ncomp_ECG','rejection /fT','dropped_epochs'))
    with open(ICA_log, 'a') as fid:
        fid.write("{}\t{}\t{}\t{:.3f}\t{:.3f}\t{}\t{}\t{}\t{.0f}\t{}\n".format(subject,state,block,raw.first_samp/raw.info['sfreq'],raw.times[-1]+raw.first_samp/raw.info['sfreq'],ica.n_components_,len(ica.labels_['eog']),len(ica.labels_['ecg']),ica.labels_['rejection']*1e15,len(ica.labels_['drop_inds_'])))
        
    return report
    

def process0(task, subject, state, block, n_components=.975, ica=None, check_ica=True, save_ica=True, overwrite_ica=False, fit_ica=False, ica_rejection={'mag':4000e-15}, notch=np.arange(50,301,50), high_pass=0.5, low_pass=None, ECG_threshold=0.25, EOG_threshold=3, custom_args=dict(), update_HPI=True, precision='0.5cm', opt='start'):
    """
    Run preprocessing and return preprocessed raw data.
    If check_ica, plot overlay and properties of ECG and EOG components (default to True).
    If save_ica, save ICA and ICA plots (deleting previously existing component properties), and write ICA log (default to True).
    If overwrite_ica, artifact scoring will be perfomed again (default to False).
    Output:
        'Analyses/<task>/meg/ICA/fixed_length/<subject>/<state>_<block>-overlay_eog.png'
        'Analyses/<task>/meg/ICA/fixed_length/<subject>/<state>_<block>-properties_eog<i>.png'
    Parameters (see mne.filter):
        ica: ICA object. If None (default), will be loaded according to previous parameters. ICA parameters:
            ica_rejection: epoch rejection threshold (default to 4000 fT for magnetometers)
            ECG_threshold: ECG artifact detection threshold (mne default to 0.25)
            ECG_threshold: ECG artifact detection threshold (mne default to 3.0)
        low_pass: frequency (in Hz) for low-pass filtering (default to None)
        high_pass: frequency (in Hz) for high-pass filtering (default to 0.5)
        notch: frequency (in Hz) or list of frequencies to notch filter (set to np.array([]) for no notch filtering)
        update_HPI:
    """
    # Load data
    data_path = op.join(Raw_data_path, task, 'meg')
    raw_fname = op.join(data_path, op.join(* get_rawpath(subject, task=task)) + block + '.ds')
    raw = mne.io.read_raw_ctf(raw_fname, preload=True)
    
    # Load ICA
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    ICA_file = op.join(ICA_path, '{}_{}_{}-raw-ica.fif'.format(subject, state, block))
    if not ica:
        if overwrite_ica or not op.isfile(ICA_file):
#            ica = run_ica(task, subject, state, block, raw=raw.copy(), n_components=n_components, save=save_ica, fit_ica=fit_ica, ica_rejection=ica_rejection, ECG_threshold=ECG_threshold, EOG_threshold=EOG_threshold, custom_args=custom_args)
            ica = raw_ica(task, subject, state, block, raw=raw.copy(), n_components=n_components, save=save_ica, fit_ica=fit_ica, ica_rejection=ica_rejection, EOG_threshold=EOG_threshold)
        else:
            ica = read_ica(ICA_file)
    
    raw.info['subject_info'].update({'sub':subject})
    raw.set_channel_types({get_chan_name(subject, 'ecg_chan', raw):'ecg', get_chan_name(subject, 'eogV_chan', raw):'eog', 'UPPT001':'stim'})
    raw.pick_types(meg=True, ecg=True, eog=True, stim=True, exclude='bads')
    
    # Update head coordinates
    if update_HPI:
        raw = HPI_update(task, subject, block, raw.copy(), precision=precision, opt=opt, reject_head_mvt=True)
    
    # Crop recording
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[-1][0]] > 300 else None
    raw.crop(tmin=start, tmax=end)
    
    # Filter
    if notch.size:
        raw.notch_filter(notch, fir_design='firwin', n_jobs=4)
    raw.filter(l_freq=high_pass, h_freq=low_pass, fir_design='firwin', n_jobs=4)
    
    # Visual check
    if check_ica:
        plot_path = op.join(Analysis_path, task, 'meg', 'Plots', 'Preprocessing', subject)
        os.makedirs(plot_path, exist_ok=True)
        
        # EOG components
        check_eog = create_eog_epochs(raw)
        ica.plot_overlay(check_eog.average(), exclude=ica.labels_['eog'])
        if save_ica:
            plt.savefig(op.join(plot_path, '{}_{}_{}-overlay_eog.png'.format(subject, state, block)))#, transparent=True, dpi=360)
            plt.close()
            for pic in glob.glob(op.join(plot_path, state+'*'+block+'*properties_eog*')):
                os.remove(pic)
        
        for comp in ica.labels_['eog']:
            ica.plot_properties(check_eog, picks=comp)
            if save_ica:
                plt.savefig(op.join(plot_path, '{}_{}_{}-properties_eog{}.png'.format(subject, state, block, comp)))#, transparent=True)
                plt.close()
    
    # Apply ICA
    ica.exclude = ica.labels_['eog']
    ica.apply(raw)
    
    # Save pre-processed data
    raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}_{}-raw.fif'.format(subject, state, block))
    os.makedirs(op.dirname(raw_file), exist_ok=True)
    raw.save(raw_file, overwrite=True)
    
    return raw


# # /!\ Custom attributes (e.g., ica.scores_) are not kept upon .save(), which calls _write_ica() whose dict ica_misc is not editable on call.
# # => exploit the attribute labels_
# # # /!\ Numpy arrays are not supported --> convert to type list with .tolist()

def raw_ica0(task, subject, state, block, raw=None, save=True, fit_ica=False, n_components=0.975, method='fastica', ica_rejection={'mag':4000e-15}, EOG_threshold=3, EOG_min=1, EOG_max=2, custom_args=dict()):
    """
    Fit ICA on raw MEG data and return ICA object.
    If save, save ICA, save ECG and EOG artifact scores plots, and write log (default to True).
    If fit_ica, fit ICA even if there is already an ICA file (default to False).
    
    Output
    ------
    'Analyses/<task>/meg/Raw/ICA/<subject>/<state>_<block>-ica.fif'
    
    'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-scores_eog.png'
    
    Log
    ---
    'Analyses/<task>/meg/Raw/ICA/ICA_log.tsv'
    
    Parameters (see mne.preprocessing.ICA)
    ----------
    raw : raw data to fit ICA on
        If None (default), will be loaded according to previous parameters.
    n_components : 
        number of components used for ICA decomposition
    method : 
        the ICA method to use
    ica_rejection : 
        epoch rejection threshold (default to 4000 fT for magnetometers)
    EOG_threshold : 
        EOG artifact detection threshold (mne default to 3.0)
    EOG_min : 
        minimum number of EOG components to exclude (default to 1)
    EOG_max : 
        maximum number of EOG components to exclude (default to 2)
    
    """
    # Load data
    data_path = op.join(Raw_data_path, task, 'meg')
    raw_fname = get_rawpath(subject, task, block)
    if not raw:
        raw = mne.io.read_raw_ctf(raw_fname, preload=True)
    
    # ICA path
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    if save and not op.isdir(ICA_path):
        os.makedirs(ICA_path)
    ICA_file = op.join(ICA_path, '{}_{}_{}-raw-ica.fif'.format(subject, state, block))
    
    # ICA log
    ICA_log = op.join(Analysis_path, task, 'meg', 'ICA', 'ICA_log.tsv')
    if save and not op.isfile(ICA_log):
        with open(ICA_log, 'w') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('date','time','subject','state','block','start','end','n_components','n_selected_comps','ncomp_EOG','rejection','dropped_epochs'))
    
    raw.set_channel_types({get_chan_name(subject, 'ecg_chan', raw):'ecg', get_chan_name(subject, 'eogV_chan', raw):'eog', 'UPPT001':'stim'})
    
    # Crop recording
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[1][0]] > 300 else None
    raw.crop(tmin=start, tmax=end)
    
    raw.pick_types(meg=True, ref_meg=False, ecg=True, eog=True, exclude='bads')
    
    # Filter for ICA
    raw.filter(l_freq=1, h_freq=40, fir_design='firwin', n_jobs=4)
    
    if ica_rejection is 'auto':
        auto_epochs = mne.Epochs(raw.copy(), make_fixed_length_events(raw.copy(), 111, duration=2, first_samp=False), tmin=0, tmax=2, baseline=None, picks=mne.pick_types(raw.info), reject=None)
        ica_rejection = get_rejection_threshold(auto_epochs)
    
    # Fit ICA
    if fit_ica or not op.isfile(ICA_file):
        ica = ICA(n_components=n_components, method=method) #create ICA object
        ica.exclude = []
        ica.drop_inds_ = []
        ica.labels_ = dict()
        
        ica.fit(raw, reject=ica_rejection, decim=6, picks=mne.pick_types(raw.info, meg=True)) #decimate: 200Hz is more than enough for ICA, saves time; picks: fit only on MEG
        ica.labels_['rejection'] = ica_rejection
        ica.labels_['drop_inds_'] = ica.drop_inds_
    
    else:
        ica = read_ica(ICA_file)
    
    # Detect EOG artifacts
    ica.labels_['eog_scores'] = ica.find_bads_eog(raw.copy(), threshold=EOG_threshold)[1].tolist()
    
    # Fix number of artifactual components
    ica.labels_['eog'] = ica.labels_['eog'][:EOG_max]
    if EOG_min and not ica.labels_['eog']:
        ica.labels_['eog'] = np.argsort(np.abs(ica.labels_['eog_scores'])).tolist()
        ica.labels_['eog'] = ica.labels_['eog'][::-1][:EOG_min]
    
    # Tag for exclusion
    ica.exclude = ica.labels_['eog']
    
    # Plot scores
    ica.plot_scores(ica.labels_['eog_scores'], exclude=ica.labels_['eog'], labels='eog')
    if save:
        plot_path = op.join(Analysis_path, task, 'meg', 'Plots', 'Preprocessing', subject)
        os.makedirs(plot_path, exist_ok=True)
        plt.savefig(op.join(plot_path, '{}_{}_{}-scores_eog.png'.format(subject, state, block)))#, transparent=True)
        plt.close()
    
    # Save ICA
    if save:
        ica.save(ICA_file)
        # Write ICA log
        with open(ICA_log, 'a') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(time.strftime('%Y_%m_%d\t%H:%M:%S',time.localtime()),subject,state,block,int(round(start)),int(round(end)) if end else int(round(raw.times[-1])),n_components,ica.n_components_,len(ica.labels_['eog']),ica.labels_['rejection'],len(ica.labels_['drop_inds_'])))
    
    return ica


def HPI_update(task, subject, block, data, precision='0.5cm', opt='start', reject_head_mvt=True):
    """
    Load *-epo.fif Epochs from Analysis/<task>/meg/Epochs/<subject>/.
    Coordinates are updated according to .hc files in Analysis/<task>/meg/Coregistration/<subject>/.
    Corresponding bad segments due to head movement are rejected.
    The average *-ave.fif Evoked response is then saved in Analysis/<task>/meg/Evoked/<subject>/.
    Parameters:
        opt: argument to choose between multiple .hc files
    """
    
    #==============================================================================
    hc_path = op.join(Analysis_path, task, 'meg', 'Coregistration', subject)
    #==============================================================================
    
    # Find the .hc file corresponding to this block
    hc_files = glob.glob(op.join(hc_path, '*{}*_{}*.hc'.format(precision,block)))
    hc_file = ''
    if len(hc_files) > 1:
        for hc in hc_files:
            if not opt and not '+' in hc:
                hc_file = hc
                continue
            if opt and opt in hc:
                hc_file = hc
                continue
        if not hc_file:
            hc_file = hc_files[0]
            warnings.warn("Multiple .hc files found, none matching the desired option. Using {}.".format(hc_file))
    else:
        hc_file = hc_files[0]
    
    bs_fname = hc_file[:-3] + '-bad.segments' + block
    bads = False
    if os.stat(bs_fname).st_size > 0:
        bad_segments = pd.read_table(bs_fname, header=None)
        bads = True
    
    # Read and load the content of the .hc file    
#==============================================================================
    with open(hc_file) as f:
        content = f.readlines()

    content = [x.strip() for x in content] 
    
    Nasion_ctf = [0,0,0,1]
    Left_ctf = [0,0,0,1]
    Right_ctf = [0,0,0,1]
    
    Nasion_ctf[0] = float(content[25].strip('x = '))
    Nasion_ctf[1] = float(content[26].strip('y = '))
    Nasion_ctf[2] = float(content[27].strip('z = '))
    
    Left_ctf[0] = float(content[29].strip('x = '))
    Left_ctf[1] = float(content[30].strip('y = '))
    Left_ctf[2] = float(content[31].strip('z = '))
    
    Right_ctf[0] = float(content[33].strip('x = '))
    Right_ctf[1] = float(content[34].strip('y = '))
    Right_ctf[2] = float(content[35].strip('z = '))
#==============================================================================
    
    
    # CTF -> Neuromag Coordinate transformation matrix
#==============================================================================
    ctf_head_trans = data.info['ctf_head_t']['trans']
    
    # New Neuromag coordinates
    
    Nasion_head = np.dot(ctf_head_trans, Nasion_ctf)*0.01
    
    Left_head = np.dot(ctf_head_trans, Left_ctf)*0.01
    
    Right_head = np.dot(ctf_head_trans, Right_ctf)*0.01
#==============================================================================
    
    
    
    # DEV <-> CTF Coordinate transformation matrix
#==============================================================================
    dev_ctf_trans = data.info['dev_ctf_t']['trans']
    ctf_dev_trans = inv(dev_ctf_trans)
    
    # New Dev coordinates                
    
    Nasion_dev = np.dot(ctf_dev_trans, Nasion_ctf)*0.01
    
    Left_dev = np.dot(ctf_dev_trans, Left_ctf)*0.01
    
    Right_dev = np.dot(ctf_dev_trans, Right_ctf)*0.01
#==============================================================================
    
    
    
    # UPDATE THE INFO 'hpi_results' DATA WITH NEW DEV COORDINATES
#==============================================================================
    for i in range(len(data.info['hpi_results'][0]['dig_points'])):
        
        #Just checking it's the right kind of digitizer ...
        if data.info['hpi_results'][0]['dig_points'][i]['kind'] == 1:
            
            #Left, Nasion, Right
            if data.info['hpi_results'][0]['dig_points'][i]['ident'] == 1:
                
                data.info['hpi_results'][0]['dig_points'][i]['r'] = Left_dev[:3]
                #data.info['hpi_results'][0]['dig_points'][i]['coord_frame'] = 1
                
            if data.info['hpi_results'][0]['dig_points'][i]['ident'] == 2:
                
                data.info['hpi_results'][0]['dig_points'][i]['r'] = Nasion_dev[:3]
                #data.info['hpi_results'][0]['dig_points'][i]['coord_frame'] = 1
                
            if data.info['hpi_results'][0]['dig_points'][i]['ident'] == 3:
                
                data.info['hpi_results'][0]['dig_points'][i]['r'] = Right_dev[:3]
                #data.info['hpi_results'][0]['dig_points'][i]['coord_frame'] = 1
#==============================================================================
                
                
    
    # UPDATE THE INFO 'dig' DATA WITH NEUROMAG HEAD COORDINATES
#==============================================================================
    for i in range(len(data.info['dig'])):
        
        #Just checking it's the right kind of digitizer ...
        if data.info['dig'][i]['kind'] == 1:
            
            #Left, Nasion, Right
            if data.info['dig'][i]['ident'] == 1:
                
                data.info['dig'][i]['r'] = Left_head[:3]
                data.info['dig'][i]['coord_frame'] = 4
                
            if data.info['dig'][i]['ident'] == 2:
                
                data.info['dig'][i]['r'] = Nasion_head[:3]
                data.info['dig'][i]['coord_frame'] = 4
                
            if data.info['dig'][i]['ident'] == 3:
                
                data.info['dig'][i]['r'] = Right_head[:3]
                data.info['dig'][i]['coord_frame'] = 4
#==============================================================================
    
    if bads and reject_head_mvt:
#        for i in range(len(bad_segments)):
        new_annot = mne.Annotations(bad_segments.get_values()[:,0], np.diff(bad_segments.get_values()).ravel(), ['bad head mvt'] * bad_segments.get_values().shape[0])
        if data.annotations:
            for onset, duration, description in zip(new_annot.onset, new_annot.duration, new_annot.description):
                data.annotations.append(onset, duration, description)
        else:
            data.annotations = new_annot
#            events = data.events[:,0]/data.info['sfreq'] #extract event timing
#            data.drop(np.logical_and(bad_segments.iloc[i][0] < events + data.times[-1], events + data.times[0] < bad_segments.iloc[i][1]), reason = 'head_movement')
#            #reject the epoch if it ends after the beginning of the bad segment, and starts before the end of the bad segment
    
    return data


def R_T_ECG_events(task, subject, state, block, raw=None, custom_args=dict(), R_id=999, T_id=333, T_window=[.15,.4], var=.05, l_freq=None, h_freq=None, save=True, reject_gaussian=True, coverage=2):
    """
    From raw data containing at least the ECG channel, returns events and event_id corresponding to both R and T peaks. If save=True, saves their timing in 'Analyses/<task>/meg/Epochs/T_timing.tsv'.
    """
    # Save T peak timing
    timing_file = op.join(Analysis_path, 'MEG', 'meta', 'T_timing-{}_{}-{}_{}-{}_var-{}_coverage.tsv'.format(l_freq, h_freq, T_window[0], T_window[1], var, coverage))
    if save and not op.isfile(timing_file):
        with open(timing_file, 'w') as fid:
            fid.write('subject\tstate\tblock\tR_peak\tT_delay\tT_included\n')
    
    # Load data
    if not raw:
        raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}_{}-raw.fif'.format(subject, state, block))
        raw = mne.io.read_raw_fif(raw_file, preload=True)
    
    # Epoch ECG on R peak
    ecg = raw.copy().pick_types(meg=False, ref_meg=False, ecg=True)
    if custom_args:
        epo, pulse = custom_ecg_epochs(ecg, custom_args, reject_by_annotation=False, event_id=R_id)
    else:
        epo = create_ecg_epochs(ecg, reject_by_annotation=False, event_id=R_id)
    
    # Save R times
    R_times_file = op.join(Analysis_path, 'MEG', 'meta', 'ECG_R_times.tsv')
    if not op.isfile(timing_file):
        with open(timing_file, 'w') as fid:
            fid.write('subject\tstate\tblock\tTime\n')
    with open(timing_file, 'a') as fid:
        for t in epo.times[:,0]:
            fid.write("{}\t{}\t{}\t{:.3f}\n".format(subject, state, block, raw.times[t]))
    np.savetxt(op.join(Analysis_path, 'MEG', 'meta', 'ECG_R_times.tsv'), raw.times[epo.events[:,0]], delimiter='\t', header='Timing /s')
    
    # Get the average ECG waveform
    epo.set_channel_types({ch:'eeg' for ch in ecg.ch_names})
    erp = epo.average()
    
    # Find R / T peak sign
    if custom_args:
        R_sign = custom_args['T_sign'] if 'T_sign' in custom_args.keys() else custom_args['R_sign']
    else:
        R_sign = np.sign(np.subtract(erp.data.ravel(), np.median(erp.data.ravel()))[erp.time_as_index(0)][0])
    if not R_sign:
        raise ValueError("Sign of the R peak could not be determined.")
    
    # Get average T peak timing
    pos_data = erp.data * R_sign
    pos_data -= np.min(pos_data)
    erp_pos = mne.EvokedArray(pos_data, erp.info, tmin=erp.times[0], nave=erp.nave, comment='Positive ECG R peak')
    ecg_channel, T_peak = erp_pos.get_peak(tmin=T_window[0], tmax=T_window[1], mode='pos')
    
    # Find event-specific T peak timing
    T_window_i = epo.time_as_index([T_peak - var, T_peak + var])
    data = epo.get_data()[:, 0, :]
    T_times_i = (R_sign*data[:, T_window_i[0]:T_window_i[1]]).argmax(axis=1) + T_window_i[0]
    
    mask = np.where(T_times_i)
    if reject_gaussian:
        mu = np.mean(T_times_i)
        sigma = np.std(T_times_i)
        filt_win = [mu - coverage * sigma, mu + coverage * sigma]
        
        mask = np.where(np.logical_and(filt_win[0] < T_times_i, T_times_i < filt_win[1]))
    
    # Create T events
    T_events = epo.copy().events
    T_events[:,0] += T_times_i - epo.time_as_index(0)
    T_events[:,2] = T_id
    T_events = T_events[mask]
    
    # Save timing of R peaks and delay until T peak (in sec)
    T_times = epo.times[T_times_i]
    included = np.zeros(T_times.shape)
    included[mask] = 1
    if save:
        for i in range(len(T_times)):
            with open(timing_file, 'a') as fid:
                fid.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(subject, state, block, np.round(raw.times[epo.events[i,0]],3), np.round(T_times[i], 3), included[i]))
    
    # Combine R and T events
    events = np.concatenate([epo.events, T_events])
    events = events[events[:,0].argsort()]
    event_id = {'R': R_id, 'T': T_id}
    
    # Save events
    event_file = op.join(Analysis_path, task, 'meg', 'Epochs', 'Events', subject, '{}_{}_{}+R_{}+T_{}.eve'.format(subject, state, block, R_id, T_id))
    os.makedirs(op.dirname(event_file), exist_ok=True)
    mne.write_events(event_file, events)
    
    return events, event_id, erp_pos


def ECG_ICA(task, subject, state, block, raw=None, events=np.array([]), event_id=None, rejection={'mag': 7000e-15}, save=True, fit_ica=False, n_components=0.975, method='picard', Rwin=-.05, Twin=0, ECG_threshold=0.25, ECG_max=3):
    """
    
    """
    # Load data
    if not raw:
        raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}_{}-raw.fif'.format(subject, state, block))
        raw = mne.io.read_raw_fif(raw_file, preload=True)
    
    # ICA path
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    os.makedirs(ICA_path, exist_ok=True)
    ICA_file = op.join(ICA_path, '{}_{}_{}-ECG-ica.fif'.format(subject, state, block))
    
    rawforica = raw.copy()
    rawforica.filter(l_freq=1, h_freq=40, fir_design='firwin', n_jobs=4)
    if not events.size and not event_id:
        #Load events and extract event ids
        event_file = glob.glob(op.join(Analysis_path, task, 'meg', 'Epochs', 'Events', subject, '{}_{}_{}+R_*+T_*.eve'.format(subject, state, block)))[0]
        events = mne.read_events(event_file)
        
        ids = op.splitext(op.basename(event_file))[0].split('+')[1:]
        event_id = dict()
        for i in ids:
            k, v = i.split('_')
            event_id[k] = int(v)
    
    if rejection is 'auto':
        auto_epochs = mne.Epochs(rawforica.copy(), make_fixed_length_events(rawforica.copy(), 111, duration=2, first_samp=False), tmin=0, tmax=2, baseline=None, picks=mne.pick_types(rawforica.info), reject=None)
        rejection = get_rejection_threshold(auto_epochs)
    
    epoforica = mne.Epochs(rawforica, events, event_id, baseline=None, reject=rejection, reject_by_annotation=True, preload=True)
    epoforica.equalize_event_counts(list(event_id.keys()))
    RT_delay_i = int(np.round(np.median(epoforica['T'].events[:,0] - epoforica['R'].events[:,0])))
    RT_delay = epoforica.times[epoforica.time_as_index(0)[0] + RT_delay_i]
    
    epoforica = epoforica['R']
    epoforica.crop(Rwin, RT_delay + Twin)
    
    if fit_ica or not op.isfile(ICA_file):
        ica = ICA(n_components=n_components, method=method) #create ICA object
        ica.exclude = []
        ica.labels_ = dict()
        
        ica.fit(epoforica, picks=mne.pick_types(epoforica.info, ref_meg=False), decim=6) #decimate: 200Hz is more than enough for ICA, saves time; picks: fit only on MEG
    
    else:
        ica = read_ica(ICA_file)
    
    # Detect ECG artifacts
    ica.labels_['ecg_scores'] = ica.find_bads_ecg(epoforica, threshold=ECG_threshold)[1].tolist()
    
    # Fix number of artifactual components
    ica.labels_['ecg'] = ica.labels_['ecg'][:ECG_max]
    
    # Tag for exclusion
    ica.exclude = ica.labels_['ecg']
    
    # Plot scores
    plot_path = op.join(Analysis_path, task, 'meg', 'Plots', 'Preprocessing', subject)
    os.makedirs(plot_path, exist_ok=True)
    
    ica.plot_scores(ica.labels_['ecg_scores'], exclude=ica.labels_['ecg'], labels='ecg', axhline=ECG_threshold)
    if save:
        plt.savefig(op.join(plot_path, '{}_{}_{}-scores_ecg.png'.format(subject, state, block)))#, transparent=True)
        plt.close()
    
    epochs = mne.Epochs(raw, events, event_id, baseline=None, reject_by_annotation=True, preload=True)['R']
    
    ica.plot_overlay(epochs.average())
    if save:
        plt.savefig(op.join(plot_path, '{}_{}_{}-overlay_ecg.png'.format(subject, state, block)))#, transparent=True, dpi=360)
        plt.close()
        for pic in glob.glob(op.join(plot_path, '*'+state+'*'+block+'*properties_ecg*')):
            os.remove(pic)
    
    for comp in ica.labels_['ecg']:
        ica.plot_properties(epochs, picks=comp)
        if save:
            plt.savefig(op.join(plot_path, '{}_{}_{}-properties_ecg{}.png'.format(subject, state, block, comp)))#, transparent=True)
            plt.close()
    
    # Save ICA
    if save:
        ica.save(ICA_file)
    
    return ica


def check_ecg(task, subject, state, block, ecg_erp, raw, events, report=None, save_report=True):
    """
    Plot ECG along with events used for epoching. If synthetic, create and plot a synthetic ECG channel as well.
    """
    if not report:
        report = Report(subject=subject, title='{} {} {} - Preprocessing report'.format(subject, state, block), image_format='svg')
    figs = dict()
    
    figs['{} ECG waveform'.format(state+block)] = ecg_erp.plot()
    
    # Load data
    if not raw:
        raw = load_preproc(task, subject, state, block, exclude_eog=True, exclude_ecg=False)
    
    ecg = raw.copy().pick_types(meg=False, ref_meg=False, ecg=True)
    
    # Add synthetic channel as in MNE's create_ecg_epochs():
    ecg_syn = raw.get_data(picks=mne.pick_types(raw.info, meg='mag', ref_meg=False)).mean(axis=0)
    ecg_raw = mne.io.RawArray(ecg_syn[None], mne.create_info(ch_names=['ECG-SYN'], sfreq=raw.info['sfreq'], ch_types=['mag']))
    ignore = ['ch_names', 'chs', 'nchan', 'bads']
    for k, v in raw.info.items():
        if k not in ignore:
            ecg_raw.info[k] = v
    
    ecg.add_channels([ecg_raw])
    
    if not events.size:
        #Load events and extract event ids
        event_file = glob.glob(op.join(epochs_path, 'Events', subject, '{}_{}_{}*.eve'.format(subject, state, block)))[0]
        events = mne.read_events(event_file)
    
    t_start = time.perf_counter()
    figs['{} Raw ECG'.format(state+block)] = ecg.plot(n_channels=len(ecg.ch_names), events=events, scalings='auto', bgcolor=None, title='Average pulse: {} bpm'.format(np.round(len(events[:,2]==999)*60/ecg.times[-1])))
#    figs['{} Raw ECG'.format(state+block)].set_size_inches(16,4)
    logger.info("Editing Report took {:.1f}s.".format(time.perf_counter() - t_start))
    
    # Save report
    report.add_figs_to_section(list(figs.values()), list(figs.keys()), state+block)
    if save_report:
        report_file = op.join(Analysis_path, task, 'meg', 'Reports', subject, '{}_ECG-report.html'.format(subject))
        os.makedirs(op.dirname(report_file), exist_ok=True)
        report.save(report_file, open_browser=False, overwrite=True)
    
    return report


def empty_room_covariance(task:str, subject:str, notch=inspect.signature(process).parameters['notch'].default, high_pass=inspect.signature(process).parameters['high_pass'].default, low_pass=inspect.signature(process).parameters['low_pass'].default):
    """
    Compute noise covariance matrix from empty room raw data.
    Parameters (notch, high_pass, low_pass) are set to the default of process(), but make sure they are the same as for preprocessing.
    Ouput: Analyses/<task>/meg/Covariance/<subject>/empty_room-cov.fif
    """
    # Load noise data
    raw_empty_room = mne.io.read_raw_ctf(get_rawpath(subject, noise=True), preload=True, clean_names=True)
    
    # Filter
    if notch.size:
        raw_empty_room.notch_filter(notch, fir_design='firwin', n_jobs=4)
    raw_empty_room.filter(l_freq=high_pass, h_freq=low_pass, fir_design='firwin', n_jobs=4)
    
    # Save noise covariance matrix
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', subject)
    if not op.exists(cov_path):
        os.makedirs(cov_path)
    
    noise_cov = mne.compute_raw_covariance(raw_empty_room, n_jobs=4)
    mne.write_cov(op.join(cov_path, 'empty_room-cov.fif'), noise_cov)


def auto_annotate(raw, window=1, overlap=0, decim=1):
    """
    Annotate Raw data based on AutoReject's global rejection threshold computation.
    'bad AutoReject' annotations are segments of length window [sec] overlapping by overlap [sec].
    
    Operates in-place.
    """
    epochs = mne.Epochs(raw.copy(), make_fixed_length_events(raw.copy(), 111, duration=window-overlap, first_samp=False), tmin=0, tmax=window, baseline=None, picks=mne.pick_types(raw.info, ref_meg=False), preload=True, reject=None)
    events = epochs.copy().events
    
    epochs_clean = epochs.copy()
    epochs_clean.drop_bad()
    already_rejected = set(np.where(epochs_clean.drop_log)[0])
    
    t_start = time.perf_counter()
    reject = get_rejection_threshold(epochs, decim=decim)
    logger.info("Rejection threshold computation took {:.1f}s.".format(time.perf_counter() - t_start))
    epochs.drop_bad(reject = reject)
    
    newly_rejected = set(np.where(epochs.drop_log)[0]) - already_rejected
    rejected = events[tuple(newly_rejected),0]
    
    new_annot = mne.Annotations(raw.times[rejected], np.full(rejected.shape, window), np.full(rejected.shape, 'bad AutoReject'))
    if raw.annotations:
        for onset, duration, description in zip(new_annot.onset, new_annot.duration, new_annot.description):
            raw.annotations.append(onset, duration, description)
    else:
        raw.annotations = new_annot
    
    return raw, reject


Pre = lambda x: detrend((x-np.mean(x))/np.std(x))

rejected_duration = lambda raw: raw.get_data().shape[-1]/raw.info['sfreq'] - raw.get_data(reject_by_annotation='omit').shape[-1]/raw.info['sfreq']
