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


def Pre(x):
    return detrend((x-np.mean(x))/np.std(x))

# # /!\ Custom attributes (e.g., ica.scores_) are not kept upon .save(), which calls _write_ica() whose dict ica_misc is not editable on call.
# # => exploit the attribute labels_
# # # /!\ Numpy arrays are not supported --> convert to type list with .tolist()

def run_ica(task, subject, state, block, raw=None, save=True, fit_ica=False, ecg_ica=True, n_components=0.975, method='fastica', ica_rejection={'mag':4000e-15}, ECG_threshold=0.25, ECG_max=3, EOG_threshold=3, EOG_min=1, EOG_max=2, custom_args=dict()):
    """
    Fit ICA on raw MEG data and return ICA object.
    If save, save ICA, save ECG and EOG artifact scores plots, and write log (default to True).
    If fit_ica, fit ICA even if there is already an ICA file (default to False).
    
    Output
    ------
    'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-ica.fif'
    
    'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_ecg.svg'
    
    'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_eog.svg'
    
    Log
    ---
    'Analyses/<task>/meg/ICA/ICA_log.tsv'
    
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
    ECG_threshold : 
        ECG artifact detection threshold (mne default to 0.25)
    EOG_threshold : 
        EOG artifact detection threshold (mne default to 3.0)
    ECG_max : 
        maximum number of ECG components to exclude, set to None for all ECG components detected (default to 3)
    EOG_min : 
        minimum number of EOG components to exclude (default to 1)
    EOG_max : 
        maximum number of EOG components to exclude (default to 2)
    
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
    
    raw.set_channel_types({get_chan_name(subject, 'ecg_chan', raw):'ecg', get_chan_name(subject, 'eogV_chan', raw):'eog', 'UPPT001':'stim'})
    
    # Crop recording
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[1][0]] > 300 else None
    raw.crop(tmin=start, tmax=end)
    
    raw.pick_types(meg=True, ref_meg=False, ecg=True, eog=True, exclude='bads')
    
    # Filter for ICA
    raw.filter(l_freq=1, h_freq=40, fir_design='firwin', n_jobs=4)
    
    # Fit ICA
    if fit_ica or not op.isfile(ICA_file):
        ica = ICA(n_components=n_components, method=method) #create ICA object
        ica.exclude = []
        ica.drop_inds_ = []
        ica.labels_ = dict()
        
        if ecg_ica:
            ecg_ica=ica.copy()
            if custom_args:
                ecg_epo, pulse = custom_ecg_epochs(raw.copy(), custom_args, tmin=-0.1, tmax=0.1)
            else:
                ecg_epo = create_ecg_epochs(raw.copy(), tmin=-0.1, tmax=0.1)
            ecg_ica.fit(ecg_epo, reject=ica_rejection, decim=6, picks=mne.pick_types(raw.info, meg=True)) #decimate: 200Hz is more than enough for ICA, saves time; picks: fit only on MEG
            ecg_ica.labels_['rejection'] = ica_rejection
            ecg_ica.labels_['drop_inds_'] = ecg_ica.drop_inds_
            
            if custom_args:
                _, scores, pulse = custom_bads_ecg(ecg_ica, ecg_epo, custom_args, threshold=ECG_threshold)
                ecg_ica.labels_['ecg_scores'] = scores.tolist()
            else:
                ecg_ica.labels_['ecg_scores'] = ecg_ica.find_bads_ecg(ecg_epo, threshold=ECG_threshold)[1].tolist()
                pulse = mne.preprocessing.find_ecg_events(raw.copy(), l_freq=8, h_freq=16, ch_name=get_chan_name(subject, 'ecg_chan', raw))[2]
            
            # Fix number of artifactual components
            ecg_ica.labels_['ecg'] = ica.labels_['ecg'][:ECG_max]
        
        ica.fit(raw, reject=ica_rejection, decim=6, picks=mne.pick_types(raw.info, meg=True)) #decimate: 200Hz is more than enough for ICA, saves time; picks: fit only on MEG
        ica.labels_['rejection'] = ica_rejection
        ica.labels_['drop_inds_'] = ica.drop_inds_
    
    else:
        ica = read_ica(ICA_file)
    
    # Detect ECG and EOG artifacts
    if custom_args:
        _, scores, pulse = custom_bads_ecg(ica, raw.copy(), custom_args, threshold=ECG_threshold)
        ica.labels_['ecg_scores'] = scores.tolist()
    else:
        ica.labels_['ecg_scores'] = ica.find_bads_ecg(raw.copy(), threshold=ECG_threshold)[1].tolist()
        pulse = mne.preprocessing.find_ecg_events(raw.copy(), l_freq=8, h_freq=16, ch_name=get_chan_name(subject, 'ecg_chan', raw))[2]
    
    ica.labels_['eog_scores'] = ica.find_bads_eog(raw.copy(), threshold=EOG_threshold)[1].tolist()
    
    # Fix number of artifactual components
    ica.labels_['ecg'] = ica.labels_['ecg'][:ECG_max]
    ica.labels_['eog'] = ica.labels_['eog'][:EOG_max]
    if EOG_min and not ica.labels_['eog']:
        ica.labels_['eog'] = np.argsort(np.abs(ica.labels_['eog_scores'])).tolist()
        ica.labels_['eog'] = ica.labels_['eog'][::-1][:EOG_min]
    
    # Tag for exclusion
    ica.exclude = ica.labels_['ecg'] + ica.labels_['eog']
    
    # Plot scores
    ica.plot_scores(ica.labels_['ecg_scores'], exclude=ica.labels_['ecg'], labels='ecg', axhline=ECG_threshold)
    if save:
        plt.savefig(op.join(ICA_path, '{}_{}-{}_components-scores_ecg.pdf'.format(state, block, n_components)), transparent=True)
        plt.close()
    
    ica.plot_scores(ica.labels_['eog_scores'], exclude=ica.labels_['eog'], labels='eog')
    if save:
        plt.savefig(op.join(ICA_path, '{}_{}-{}_components-scores_eog.pdf'.format(state, block, n_components)), transparent=True)
        plt.close()
    
    # Save ICA
    if save:
        ica.save(ICA_file)
        # Write ICA log
        with open(ICA_log, 'a') as fid:
            fid.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(time.strftime('%Y_%m_%d\t%H:%M:%S',time.localtime()),subject,state,block,int(round(start)),int(round(end)) if end else int(round(raw.times[-1])),n_components,ica.n_components_,len(ica.labels_['ecg']),int(round(pulse)),len(ica.labels_['eog']),ica.labels_['rejection'],len(ica.labels_['drop_inds_'])))
    
    return ica


def process(task, subject, state, block, n_components=.975, ica=None, check_ica=True, save_ica=True, overwrite_ica=False, fit_ica=False, ica_rejection={'mag':4000e-15}, notch=np.arange(50,301,50), high_pass=0.5, low_pass=None, ECG_threshold=0.25, EOG_threshold=3, custom_args=dict(), update_HPI=True, precision='0.5cm', opt='start'):
    """
    Run preprocessing and return preprocessed raw data.
    If check_ica, plot overlay and properties of ECG and EOG components (default to True).
    If save_ica, save ICA and ICA plots (deleting previously existing component properties), and write ICA log (default to True).
    If overwrite_ica, artifact scoring will be perfomed again (default to False).
    Output:
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_ecg.png'
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-scores_eog.png'
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-<properties>_ecg<i>.pdf'
        'Analyses/<task>/meg/ICA/<subject>/<state>_<block>-<n>_components-<properties>_eog<j>.pdf'
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
    ICA_file = op.join(ICA_path, '{}_{}-{}_components-ica.fif'.format(state, block, n_components))
    if not ica:
        if overwrite_ica or not op.isfile(ICA_file):
            ica = run_ica(task, subject, state, block, raw=raw.copy(), n_components=n_components, save=save_ica, fit_ica=fit_ica, ica_rejection=ica_rejection, ECG_threshold=ECG_threshold, EOG_threshold=EOG_threshold, custom_args=custom_args)
        else:
            ica = read_ica(ICA_file)
    
    raw.info['subject_info'].update({'sub':subject})
    raw.set_channel_types({get_chan_name(subject, 'ecg_chan', raw):'ecg', get_chan_name(subject, 'eogV_chan', raw):'eog', 'UPPT001':'stim'})
    raw.pick_types(meg=True, ecg=True, eog=True, stim=True, exclude='bads')
    
    # Crop recording
    events = mne.find_events(raw)
    start = raw.times[events[0][0]] if raw.times[events[0][0]] < 120 else 0
    end = raw.times[events[-1][0]] if len(events) > 1 and raw.times[events[-1][0]] > 300 else None
    raw.crop(tmin=start, tmax=end)
    
    # Filter
    if notch.any():
        raw.notch_filter(notch, fir_design='firwin', n_jobs=4)
    raw.filter(l_freq=high_pass, h_freq=low_pass, fir_design='firwin', n_jobs=4)
    
    # Visual check
    if check_ica:
        # ECG components
        if custom_args:
            check_ecg, pulse = custom_ecg_epochs(raw, custom_args, reject=ica_rejection)
        else:
            check_ecg = create_ecg_epochs(raw, reject=ica_rejection)
        ica.plot_overlay(check_ecg.average(), exclude=ica.labels_['ecg'])
        if save_ica:
            plt.savefig(op.join(ICA_path, '{}_{}-{}_components-overlay_ecg.png'.format(state, block, n_components)), transparent=True, dpi=360)
            plt.close()
            for pic in glob.glob(op.join(ICA_path, state+'*'+block+'*properties_ecg*')):
                os.remove(pic)
        
        for comp in ica.labels_['ecg']:
            ica.plot_properties(check_ecg, picks=comp)
            if save_ica:
                plt.savefig(op.join(ICA_path, '{}_{}-{}_components-properties_ecg{}.pdf'.format(state, block, n_components, comp)), transparent=True)
                plt.close()
        
        # EOG components
        check_eog = create_eog_epochs(raw, reject=ica_rejection)
        ica.plot_overlay(check_eog.average(), exclude=ica.labels_['eog'])
        if save_ica:
            plt.savefig(op.join(ICA_path, '{}_{}-{}_components-overlay_eog.png'.format(state, block, n_components)), transparent=True, dpi=360)
            plt.close()
            for pic in glob.glob(op.join(ICA_path, state+'*'+block+'*properties_eog*')):
                os.remove(pic)
        
        for comp in ica.labels_['eog']:
            ica.plot_properties(check_eog, picks=comp)
            if save_ica:
                plt.savefig(op.join(ICA_path, '{}_{}-{}_components-properties_eog{}.pdf'.format(state, block, n_components, comp)), transparent=True)
                plt.close()
    
    # Save pre-processed data with updated head coordinates
    if update_HPI:
        raw = HPI_update(task, subject, block, raw.copy(), precision=precision, opt=opt, reject_head_mvt=True)
        raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}-raw.fif'.format(state, block))
        if not op.isdir(op.dirname(raw_file)):
            os.makedirs(op.dirname(raw_file))
        raw.save(raw_file, overwrite=True)
        
        return raw
    
    # Apply ICA
    raw_ECG = raw.copy()
    ica.exclude = ica.labels_['eog']
    ica.apply(raw_ECG)
    
    ica.exclude = ica.labels_['ecg'] + ica.labels_['eog']
    ica.apply(raw)
    
    return {'ECG_exlcuded': raw, 'ECG_included': raw_ECG}


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
        for i in range(len(bad_segments)):
            data.annotations = mne.Annotations(bad_segments.get_values()[:,0], np.diff(bad_segments.get_values()).ravel(), ['bad head mvt'] * bad_segments.get_values().shape[0])
#            events = data.events[:,0]/data.info['sfreq'] #extract event timing
#            data.drop(np.logical_and(bad_segments.iloc[i][0] < events + data.times[-1], events + data.times[0] < bad_segments.iloc[i][1]), reason = 'head_movement')
#            #reject the epoch if it ends after the beginning of the bad segment, and starts before the end of the bad segment
    
    return data


def epoch(task, subject, state, block, raw, name, events, event_id, save=True, rejection=None, save_t_timing=False, sliding=False, tmin=-.5, tmax=.8, baseline=None, update_HPI=True, precision='0.5cm', opt='start', reject_head_mvt=True, check_ica=False, overwrite_ica=False, fit_ica=False, ica_rejection={'mag':4000e-15}, notch=np.arange(50,301,50), high_pass=0.5, low_pass=None, ECG_threshold=0.25, EOG_threshold=3, custom_ecg=dict()):
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
        names: list of events to epoch
        tmin, tmax, baseline: epoching parameters
    """
    # Save path
    epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs', subject)
    if not op.exists(epochs_path):
        os.makedirs(epochs_path)
    epochs_file = op.join(epochs_path, '{}_{}-{}-epo.fif'.format(state, block, name))
    
#    check_ecg_epoch(task, subject, state, block, raw_ECG.copy(), events=np.concatenate([R_epochs.events, T_events]), save=save)
    
    picks = mne.pick_types(raw.info, meg=True, ecg=True, eog=True, stim=True, exclude='bads')    
    epochs = Epochs(raw, events, event_id, tmin=tmin, tmax=tmax, reject=rejection, baseline=baseline, preload=True, picks=picks)
    
    if rejection:
        epochs.drop_bad()
        epochs.plot_drop_log()
        if save:
            plt.savefig(op.join(epochs_path, '{}_{}-{}-drop_log.pdf'.format(state, block, name)), transparent=True)
            plt.close()
            drop_log = op.join(Analysis_path, task, 'meg', 'Epochs', 'drop_log.txt')
            with open(drop_log, 'a') as fid:
                fid.write('{} {} epochs dropped\t{}\n'.format(epochs_file.split('/')[-2:], len(np.array(epochs[name].drop_log)[np.where(epochs[name].drop_log)]), rejection))
    
    if update_HPI:
        epochs = HPI_update(task, subject, block, data=epochs.copy(), precision=precision, opt=opt, reject_head_mvt=reject_head_mvt)
    
    # Save epochs
    if save:
        epochs.save(epochs_file)
    print(colored (epochs_file, 'green'))
        
    return epochs


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
        raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}-raw.fif'.format(state, block))
        raw = mne.io.read_raw_fif(raw_file, preload=True)
    
    # Epoch ECG on R peak
    ecg = raw.copy().pick_types(meg=False, ref_meg=False, ecg=True)
    if custom_args:
        epo, pulse = custom_ecg_epochs(ecg, custom_args, reject_by_annotation=False, event_id=R_id)
    else:
        epo = create_ecg_epochs(ecg, reject_by_annotation=False, event_id=R_id)
    
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
    event_file = op.join(Analysis_path, task, 'meg', 'Epochs', 'Events', subject, '{}_{}+R_{}+T_{}.eve'.format(state, block, R_id, T_id))
    os.makedirs(op.dirname(event_file), exist_ok=True)
    mne.write_events(event_file, events)
    
    return events, event_id


def ECG_ICA(task, subject, state, block, raw=None, events=None, event_id=None, rejection={'mag': 7000e-15}, Rwin=-.1, Twin=.1):
    """
    
    """
    # Load data
    if not raw:
        raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}-raw.fif'.format(state, block))
        raw = mne.io.read_raw_fif(raw_file, preload=True)
    
    rawforica = raw.copy()
    rawforica.filter(l_freq=1, h_freq=40, fir_design='firwin', n_jobs=4)
    if not events or not event_id:
        #Load events and extract event ids
        event_file = glob.glob(op.join(Analysis_path, task, 'meg', 'Epochs', 'Events', subject, '{}_{}+R_*+T_*.eve'.format(state, block)))[0]
        events = mne.read_events(event_file)
        
        ids = op.splitext(op.basename(event_file))[0].split('+')[1:]
        event_id = dict()
        for i in ids:
            k, v = i.split('_')
            event_id[k] = int(v)
    
    epoforica = mne.Epochs(rawforica, events, event_id, baseline=None, reject=rejection, reject_by_annotation=True, preload=True)
    epoforica.equalize_event_counts(list(event_id.keys()))
    RT_delay_i = int(np.round(np.median(epoforica['T'].events[:,0] - epoforica['R'].events[:,0])))
    RT_delay = epoforica.times[epoforica.time_as_index(0)[0] + RT_delay_i]
    
    epoforica = epoforica['R']
    epoforica.crop(Rwin, RT_delay + Twin)
    epoforica.average().plot_joint()
    
    ica = ICA(n_components=n_components, method=method) #create ICA object
    ica.exclude = []
    ica.labels_ = dict()
    
    ica.fit(epoforica, picks=mne.pick_types(epoforica.info, ref_meg=False), decim=6) #decimate: 200Hz is more than enough for ICA, saves time; picks: fit only on MEG
    
    # Detect ECG artifacts
    ica.labels_['ecg_scores'] = ica.find_bads_ecg(epoforica, threshold=0)[1].tolist()
    
    # Fix number of artifactual components
    ica.labels_['ecg'] = ica.labels_['ecg'][:3]
    
    # Tag for exclusion
    ica.exclude = ica.labels_['ecg']
    
    # Plot scores
    ica.plot_scores(ica.labels_['ecg_scores'], exclude=ica.labels_['ecg'], labels='ecg', axhline=ECG_threshold)
    if save:
        plt.savefig(op.join(ICA_path, '{}_{}-{}_components-scores_ecg.pdf'.format(state, block, n_components)), transparent=True)
        plt.close()
    
    epochs = mne.Epochs(raw, events, event_id, baseline=None, reject=rejection, reject_by_annotation=True, preload=True)['R']
    ica.plot_overlay(epochs.average())
    


def t_detector(task, subject, state, block, raw, R_sign=0, event_id=333, l_freq=5, h_freq=35, T_window=[.2,.35], save=True, custom_ecg=dict()):
    """
    From raw data containing at least the ECG channel, returns events corresponding to T peaks. If save=True, saves their timing in 'Analyses/<task>/meg/Epochs/T_timing.tsv'.
    """
    # Save T peak timing
    timing_file = op.join(Analysis_path, 'MEG', 'meta', 'T_timing-{}_{}-{}.tsv'.format(l_freq, h_freq, T_window))
    if save and not op.isfile(timing_file):
        with open(timing_file, 'w') as fid:
            fid.write('subject\tstate\tblock\tR_peak\tT_delay\n')
    
    # Pick ECG
    ECG_channel = get_chan_name(subject, 'ecg_chan', raw)
    raw.pick_channels([ECG_channel])
    
    # Set ECG_channel as EEG data and its copy 'ECG' as ECG
    ecg = raw.copy()
    raw.set_channel_types({ECG_channel:'eeg'})
    ecg.rename_channels({ECG_channel:'ECG'})
    ecg.set_channel_types({'ECG':'ecg'})
    raw.add_channels([ecg])
    
    # Filter EEG as in qrs_detector
    if l_freq and h_freq:
        raw.filter(l_freq=l_freq, h_freq=h_freq, filter_length='10s', l_trans_bandwidth=.5, h_trans_bandwidth=.5, phase='zero-double', fir_window='hann', fir_design='firwin2')
    
    # Find R peak events
    if subject in custom_ecg.keys():
        custom_args = custom_ecg[subject].copy()
        if 'tstart' in custom_args.keys():
            custom_args['tstart'] = custom_args['tstart'][state+block]
        R_sign = custom_args['T_sign'] if 'T_sign' in custom_args.keys() else custom_args['R_sign']
        R_epochs, pulse = custom_ecg_epochs(raw, custom_args)
    else:
        R_epochs = create_ecg_epochs(raw)
    
    # Find T peak timing
    R_epochs.pick_channels([ECG_channel])
    data = R_epochs.get_data()[:, 0, :] #The indices delete the axis=1 that is not used.
    
    R_time_i = R_epochs.time_as_index(0)
    if not R_sign:
        R_sign = np.median(np.sign(Pre(data)[:, R_time_i]))
    print('R_sign =', R_sign)
    T_window_i = R_epochs.time_as_index(T_window)
    
    T_times_i = (R_sign*data[:, T_window_i[0]:T_window_i[1]]).argmax(axis=1) + T_window_i[0]
    T_times = R_epochs.times[T_times_i]
    
    # Save timing of R peaks and delay until T peak (in sec)
    if save:
        for i in range(len(T_times)):
            with open(timing_file, 'a') as fid:
                fid.write("{}\t{}\t{}\t{}\t{}\n".format(subject, state, block, np.round(raw.times[R_epochs.events[i,0]],3), np.round(T_times[i], 3)))
    
    # Return T peak events
    T_events = R_epochs.copy().events
    T_events[:,0] += T_times_i - R_time_i
    T_events[:,2] = event_id
    
    return T_events,T_times,R_epochs,raw


def t_detector_sliding(task, subject, state, block, raw, R_sign=0, event_id=333, l_freq=5, h_freq=35, T_window=[.2,.35], slide=[.1,.5], step=.05, save=False, custom_ecg=dict()):
    """
    From raw data containing at least the ECG channel, returns events corresponding to T peaks. If save=True, saves their timing in 'Analyses/<task>/meg/Epochs/T_timing.tsv'.
    """
    # Save T peak timing
    timing_file = op.join(Analysis_path, 'MEG', 'meta', 'T_timing-{}_{}-{}_sliding.tsv'.format(l_freq, h_freq, T_window))
    if save and not op.isfile(timing_file):
        with open(timing_file, 'w') as fid:
            fid.write('subject\tstate\tblock\tR_peak\tT_delay\n')
    
    # Pick ECG
    ECG_channel = get_chan_name(subject, 'ecg_chan', raw)
    raw.pick_channels([ECG_channel])
    
    # Set ECG_channel as EEG data and its copy 'ECG' as ECG
    ecg = raw.copy()
    raw.set_channel_types({ECG_channel:'eeg'})
    ecg.rename_channels({ECG_channel:'ECG'})
    ecg.set_channel_types({'ECG':'ecg'})
    raw.add_channels([ecg])
    
    # Filter EEG as in qrs_detector
    if l_freq and h_freq:
        raw.filter(l_freq=l_freq, h_freq=h_freq, filter_length='10s', l_trans_bandwidth=.5, h_trans_bandwidth=.5, phase='zero-double', fir_window='hann', fir_design='firwin2')
    
    # Find R peak events
    if subject in custom_ecg.keys():
        custom_args = custom_ecg[subject].copy()
        if 'tstart' in custom_args.keys():
            custom_args['tstart'] = custom_args['tstart'][state+block]
        R_sign = custom_args['T_sign'] if 'T_sign' in custom_args.keys() else custom_args['R_sign']
        R_epochs, pulse = custom_ecg_epochs(raw, custom_args)
    else:
        R_epochs = create_ecg_epochs(raw)
    
    # Find T peak timing
    R_epochs.pick_channels([ECG_channel])
    data = R_epochs.get_data()[:, 0, :] #The indices delete the axis=1 that is not used.
    
    R_time_i = R_epochs.time_as_index(0)
    if not R_sign:
        R_sign = np.median(np.sign(Pre(data)[:, R_time_i]))
    print('R_sign =', R_sign)
    
    win_size = np.diff(T_window)
    peaks = []
    windows = []
    while slide[0] < T_window[0]:
        T_window_i = R_epochs.time_as_index([slide[0], T_window[1]])
        print(T_window_i)
        T_times_i = (R_sign*data[:, T_window_i[0]:T_window_i[1]]).argmax(axis=1) + T_window_i[0]
        print(T_times_i.size)
        peaks.append(T_times_i.mean())
#        peaks.append(np.median(T_times_i))
        windows.append(T_window_i)
        slide[0] += step
    
    while T_window[1] <= slide[1]:
        T_window_i = R_epochs.time_as_index([T_window[0], T_window[1]])
        print(T_window_i)
        T_times_i = (R_sign*data[:, T_window_i[0]:T_window_i[1]]).argmax(axis=1) + T_window_i[0]
        print(T_times_i.size)
        peaks.append(T_times_i.mean())
#        peaks.append(np.median(T_times_i))
        windows.append(T_window_i)
        T_window[1] += step
    
    opti_peak = np.median(peaks)
    print(peaks, opti_peak)
    
    T_window_i = windows[np.abs(np.subtract(peaks, opti_peak)).argmin()]
    T_times_i = (R_sign*data[:, T_window_i[0]:T_window_i[1]]).argmax(axis=1) + T_window_i[0]
    T_times = R_epochs.times[T_times_i]
    
    # Save timing of R peaks and delay until T peak (in sec)
    if save:
        for i in range(len(T_times)):
            with open(timing_file, 'a') as fid:
                fid.write("{}\t{}\t{}\t{}\t{}\n".format(subject, state, block, np.round(raw.times[R_epochs.events[i,0]],3), np.round(T_times[i], 3)))
    
    # Return T peak events
    T_events = R_epochs.copy().events
    T_events[:,0] += T_times_i - R_time_i
    T_events[:,2] = event_id
    
    return T_events,T_times,R_epochs,raw


def check_ecg_epoch(task, subject, state, block, raw=None, events=np.array([]), n_components=.975, synthetic=True, save=False):
    """
    Plot ECG along with events used for epoching. If synthetic, create and plot a synthetic ECG channel as well.
    """
    epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs')
    os.makedirs(op.join(epochs_path, subject), exist_ok=True)
    
    # Load data
    if not raw:
        raw_file = op.join(Analysis_path, task, 'meg', 'Raw', subject, '{}_{}-raw.fif'.format(state, block))
        raw = mne.io.read_raw_fif(raw_file, preload=True)
    
    ecg = raw.copy().pick_types(meg=False, ref_meg=False, ecg=True)
    
    if synthetic:
        ica = read_ica(op.join(Analysis_path, task, 'meg', 'ICA', subject, '{}_{}-{}_components-ica.fif'.format(state, block, n_components)))
        ica.exclude = ica.labels_['eog']
        ica.apply(raw)
        
        # Add synthetic channel as in MNE's create_ecg_epochs():
        ecg_syn = raw.get_data(picks=mne.pick_types(raw.info, meg='mag', ref_meg=False)).mean(axis=0)
        ecg_raw = mne.io.RawArray(ecg_syn[None], mne.create_info(ch_names=['ECG-SYN'], sfreq=raw.info['sfreq'], ch_types=['mag']))
        ignore = ['ch_names', 'chs', 'nchan', 'bads']
        for k, v in raw.info.items():
            if k not in ignore:
                ecg_raw.info[k] = v
        
        ecg.add_channels([ecg_raw])
    
    if not events.any():
        #Load events and extract event ids
        event_file = glob.glob(op.join(epochs_path, 'Events', subject, '{}_{}*.eve'.format(state, block,)))[0]
        events = mne.read_events(event_file)
    
    ecg.plot(n_channels=len(ecg.ch_names), events=events, scalings='auto', bgcolor=None)
    if save:
        plt.savefig(op.join(epochs_path, subject, '{}_{}-ECG.pdf'.format(state, block)), transparent=True)
        plt.close()
    
    return ecg


def empty_room_covariance(task:str, subject:str, notch=inspect.signature(process).parameters['notch'].default, high_pass=inspect.signature(process).parameters['high_pass'].default, low_pass=inspect.signature(process).parameters['low_pass'].default):
    """
    Compute noise covariance matrix from empty room raw data.
    Parameters (notch, high_pass, low_pass) are set to the default of process(), but make sure they are the same as for preprocessing.
    Ouput: Analyses/<task>/meg/Covariance/<subject>/empty_room-cov.fif
    """
    # Load noise data
    noise_path = op.join(Raw_data_path, 'MEG_noise')
    raw_empty_room = mne.io.read_raw_ctf(op.join(noise_path, * get_rawpath(subject, noise=1)), preload=True)
    
    # Filter
    if notch.any():
        raw_empty_room.notch_filter(notch, fir_design='firwin', n_jobs=4)
    raw_empty_room.filter(l_freq=high_pass, h_freq=low_pass, fir_design='firwin', n_jobs=4)
    
    # Save noise covariance matrix
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', subject)
    if not op.exists(cov_path):
        os.makedirs(cov_path)
    
    noise_cov = mne.compute_raw_covariance(raw_empty_room, n_jobs=4)
    mne.write_cov(op.join(cov_path, 'empty_room-cov.fif'), noise_cov)
