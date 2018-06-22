# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 17:42:45 2017

@author: oussama.abdoun
"""

#==============================================================================
from header import *
from tqdm import tqdm
#==============================================================================
warnings.filterwarnings("ignore",category=DeprecationWarning)

# MANUAL CHECK
#==============================================================================
#task='SMEG'; names=['ECG_included','ECG_excluded']; keys=['R','T']
#subject='069'; evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
#state='FA'; block='02'; name=names[0]; k=keys[0]
#evoked = mne.Evoked(op.join(evoked_path, '{}_{}-{}-ave.fif'.format(state, block, name)), condition = k)
#bad_chan = get_chan_name(subject, 'bad', data=evoked)
#evoked.drop_channels(bad_chan).plot_joint()
#evoked.plot_sensors()
#
## Save ERP
#for state in ['RS','FA','OM']:
#    for block in get_blocks(subject, task=task, state=state):
#        mne.read_evokeds(op.join(evoked_path, '{}_{}-{}-ave.fif'.format(state, block, name)))[0].drop_channels(bad_chan).plot_joint(title = 'Subject {} - {} {}\n{}\n{} channels rejected'.format(subject, state, block, name, len(bad_chan)))
#        plt.savefig(op.join(evoked_path, '{}_{}-{}-ERP_rejected.pdf'.format(state, block, name)), transparent=True)
#        plt.close()
#==============================================================================


def ERP(task, subject, state, block_group, name='ECG_included', cov_keys=['R','T'], rejection={'mag':3500e-15}, baseline={'R':(None,-.225), 'T':(None,None)}, window={'R':(-.35,.425), 'T':(-.175,.275)}, data_win={'R':(.35,None), 'T':(.1,.25)}, n_components=.975):
    """
    Returns baseline noise covariance for the block_group and a dict containing a list of Evoked (for each block in block_group) assigned to their name.
    Parameters:
        block_group: list of coregistered blocks of the same state to be combined
        rejection: epoch rejection threshold (default to 2500 fT for magnetometers)
        baseline: baseline to apply as a tuple (if not provided, (None,0) will be used)
        t_delay: for T peak Epochs, shift the baseline by subtracting t_delay
        names: list of Epoch name
    Output:
        Analyses/<task>/meg/Covariance/<subject>/<name>-<state>_<block_group>-cov.fif (Baseline covariance for the block_group)
        Analyses/<task>/meg/Evoked/<subject>/<name>-<state>_<block>-ave.fif (Evoked per block)
    """
    # Define pathes
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', subject)
    if not op.exists(cov_path):
        os.makedirs(cov_path)
    raw_path = op.join(Analysis_path, task, 'meg', 'Raw', subject)
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', subject)
    epochs_path = op.join(Analysis_path, task, 'meg', 'Epochs')
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    if not op.exists(evoked_path):
        os.makedirs(evoked_path)
#    ICA_file = op.join(Analysis_path, task, 'meg', 'ICA', subject, '{}_{}-{}_components-ica.fif'.format(state, block, n_components))
#    ica = read_ica(ICA_file)
    
    epochs_list = []
    evoked_list = dict()
    evoked = dict()
    
    for b,block in enumerate(block_group):
        # Load raw and apply ica
        raw = mne.io.read_raw_fif(op.join(raw_path, '{}_{}-raw.fif'.format(state, block)), preload=True)
        ica = read_ica(op.join(ICA_path, '{}_{}-{}_components-ica.fif'.format(state, block, n_components)))
        ica.exclude = ica.labels_['eog']
        if name != 'ECG_included':
            ica.exclude += ica.labels_['ecg']
        ica.apply(raw)
        
        #Load events and extract event ids
        event_file = glob.glob(op.join(epochs_path, 'Events', subject, '{}_{}+{}_*.eve'.format(state, block, '_*+'.join(cov_keys))))[0]
        events = mne.read_events(event_file)
        
        ids = op.splitext(op.basename(event_file))[0].split('+')[1:]
        event_id = dict()
        for i in ids:
            k, v = i.split('_')
            event_id[k] = int(v)
        
#        #Epoch
#        epochs = mne.Epochs(raw, events, event_id, tmin, tmax, baseline=None, reject=rejection, reject_by_annotation=True, preload=True)
#        epochs_list.append(epochs)
        
        for k in event_id.keys():
            if not k in evoked_list.keys(): evoked_list[k] = []
            #Epoch
            epochs = mne.Epochs(raw, events, event_id, *window[k], baseline=baseline[k], reject=rejection, reject_by_annotation=True, preload=True)
            data_cov_file = op.join(cov_path, '{}_{}_{}-{}_{}-data-cov.fif'.format(subject, state, '_'.join(block_group), k, name))
            mne.write_cov(data_cov_file, mne.compute_covariance(epochs, True, *data_win[k], method='auto', n_jobs=4))
            
            if k in cov_keys: epochs_list.append(epochs)
            
#            # Apply baseline
#            epochs[k].apply_baseline(baseline[k])
            
            # Save Evoked
            evoked[k] = epochs[k].average()
            evoked_list[k].append(evoked[k])
        
        evoked_file = op.join(evoked_path, '{}_{}_{}-{}-ave.fif'.format(subject, state, block, name))
        mne.write_evokeds(evoked_file, [evo for evo in evoked.values()])
        
        # Plot and save drop log
        epochs.plot_drop_log(subject=subject)
        plt.savefig(op.join(epochs_path, subject, '{}_{}_{}-{}-drop_log.pdf'.format(subject, state, block, name)), transparent=True)
        plt.close()
        
        drop_log = op.join(Analysis_path, task, 'meg', 'Epochs', 'drop_log.txt')
        with open(drop_log, 'a') as fid:
            fid.write('{} {} epochs dropped\t{}\n'.format(evoked_file.split('/')[-2:], len(np.array(epochs.drop_log)[np.where(epochs.drop_log)]), rejection))
        
        for k in event_id.keys():
            # Plot ERP time course
            evoked[k].plot_joint(title = 'Subject {} - {} {}\n{} {}'.format(subject, state, block, k, name))
            plt.savefig(op.join(evoked_path, '{}_{}_{}-{}_{}-ERP.pdf'.format(subject, state, block, k, name)), transparent=True)
            plt.close()
            
            # Same without bad channels
            bad_chan = list(set(evoked[k].info['bads']) | set(get_chan_name(subject, 'bad', evoked[k])))
            if bad_chan:
                evo = mne.pick_channels_evoked(evoked[k], exclude=bad_chan)
                evo.plot_joint(title = 'Subject {} - {} {}\n{} {}\n{} channels rejected'.format(subject, state, block, k, name, len(bad_chan)))
                plt.savefig(op.join(evoked_path, '{}_{}_{}-{}_{}-ERP_rejected.pdf'.format(subject, state, block, k, name)), transparent=True)
                plt.close()
    
    # Compute and save noise covariance
    all_epochs = mne.concatenate_epochs(epochs_list)
    noise_cov = mne.compute_covariance(all_epochs, keep_sample_mean=False, n_jobs=4)
    
    noise_cov_file = op.join(cov_path, '{}_{}_{}-{}-baseline_noise-cov.fif'.format(subject, state, '_'.join(block_group), name))
    mne.write_cov(noise_cov_file, noise_cov)
    
    return noise_cov, evoked_list


def src_rec(task, subject, state, block_group, evoked=dict(), noise_cov=None, keys=['R','T'], name='ECG_included', window={'R':(.35,.5), 'T':(.1,.25)}, surface='ico4', volume=6.2, bem_spacing='ico4', minimu_norm=True, compute_stc=True, fsaverage=True, baseline_cov=True, mindist=5, method:"'MNE', 'dSPM', 'sLORETA', 'eLORETA'"="dSPM"):
    """
    If compute_fwd=True, compute and save forward solution (overwriting previously existing one).
    If compute_inv=True, compute and save inverse solution.
    If compute_stc=True, compute and save SourceEstimate.
    If fsaverage=True, morph SourceEstimate to fsaverage and save it.
    Do this for surface if provided and for volume if provided.
    Parameters:
        block_group: list of coregistered blocks of the same state to be combined
        evoked: dict of Evoked assigned to their name as returned by baseline_covariance() (will be loaded if not provided)
        noise_cov: if baseline_cov=True, should be baseline covariance, else should be empty room noise covariance (will be loaded if not provided)
        bem_spacing: BEM surface downsampling as specified for anat.BEM()
        surface: source space subdivision as specified for anat.src_space()
        volume: distance between volume sources as `pos` specified for anat.src_space()
        names: list of Epoch name
        mindist: minimum distance (in mm) of sources from inner skull surface (see mne.make_forward_solution)
        method: source reconstruction method (see mne.minimum_norm.apply_inverse)
    Output (Analyses/<task>/meg/SourceEstimate/<subject>/):
        *'-fwd.fif' (forward model)
        *'-inv.fif' (inverse model)
        *'-lh.stc' (left hemisphere SourceEstimate)
        *'-rh.stc' (right hemisphere SourceEstimate)
        *'-fsaverage-lh.stc' (left hemisphere morphed SourceEstimate)
        *'-fsaverage-rh.stc' (right hemisphere morphed SourceEstimate)
    """
    # Define pathes
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', subject)
    stc_path = op.join(Analysis_path, task, 'meg', 'SourceEstimate', subject)
    if not op.exists(stc_path):
        os.makedirs(stc_path)
    trans_file = glob.glob(op.join(Analysis_path, task, 'meg', 'Coregistration', subject, '*'+block_group[0]+'*-trans.fif'))[0]
    
    # Load source space
    bem_sol = mne.read_bem_solution(op.join(os.environ['SUBJECTS_DIR'], subject, 'bem', '{}_{}_bem_solution.fif'.format(subject, (bem_spacing if bem_spacing else 'full'))))
    if surface:
        src_surf = mne.read_source_spaces(op.join(os.environ['SUBJECTS_DIR'], subject, 'src', '{}_{}_surface-src.fif'.format(subject, surface)))
    if volume:
        src_vol = mne.read_source_spaces(op.join(os.environ['SUBJECTS_DIR'], subject, 'src', '{}_{}_volume-src.fif'.format(subject, volume)))
    
    if not noise_cov:
        # Load baseline covariance if appropriate
        if baseline_cov:
            noise_cov = mne.read_cov(op.join(cov_path, '{}_{}_{}-{}-baseline_noise-cov.fif'.format(subject, state, '_'.join(block_group), name)))
        # Load empty room noise covariance if appropriate
        else:
            noise_cov = mne.read_cov(op.join(cov_path, 'empty_room-cov.fif'))
    
    # Evokeds will be loaded if not provided
    load_evoked = True if not evoked else False
    
    for k in keys:
        # Load Evokeds if not provided and combine the block_group
        if load_evoked:
            evoked[k] = mne.combine_evoked([mne.Evoked(op.join(evoked_path, '{}_{}_{}-{}-ave.fif'.format(subject, state, block, name)), condition = k) for block in block_group], 'nave')
            evoked[k].crop(*window[k])
    
    k0 = keys[0]
    
    # Bad channels will be rejected when computing the inverse operator and SourceEstimate
    bad_chan = list(set(evoked[k0].info['bads']) | set(get_chan_name(subject, 'bad', evoked[k0])))
    
    # Do surface SourceEstimate
    if surface:
        # Compute forward solution
        fwd_surf = mne.make_forward_solution(evoked[k0].info, trans=trans_file, src=src_surf, bem=bem_sol, meg=True, eeg=False, mindist=mindist, n_jobs=4)
        
        # Compute inverse solution
        if method is 'beamformer':
            for k in keys:
                data_cov = mne.read_cov(op.join(cov_path, '{}_{}_{}-{}_{}-data-cov.fif'.format(subject, state, '_'.join(block_group), k, name)))
                filters = mne.beamformer.make_lcmv(mne.pick_channels_evoked(evoked[k], exclude=bad_chan).info, mne.pick_channels_forward(fwd_surf, exclude=bad_chan), mne.pick_channels_cov(data_cov, exclude=bad_chan), noise_cov=mne.pick_channels_cov(noise_cov, exclude=bad_chan), pick_ori='vector')
                
                # Compute SourceEstimate
                stc_surf_file = op.join(stc_path, '{}_{}_{}-{}_{}-{}-surface_{}-{}_{}-{}-stc.h5'.format(subject, state, '_'.join(block_group), k, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window[k], method))
                stc_surf = mne.beamformer.apply_lcmv(mne.pick_channels_evoked(evoked[k], exclude=bad_chan), filters)
                stc_surf.save(stc_surf_file)
        
        else:
            inv_surf = make_inverse_operator(mne.pick_channels_evoked(evoked[k0], exclude=bad_chan).info, mne.pick_channels_forward(fwd_surf, exclude=bad_chan), mne.pick_channels_cov(noise_cov, exclude=bad_chan), loose=.2)
            
            for k in keys:
                # Compute SourceEstimate
                stc_surf_file = op.join(stc_path, '{}_{}_{}-{}_{}-{}-surface_{}-{}_{}-{}-stc.h5'.format(subject, state, '_'.join(block_group), k, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window[k], method))
                stc_surf = apply_inverse(mne.pick_channels_evoked(evoked[k], exclude=bad_chan), inv_surf, method=method, pick_ori='vector')
                stc_surf.save(stc_surf_file)
    
    # Do volume SourceEstimate
    if volume:
        # Compute forward solution
            fwd_vol = mne.make_forward_solution(evoked[k0].info, trans=trans_file, src=src_vol, bem=bem_sol, meg=True, eeg=False, mindist=mindist, n_jobs=4)
        
        # Compute inverse solution
        if method is 'beamformer':
            for k in keys:
                data_cov = mne.read_cov(op.join(cov_path, '{}_{}_{}-{}_{}-data-cov.fif'.format(subject, state, '_'.join(block_group), k, name)))
                filters = mne.beamformer.make_lcmv(mne.pick_channels_evoked(evoked[k], exclude=bad_chan).info, mne.pick_channels_forward(fwd_surf, exclude=bad_chan), mne.pick_channels_cov(data_cov, exclude=bad_chan), noise_cov=mne.pick_channels_cov(noise_cov, exclude=bad_chan), pick_ori='max-power', weight_norm='nai')
                
                # Compute SourceEstimate
                stc_vol_file = op.join(stc_path, '{}_{}_{}-{}_{}-{}-surface_{}-{}_{}-{}'.format(subject, state, '_'.join(block_group), k, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window[k], method))
                stc_vol = mne.beamformer.apply_lcmv(mne.pick_channels_evoked(evoked[k], exclude=bad_chan), filters)
                stc_vol.save(stc_vol_file)
        
        else:
            inv_vol = make_inverse_operator(mne.pick_channels_evoked(evoked[k0], exclude=bad_chan).info, mne.pick_channels_forward(fwd_vol, exclude=bad_chan), mne.pick_channels_cov(noise_cov, exclude=bad_chan), loose=1)
            
            for k in keys:
                # Compute SourceEstimate
                stc_vol_file = op.join(stc_path, '{}_{}_{}-{}_{}-{}-volume_{}-{}_{}-{}-stc.h5'.format(subject, state, '_'.join(block_group), k, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), volume, *window[k], method))
                stc_vol = apply_inverse(mne.pick_channels_evoked(evoked[k], exclude=bad_chan), inv_vol, method=method, pick_ori='vector')
                stc_vol.save(stc_vol_file)


def src_rec0(task, subject, state, block_group, evoked=dict(), noise_cov=None, keys=['R','T'], name='ECG_included', window={'R':(.35,.5), 'T':(.1,.25)}, surface='ico4', volume=6.2, bem_spacing='ico4', compute_fwd=True, compute_inv=True, compute_stc=True, fsaverage=True, baseline_cov=True, mindist=5, method="dSPM"):
    """
    If compute_fwd=True, compute and save forward solution (overwriting previously existing one).
    If compute_inv=True, compute and save inverse solution.
    If compute_stc=True, compute and save SourceEstimate.
    If fsaverage=True, morph SourceEstimate to fsaverage and save it.
    Do this for surface if provided and for volume if provided.
    Parameters:
        block_group: list of coregistered blocks of the same state to be combined
        evoked: dict of Evoked assigned to their name as returned by baseline_covariance() (will be loaded if not provided)
        noise_cov: if baseline_cov=True, should be baseline covariance, else should be empty room noise covariance (will be loaded if not provided)
        bem_spacing: BEM surface downsampling as specified for anat.BEM()
        surface: source space subdivision as specified for anat.src_space()
        volume: distance between volume sources as `pos` specified for anat.src_space()
        names: list of Epoch name
        mindist: minimum distance (in mm) of sources from inner skull surface (see mne.make_forward_solution)
        method: source reconstruction method (see mne.minimum_norm.apply_inverse)
    Output (Analyses/<task>/meg/SourceEstimate/<subject>/):
        *'-fwd.fif' (forward model)
        *'-inv.fif' (inverse model)
        *'-lh.stc' (left hemisphere SourceEstimate)
        *'-rh.stc' (right hemisphere SourceEstimate)
        *'-fsaverage-lh.stc' (left hemisphere morphed SourceEstimate)
        *'-fsaverage-rh.stc' (right hemisphere morphed SourceEstimate)
    """
    # Define pathes
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', subject)
    stc_path = op.join(Analysis_path, task, 'meg', 'SourceEstimate', subject)
    if not op.exists(stc_path):
        os.makedirs(stc_path)
    trans_file = glob.glob(op.join(Analysis_path, task, 'meg', 'Coregistration', subject, '*'+block_group[0]+'*-trans.fif'))[0]
    
    # Load source space
    bem_sol = mne.read_bem_solution(op.join(os.environ['SUBJECTS_DIR'], subject, 'bem', '{}_{}_bem_solution.fif'.format(subject, (bem_spacing if bem_spacing else 'full'))))
    if surface:
        src_surf = mne.read_source_spaces(op.join(os.environ['SUBJECTS_DIR'], subject, 'src', '{}_{}_surface-src.fif'.format(subject, surface)))
    if volume:
        src_vol = mne.read_source_spaces(op.join(os.environ['SUBJECTS_DIR'], subject, 'src', '{}_{}_volume-src.fif'.format(subject, volume)))
    
    # Load empty room noise covariance if appropriate
    if not noise_cov and not baseline_cov:
        noise_cov = mne.read_cov(op.join(cov_path, 'empty_room-cov.fif'))
    
    # Evokeds will be loaded if not provided
    load_evoked = True if not evoked else False
    
    # Initialise SourceEstimates dict
    stc_surf = dict()
    stc_vol = dict()
    fwd_surf = None
    inv_surf = None
    fwd_vol = None
    inv_vol = None

    for k in keys:
        # Load Evokeds if not provided and combine the block_group
#        if load_evoked:
#            evoked_list = []
#            for block in block_group:
#                evoked_list.append(mne.Evoked(op.join(evoked_path, '{}_{}-{}-ave.fif'.format(state, block, name)), condition = k))
        
        evoked[k] = mne.combine_evoked([mne.Evoked(op.join(evoked_path, '{}_{}-{}-ave.fif'.format(state, block, name)), condition = k) for block in block_group], 'nave')
        evoked[k].crop(*window[k])
    
    k0 = keys[0]
    
    # Bad channels will be rejected when computing the inverse operator and SourceEstimate
    bad_chan = list(set(evoked[k0].info['bads']) | set(get_chan_name(subject, 'bad', evoked[k0])))
    
    # Save files
    fwd_surf_file = op.join(stc_path, '{}_{}-{}-surface_{}-fwd.fif'.format(state, '_'.join(block_group), name, surface))
    fwd_vol_file = op.join(stc_path, '{}_{}-{}-volume_{}-fwd.fif'.format(state, '_'.join(block_group), name, volume))
    inv_surf_file = op.join(stc_path, '{}_{}-{}-{}-surface_{}-{}_{}-inv.fif'.format(state, '_'.join(block_group), name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window[k]))
    inv_vol_file = op.join(stc_path, '{}_{}-{}-{}-volume_{}-{}_{}-inv.fif'.format(state, '_'.join(block_group), name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), volume, *window[k]))
    
    # Do surface SourceEstimate
    if surface:
        # Compute forward solution
        if compute_fwd:
            fwd_surf = mne.make_forward_solution(evoked[k0].info, trans=trans_file, src=src_surf, bem=bem_sol, meg=True, eeg=False, mindist=mindist, n_jobs=4)
            mne.write_forward_solution(fwd_surf_file, fwd_surf, overwrite = True)
        
        # Compute inverse solution
        if compute_inv:
            # Load baseline covariance if appropriate
            if baseline_cov and not noise_cov:
                noise_cov = mne.read_cov(op.join(cov_path, '{}_{}_{}-{}-cov.fif'.format(subject, state, '_'.join(block_group), name)))
            
            # Load forward solution
            if not fwd_surf:
                fwd_surf = mne.read_forward_solution(fwd_surf_file)
            inv_surf = make_inverse_operator(mne.pick_channels_evoked(evoked[k0], exclude=bad_chan).info, mne.pick_channels_forward(fwd_surf, exclude=bad_chan), mne.pick_channels_cov(noise_cov, exclude=bad_chan))
            write_inverse_operator(inv_surf_file, inv_surf)
        
        for k in keys:
            # Compute SourceEstimate
            if compute_stc:
                stc_surf_file = op.join(stc_path, '{}_{}-{}_{}-{}-surface_{}-{}_{}'.format(state, '_'.join(block_group), k, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window[k]))
                # Load inverse solution
                if not inv_surf:
                    inv_surf = read_inverse_operator(inv_surf_file)
                stc_surf[k] = apply_inverse(mne.pick_channels_evoked(evoked[k], exclude=bad_chan), inv_surf, method=method)
                stc_surf[k].save(stc_surf_file)
            
            if fsaverage:
                if not stc_surf[k]:
                    stc_surf[k] = mne.read_source_estimate(stc_surf_file)
                fs_stc = mne.morph_data(subject, 'fsaverage', stc_surf[k], n_jobs=4)
                fs_stc.save(stc_surf_file + '-fsaverage')
    
    # Do volume SourceEstimate
    if volume:
        # Compute forward solution if not do_morphing else '')+'-lh.stc'
        if compute_fwd:
            fwd_vol = mne.make_forward_solution(evoked[k0].info, trans=trans_file, src=src_vol, bem=bem_sol, meg=True, eeg=False, mindist=mindist, n_jobs=4)
            mne.write_forward_solution(fwd_vol_file, fwd_vol, overwrite = True)
        
        # Compute inverse solution
        if compute_inv:
            # Load baseline covariance if appropriate
            if baseline_cov and not noise_cov:
                noise_cov = mne.read_cov(op.join(cov_path, '{}_{}_{}-{}-cov.fif'.format(subject, state, '_'.join(block_group), name)))
            
            # Load forward solution
            if not fwd_vol:
                fwd_vol = mne.read_forward_solution(fwd_vol_file)
            inv_vol = make_inverse_operator(mne.pick_channels_evoked(evoked[k0], exclude=bad_chan).info, mne.pick_channels_forward(fwd_vol, exclude=bad_chan), mne.pick_channels_cov(noise_cov, exclude=bad_chan), loose=1)
            write_inverse_operator(inv_vol_file, inv_vol)
        
        for k in keys:
            # Compute SourceEstimate
            if compute_stc:
                stc_vol_file = op.join(stc_path, '{}_{}-{}_{}-{}-volume_{}-{}_{}'.format(state, '_'.join(block_group), k, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), volume, *window[k]))
                # Load inverse solution
                if not inv_vol:
                    inv_vol = read_inverse_operator(inv_vol_file)
                stc_vol[k] = apply_inverse(mne.pick_channels_evoked(evoked[k], exclude=bad_chan), inv_vol, method=method)
                stc_vol[k].save(stc_vol_file)
    
    if compute_stc: return stc_surf,stc_vol


def fs_average(task, state, name, key, subjects=None, do_morphing=False, overwrite=False, surface='ico4', baseline_cov=True, window=(None,None)):
    """
    If do_morphing=True, morph all SourceEstimates to fsaverage, else expecting morphing to be already performed.
    For each subject, combine all SourceEstimates for the input state using the weighting of mne.evoked.combine_evoked.
    Finally, compute the grand average of all subjects.
    Parameters:
        subjects: list of subjects to include. If not provided (default), set to all subjects of the task
        overwrite: if True, combining and weighting will be performed again even if the corresponding file already exists.
        surface: source space subdivision as specified for src_rec()
        baseline_cov: if True, expecting baseline noise covariance, else empty room noise covariance
    Output (Analyses/<task>/meg/SourceEstimate/fsaverage/):
        <subject>/<state>-<name>-*-<subject>_to_fs-lh.stc (morphed and combined left hemisphere SourceEstimate)
        <subject>/*-<subject>_to_fs-rh.stc (likewise for right hemisphere)
        <state>-<name>-*-lh.stc' (grand average left hemisphere SourceEstimate)
        <state>-<name>-*-rh.stc' (likewise for right hemisphere)
    """
    mne.set_log_level(False) # Don't pollute the terminal each time a file is loaded.
    
    # Define pathes
    stc_path = op.join(Analysis_path, task, 'meg', 'SourceEstimate')
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked')
    gd_average_file = op.join(stc_path, 'fsaverage', 'Grand_average', '{}-{}_{}-{}-surface_{}-{}_{}'.format(state, key, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window))
    evo_average_file = op.join(evoked_path, 'Grand_average', '{}-{}_{}-ave.fif'.format(state, key, name))
    os.makedirs(op.dirname(gd_average_file), exist_ok=True)
    os.makedirs(op.dirname(evo_average_file), exist_ok=True)
    
    evoked_all = []
    fs_stc_all = []
    if not subjects:
        subjects = get_subjlist(task)
    print(colored('Averaging '+('and morphing ' if do_morphing else '')+key+'_'+name+' for state '+state, 'cyan'))
    
    for s,sub in enumerate(tqdm(subjects)):
        # Define files
        stc_file = op.join(stc_path, sub, '{}_{}_*-{}_{}-{}-surface_{}-{}_{}{}-lh.stc'.format(sub, state, key, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window, ('-fsaverage' if not do_morphing else '')))
        fs_file = op.join(stc_path, 'fsaverage', sub, '{}_{}-{}_{}-{}-surface_{}-{}_{}-{}_to_fs'.format(sub, state, key, name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface, *window, sub))
        os.makedirs(op.dirname(fs_file), exist_ok=True)
        
        # If the combined SourceEstimate already exists, just load it (if you do not wish to overwrite it)
        if op.isfile(fs_file) and not overwrite:
            fs_stc_all.append(mne.read_source_estimate(fs_file))
            continue
        
        files = glob.glob(stc_file) # All SourceEstimates of this subject for this state
        stc_all = []
        nave_all = []
        
        for file in files:
            if len(files) > 1: print(colored(op.basename(file[:-7]), 'yellow'))
            stc = mne.read_source_estimate(file[:-7])
            
            # Morph to fsaverage
            if do_morphing:
                stc_all.append(mne.morph_data(subject, 'fsaverage', stc, n_jobs=4, verbose=True))
            else:
                stc_all.append(stc.copy())
            
            # If there are several SourceEstimate for this state, combine them using the same weighting as mne.evoked.combine_evoked
            if len(files) > 1:
                blocks = file[file.index(state+'_'):file.index('-'+key+'_'+name)].split('_')[1:]
                nave_blk = []
                
                # To weight, retrieve the number of events (nave) from the Evoked
                evoked_list = [mne.Evoked(op.join(evoked_path, sub, '{}_{}_{}-{}-ave.fif'.format(sub, state, blk, name)), condition=key) for blk in blocks]
                nave_blk.extend([evo.nave for evo in evoked_list])
                evoked_all.append(mne.combine_evoked(evoked_list, 'nave'))
                
                # If there are several blocks corresponding to this single SourceEstimate, then the Evoked have already been combined.
                # Thus, nave should be reset according to mne.evoked.combine_evoked:
                if len(nave_blk) > 1:
                    print(colored('Re-weighting blocks {}...'.format(blocks), 'yellow'))
                    weights = np.array([n for n in nave_blk], float)
                    weights /= weights.sum()
                    nave = max(int(round(1. / sum( w ** 2 / n for w, n in zip(weights, nave_blk)))), 1)
                else:
                    nave = nave_blk[0]
                
                nave_all.append(nave)
            else:
                file = files[0]
                blk = file[file.index(state+'_'):file.index('-'+key+'_'+name)].split('_')[1]
                evoked_all.append(mne.Evoked(op.join(evoked_path, sub, '{}_{}_{}-{}-ave.fif'.format(sub, state, blk, name)), condition=key))
        
        stc = stc_all[0].copy()
        if len(files) > 1:
            # Weight and combine:
            weights = np.array([n for n in nave_blk], float)
            weights /= weights.sum()
            stc.data = sum(w * s.data for w, s in zip(weights, stc_all))
        
        # Save this subject's morphed SourceEstimate, all blocks of this state combined
        stc.save(fs_file)
        fs_stc_all.append(stc.copy())
    
    # Grand average
    evoked = mne.combine_evoked(evoked_all, 'nave')
    fs_stc = fs_stc_all[0].copy()
    data = [s.data for s in fs_stc_all]
    fs_stc.data = np.mean(data, axis=0)
    
    mne.write_evokeds(evo_average_file, evoked)
    fs_stc.save(gd_average_file)
    print(colored(gd_average_file, 'green'))
    mne.set_log_level(True)
