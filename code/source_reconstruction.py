# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 17:42:45 2017

@author: oussama.abdoun
"""

#==============================================================================
from header import *
#==============================================================================
warnings.filterwarnings("ignore",category=DeprecationWarning)


def baseline_covariance(task, subject, state, block_group, baseline=(-.4,-.25), t_delay=.3, names=['R_ECG_included','R_ECG_excluded','T_ECG_included','T_ECG_excluded']):
    """
    Returns baseline noise covariance for the block_group and a dict containing a list of Evoked (for each block in block_group) assigned to their name.
    Parameters:
        block_group: list of coregistered blocks of the same state to be combined
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
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    if not op.exists(evoked_path):
        os.makedirs(evoked_path)
    
    # Initialise Epochs and Evokeds dict
    epochs = dict()
    evoked = dict()
    
    for name in names:
        # Use the proper baseline according the the Epochs name
        if not baseline:
            use_baseline = (None,0)
        elif 'R_ECG' in name:
            use_baseline = baseline
        elif 'T_ECG' in name:
            use_baseline = (baseline[0]-t_delay, baseline[1]-t_delay)
        else:
            use_baseline = baseline
        
        epochs[name] = []
        evoked[name] = []
        
        for b,block in enumerate(block_group):
            # Apply baseline
            epochs[name].append(mne.read_epochs(op.join(Analysis_path, task, 'meg', 'Epochs', subject, '{}-{}_{}-epo.fif'.format(name, state, block))))
            epochs[name][b].apply_baseline(use_baseline)
            
            # Save Evoked
            evoked[name].append(epochs[name][b].average())
            evoked[name][b].save(op.join(evoked_path, '{}-{}_{}-ave.fif'.format(name, state, block)))
        
        # Compute and save noise covariance
        epochs[name] = mne.concatenate_epochs(epochs[name])
        noise_cov = mne.compute_covariance(epochs[name], n_jobs=4)
        
        cov_file = op.join(cov_path, '{}-{}_{}-cov.fif'.format(name, state, '_'.join(block_group)))
        mne.write_cov(cov_file, noise_cov)
    
    return noise_cov,evoked


def src_rec(task, subject, state, block_group, evoked=None, noise_cov=None, surface='ico4', volume=6.2, bem_spacing='ico4', compute_fwd=True, compute_inv=True, compute_stc=True, baseline_cov=True, names=['R_ECG_included','R_ECG_excluded','T_ECG_included','T_ECG_excluded'], mindist=5, method="dSPM"):
    """
    If compute_fwd=True, compute and save forward solution (overwriting previously existing one).
    If compute_inv=True, compute and save inverse solution.
    If compute_stc=True, compute and save SourceEstimate.
    Do this for surface if provided and for volume if provided.
    Parameters:
        block_group: list of coregistered blocks of the same state to be combined
        evoked: dict of Evoked assigned to their name as returned by baseline_covariance() (will be loaded if not provided)
        noise_cov: if baseline_cov=True, should be empty room noise covariance, else should be baseline covariance (will be loaded if not provided)
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
    
    # Initialise Evokeds dict if not provided
    load_evoked = False
    if not evoked:
        evoked = {name:[] for name in names}
        load_evoked = True
    
    # Initialise SourceEstimates dict
    stc_surf = {name:[] for name in names}
    stc_vol = {name:[] for name in names}
    
    for name in names:
        # Load Evokeds if not provided and combine the block_group
        for block in block_group:
            if load_evoked:
                evoked[name].extend(mne.read_evokeds(op.join(evoked_path, '{}-{}_{}-ave.fif'.format(name, state, block))))
        evoked[name] = mne.combine_evoked(evoked[name], 'nave')
        
        # Save files
        fwd_surf_file = op.join(stc_path, '{}_{}-{}-surface_{}-fwd.fif'.format(state, '_'.join(block_group), name, surface))
        fwd_vol_file = op.join(stc_path, '{}_{}-{}-volume_{}-fwd.fif'.format(state, '_'.join(block_group), name, volume))
        inv_surf_file = op.join(stc_path, '{}_{}-{}-{}-surface_{}-inv.fif'.format(state, '_'.join(block_group), name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface))
        inv_vol_file = op.join(stc_path, '{}_{}-{}-{}-volume_{}-inv.fif'.format(state, '_'.join(block_group), name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), volume))
        stc_surf_file = op.join(stc_path, '{}_{}-{}-{}-surface_{}'.format(state, '_'.join(block_group), name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface))
        stc_vol_file = op.join(stc_path, '{}_{}-{}-{}-volume_{}'.format(state, '_'.join(block_group), name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), volume))
        
        # Do surface SourceEstimate
        if surface:
            # Compute forward solution
            fwd_surf = None
            if compute_fwd:
                fwd_surf = mne.make_forward_solution(evoked[name].info, trans=trans_file, src=src_surf, bem=bem_sol, meg=True, eeg=False, mindist=mindist, n_jobs=4)
                mne.write_forward_solution(fwd_surf_file, fwd_surf, overwrite = True)
            
            # Compute inverse solution
            inv_surf = None
            if compute_inv:
                # Load baseline covariance if appropriate
                if baseline_cov and not noise_cov:
                    noise_cov = mne.read_cov(op.join(cov_path, '{}-{}_{}-cov.fif'.format(name, state, '_'.join(block_group))))
                
                # Load forward solution
                if not fwd_surf:
                    fwd_surf = mne.read_forward_solution(fwd_surf_file)
                inv_surf = make_inverse_operator(evoked[name].info, fwd_surf, noise_cov)
                write_inverse_operator(inv_surf_file, inv_surf)
            
            # Compute SourceEstimate
            if compute_stc:
                # Load inverse solution
                if not inv_surf:
                    inv_surf = read_inverse_operator(inv_surf_file)
                stc_surf[name] = apply_inverse(evoked[name], inv_surf, method=method)
                stc_surf[name].subject = subject
                stc_surf[name].save(stc_surf_file)
        
        # Do volume SourceEstimate
        if volume:
            # Compute forward solution
            fwd_vol = None
            if compute_fwd:
                fwd_vol = mne.make_forward_solution(evoked[name].info, trans=trans_file, src=src_vol, bem=bem_sol, meg=True, eeg=False, mindist=mindist, n_jobs=4)
                mne.write_forward_solution(fwd_vol_file, fwd_vol, overwrite = True)
            
            # Compute inverse solution
            inv_vol = None
            if compute_inv:
                # Load baseline covariance if appropriate
                if baseline_cov and not noise_cov:
                    noise_cov = mne.read_cov(op.join(cov_path, '{}-{}_{}-cov.fif'.format(name, state, '_'.join(block_group))))
                
                # Load forward solution
                if not fwd_vol:
                    fwd_vol = mne.read_forward_solution(fwd_vol_file)
                inv_vol = make_inverse_operator(evoked[name].info, fwd_vol, noise_cov, loose=1)
                write_inverse_operator(inv_vol_file, inv_vol)
            
            # Compute SourceEstimate
            if compute_stc:
                # Load inverse solution
                if not inv_vol:
                    inv_vol = read_inverse_operator(inv_vol_file)
                stc_vol[name] = apply_inverse(evoked[name], inv_vol, method=method)
                stc_vol[name].subject = subject
                stc_vol[name].save(stc_vol_file)
    
    return stc_surf,stc_vol


def fs_average(task, subject, state, block_group, stc=None, surface='ico4', names=['R_ECG_included','R_ECG_excluded','T_ECG_included','T_ECG_excluded'], baseline_cov=True):
    """
    Morph surface SourceEstimate to fsaverage (does not work with volume).
    Parameters:
        block_group: list of coregistered blocks of the same state to be combined
        surface: source space subdivision as specified for src_rec()
        names: list of Epoch name
    Output (Analyses/<task>/meg/SourceEstimate/<subject>/):
        *'-fsaverage-lh.stc' (left hemisphere SourceEstimate)
        *'-fsaverage-rh.stc' (right hemisphere SourceEstimate)
    """
    # Define pathes
    stc_path = op.join(Analysis_path, task, 'meg', 'SourceEstimate', subject)
    
    # Initialise SourceEstimates dict if not provided
    load_stc = False
    if not stc:
        stc = {name:[] for name in names}
        load_stc = True
    
    for name in names:
        stc_file = op.join(stc_path, '{}_{}-{}-{}-surface_{}'.format(state, '_'.join(block_group), name, ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface))
        
        # Load SourceEstimates if not provided
        if load_stc:
            stc[name] = mne.read_source_estimate(stc_file)
        
        # Morph and save SourceEstimate
        fs_stc = mne.morph_data(subject, 'fsaverage', stc, n_jobs=4)
        fs_stc.subject = subject + '_fsaverage
        fs_stc.save(stc_file + '-fsaverage')
