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
    
    """
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', subject)
    if not op.exists(cov_path):
        os.makedirs(cov_path)
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    if not op.exists(evoked_path):
        os.makedirs(evoked_path)
    
    epochs = dict.fromkeys(names, [])
    evoked = dict.fromkeys(names, [])
    
    for name in names:
        if not baseline:
            use_baseline = (None,0)
        elif 'R_ECG' in name:
            use_baseline = baseline
        elif 'T_ECG' in name:
            use_baseline = (baseline[0]-t_delay, baseline[1]-t_delay)
        else:
            use_baseline = baseline
        
        for b,block in enumerate(block_group):
            epochs[name].append(mne.read_epochs(op.join(Analysis_path, task, 'meg', 'Epochs', subject, '{}-{}_{}-epo.fif'.format(name, state, block))))
            epochs[name][b].apply_baseline(use_baseline)
            
            evoked[name].append(epochs[name][b].average())
            evoked[name][b].save(op.join(evoked_path, '{}-{}_{}-ave.fif'.format(name, state, block)))
        
        epochs[name] = mne.concatenate_epochs(epochs[name])
        
        noise_cov = mne.compute_covariance(epochs[name], n_jobs=4)
        cov_file = op.join(cov_path, '{}-{}_{}-cov.fif'.format(name, state, '_'.join(block_group)))
        mne.write_cov(cov_file, noise_cov)
    
    return noise_cov,evoked


def src_rec(task, subject, state, block_group, evoked=None, noise_cov=None, surface='ico4', volume=6.2, bem_spacing='ico4', compute_fwd=True, compute_inv=True, compute_stc=True, baseline_cov=True, names=['R_ECG_included','R_ECG_excluded','T_ECG_included','T_ECG_excluded'], method="dSPM"):
    """
    Output (Source_Rec directory): *'-fwd.fif' (forward model), *'-inv.fif' (inverse model), *'-rh.stc' (right hemisphere source estimates), and *'-lh.stc' (left hemisphere).
    Parameters:
        states: list of the states.
        subjects: list of the subjects, default to all subjects available for the task.
        spacing: 'oct5', 'ico4', 'oct6', or 'ico5' according to mne source_space. Should be as specified for anatomy.BEM().
        method: source reconstruction method.
    """
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
    cov_path = op.join(Analysis_path, task, 'meg', 'Covariance', subject)
    stc_path = op.join(Analysis_path, task, 'meg', 'SourceEstimate', subject)
    if not op.exists(stc_path):
        os.makedirs(stc_path)
    trans_file = glob.glob(op.join(Analysis_path, task, 'meg', 'Coregistration', '*'+block_group[0]+'*-trans.fif'))[0]
    
    if not noise_cov and not baseline_cov:
        noise_cov = mne.read_cov(op.join(cov_path, 'empty_room-cov.fif'))
    
    bem_sol = mne.read_bem_solution(op.join(os.environ['SUBJECTS_DIR'], subject, 'bem', '{}_{}_bem_solution.fif'.format(subject, bem_spacing)))
    if surface:
        src_surf = mne.read_source_spaces(op.join(os.environ['SUBJECTS_DIR'], subject, 'src', '{}_{}_surface-src.fif'.format(subject, surface)))
    if volume:
        src_vol = mne.read_source_spaces(op.join(os.environ['SUBJECTS_DIR'], subject, 'src', '{}_{}_volume-src.fif'.format(subject, volume)))
    
    load_evoked = False
    if not evoked:
        evoked = dict.fromkeys(names, [])
        load_evoked = True
    
    stc_surf = dict.fromkeys(names, [])
    stc_vol = dict.fromkeys(names, [])
    
    for n,name in enumerate(names):
        for b,block in enumerate(block_group):
            if load_evoked:
                evoked[name].extend(mne.read_evokeds(op.join(evoked_path, '{}-{}_{}-epo.fif'.format(name, state, block))))
        evoked[name] = mne.combine_evoked(evoked[name], 'nave')
        
        if surface:
            fwd_surf = None
            if compute_fwd and not n:
                fwd_surf = mne.make_forward_solution(evoked[name].info, trans=trans_file, src=src_surf, bem=bem_sol, meg=True, eeg=False, mindist=5.0, n_jobs=4)
                mne.write_forward_solution(op.join(stc_path, '{}_{}-{}_surface-fwd.fif'.format(state, '_'.join(block_group), surface)), fwd_surf, overwrite = True)
            
            inv_surf = None
            if compute_inv and not n:
                if baseline_cov and not noise_cov:
                    noise_cov = mne.read_cov(op.join(cov_path, '{}-{}_{}-cov.fif'.format(name, state, '_'.join(block_group))))
                
                if not fwd_surf:
                    fwd_surf = mne.read_forward_solution(op.join(stc_path, '{}_{}-{}_surface-fwd.fif'.format(state, '_'.join(block_group), surface)))
                inv_surf = make_inverse_operator(evoked[name].info, fwd_surf, noise_cov)
                write_inverse_operator(op.join(stc_path, '{}_{}-{}-{}_surface-inv.fif'.format(state, '_'.join(block_group), ('baseline' if baseline_cov else 'empty_room'), surface)), inv_surf)
            
            if compute_stc:
                if not inv_surf:
                    inv_surf = read_inverse_operator(op.join(stc_path, '{}_{}-{}-{}_surface-inv.fif'.format(state, '_'.join(block_group), ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface)))
                stc_surf[name] = apply_inverse(evoked[name], inv_surf, method=method)
                stc_surf[name].save(op.join(stc_path, '{}-{}_{}-{}-{}_surface'.format(name, state, '_'.join(block_group), ('baseline_cov' if baseline_cov else 'empty_room_cov'), surface)))
        
        if volume:
            fwd_vol = None
            if compute_fwd and not n:
                fwd_vol = mne.make_forward_solution(evoked[name].info, trans=trans_file, src=src_vol, bem=bem_sol, meg=True, eeg=False, mindist=5.0, n_jobs=4)
                mne.write_forward_solution(op.join(stc_path, '{}_{}-{}_volume-fwd.fif'.format(state, '_'.join(block_group), volume)), fwd_vol, overwrite = True)
            
            inv_vol = None
            if compute_inv and not n:
                if baseline_cov and not noise_cov:
                    noise_cov = mne.read_cov(op.join(cov_path, '{}-{}_{}-cov.fif'.format(name, state, '_'.join(block_group))))
                
                if not fwd_vol:
                    fwd_vol = mne.read_forward_solution(op.join(stc_path, '{}_{}-{}_volume-fwd.fif'.format(state, '_'.join(block_group), volume)))
                inv_vol = make_inverse_operator(evoked[name].info, fwd_vol, noise_cov)
                write_inverse_operator(op.join(stc_path, '{}_{}-{}-{}_volume-inv.fif'.format(state, '_'.join(block_group), ('baseline' if baseline_cov else 'empty_room'), volume)), inv_vol)
            
            if compute_stc:
                if not inv_vol:
                    inv_vol = read_inverse_operator(op.join(stc_path, '{}_{}-{}-{}_volume-inv.fif'.format(state, '_'.join(block_group), ('baseline_cov' if baseline_cov else 'empty_room_cov'), volume)))
                stc_vol[name] = apply_inverse(evoked[name], inv_vol, method=method)
                stc_vol[name].save(op.join(stc_path, '{}-{}_{}-{}-{}_volume'.format(name, state, '_'.join(block_group), ('baseline_cov' if baseline_cov else 'empty_room_cov'), volume)))
    
    return stc_surf,stc_vol
