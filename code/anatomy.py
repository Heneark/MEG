# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:16:49 2017

@author: oussama.abdoun
"""

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/MEG/code/")
import matplotlib
matplotlib.use('Agg') #not to show BEM surfaces and save them instead
from header import *
#==============================================================================

#mne.set_log_level(verbose = 'DEBUG')


def BEM(subject, spacing='ico4', watershed=True):
    """
    Creates bem/watershed/ with surface files in FreeSurfer's SUBJECTS_DIR directory.
    Parameters:
        spacing: BEM surface downsampling = 'ico3', 'ico4' (default), 'ico5', or None for full model (no subsampling).
    Output:
        SUBJECTS_DIR/<subject>/bem/<subject>_<spacing>_bem_solution.fif
    """
    
    output_path = op.join(os.environ['SUBJECTS_DIR'], subject, 'bem')
    
    # BEM surfaces
    if watershed:
        mne.bem.make_watershed_bem(subject, overwrite=True, show=True) # does not work with IPython
        plt.savefig(op.join(output_path,'BEM_surfaces.png'))
        plt.close()
    
    # BEM model
    surfaces = mne.make_bem_model(subject, ico = (int(spacing[3]) if spacing else None))
    
    # BEM solution
    bem_sol = mne.make_bem_solution(surfaces)
    mne.write_bem_solution(op.join(output_path, '{sub}_{ico}_bem_solution.fif'.format(sub=subject, ico=(spacing if spacing else 'full'))), bem_sol)
    
    # Visualize results
    #mne.viz.plot_bem(subject=subject, src=None, orientation='coronal')


def src_space(subject, spacing='ico4', pos=6.2):
    """
    Parameters:
        spacing: 'oct5', 'ico4' (default), 'oct6' (default of mne.setup_source_space), 'ico5', or 'all' (no subdivision of the source space).
        pos: distance (in mm) between sources (arranged as a grid).
    Output:
        SUBJECTS_DIR/<subject>/src/<subject>_<spacing>_surface-src.fif
        SUBJECTS_DIR/<subject>/src/<subject>_<pos>_volume-src.fif
        SUBJECTS_DIR/<subject>/src/<subject>_<spacing>_<pos>_mixed-src.fif
    
    MNE naming conventions: all source space files should end with -src.fif or -src.fif.gz.
    """
    
    output_path = op.join(os.environ['SUBJECTS_DIR'], subject, 'src')
    mri_path = op.join(os.environ['SUBJECTS_DIR'], subject, 'mri', 'brain.mgz')
    surf_path = op.join(os.environ['SUBJECTS_DIR'], subject, 'bem', 'brain.surf')
    
    if not op.exists(output_path):
        os.makedirs(output_path)
    
    # Surface source space
    src_surf = mne.setup_source_space(subject, spacing=spacing, n_jobs=2)
    src_surf.subject = subject
    mne.write_source_spaces(op.join(output_path, '{}_{}_surface-src.fif'.format(subject, spacing)), src_surf, overwrite=True)
    
    # Volume source space
    src_vol = mne.setup_volume_source_space(subject, pos=pos, mri=mri_path, surface=surf_path)
    src_vol.subject = subject
    mne.write_source_spaces(op.join(output_path, '{}_{}_volume-src.fif'.format(subject, pos)), src_vol, overwrite=True)
    
    # Mixed source space
#    src_mix = src_surf + src_vol
#    src_mix.subject = subject
#    mne.write_source_spaces(op.join(output_path, '{}_{}_{}_mixed-src.fif'.format(subject, spacing, pos)), src_mix, overwrite=True)
