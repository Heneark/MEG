# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:16:49 2017

@author: oussama.abdoun
"""

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from header import *
#==============================================================================

mne.set_log_level(verbose = 'DEBUG')


task = 'SMEG'
subjects =  ['028','030','032', '037', '040', '042']


for subj in subjects:


    # PATHS 
    #==============================================================================    
    data_path = op.join(Raw_data_path, task, 'meg')
    output_path = op.join(Analysis_path, task, 'meg', 'Source_Rec', get_id(subj))
    
    if not op.exists(output_path):
        os.makedirs(output_path)
    #==============================================================================
    
            
    # BEM MODEL AND SOLUTION        
    #==============================================================================    
    # Surfaces & model      
    mne.bem.make_watershed_bem(subj, overwrite=True, show=True) # does not work with IPython
    
    model = mne.make_bem_model(subj)
    mne.write_bem_surfaces(op.join(output_path, subj + '_bem_model-bem.fif'), model)
    
    # Solution
    bem_sol = mne.make_bem_solution(model)
    mne.write_bem_solution(op.join(output_path, subj + '_bem_solution-sol.fif'), bem_sol)
    
    # Visualize results
    #mne.viz.plot_bem(subject=subj, src=None, orientation='coronal')
    #==============================================================================
    
    
    
    # SOURCE SPACE
    #==============================================================================
    spacing = 'oct6'
    
    vol = 0
    
    src = mne.setup_source_space(subj, spacing=spacing, overwrite = True)
    
    #vol_src = mne.setup_volume_source_space(subj, pos=3.0, mri = '/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/012/mri/brain.mgz', surface = '/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/012/bem/brain.surf')
    #
    #src += vol_src
    
    mne.write_source_spaces(op.join(output_path, subj +'_source_space_' + spacing + '-src.fif'), src, overwrite=True)
    #==============================================================================
    
    
    #==============================================================================
    #==============================================================================
    # # COREGISTRATION USING mne_analyze / mne.gui.coregistration() / $ mne coreg
    #==============================================================================
    #==============================================================================

