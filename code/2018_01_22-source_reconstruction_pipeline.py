# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:30:25 2018

@author: oussama.abdoun
"""

#IMPORT PACKAGES, DEFINE ENVIRONMENT VARIABLES, CREATE PATHES
#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/MEG/code/")
import matplotlib
matplotlib.use('Agg') #not to display figures and run all subjects at once
from header import *
#==============================================================================


#PARAMETERS
#==============================================================================
task = 'SMEG'
states = ['RS','FA','OM']
subjects = get_subjlist(task)
subjects = subjects[subjects.index('089'):]
# # Last subject preprocessed: 093
# # Future subjects list:
#subjects = subjects[subjects.index('093')+1:]
subjects.sort()
#==============================================================================
#subjects = subjects[subjects.index('068')+1:]
#subjects =  ['028','030','032', '037', '040', '042']


#PROCESSING SCRIPTS
#==============================================================================
import preprocessing
import anatomy
import HPI_update
import covariance_matrix
import source_reconstruction
import morphing
#==============================================================================


#preprocessing.process(tasks=[task], states=states, subjects=subjects, run_ICA=False)

# ANATOMICAL RECONSTRUCTION: FREESURFER
#==============================================================================
# recon-all -i <sub_T1.nii> -s <sub> -all
# # INPUT: MRI T1 raw data
# # OUPUT in FreeSurfer SUBJECTS_DIR
#==============================================================================

#anatomy.BEM(task=task, subjects=subjects)

#COREGISTRATION: MATLAB
#==============================================================================
# Adjust_Head_Pos_3DirbyCoil.m
# # Select which block to coregister --> new .hc files
#==============================================================================

#HPI_update.update(task=task, subjects=subjects)

#%gui wx
#mne.gui.coregistration()
# https://www.slideshare.net/mne-python/mnepython-coregistration
# Files corresponding to the coregistration: *'-trans.fif' (Source_Rec directory)

#covariance_matrix.covariance(task=task, states=states, subjects=subjects, do_epochs_cov=True)

#source_reconstruction.src_rec(task=task, states=states, subjects=subjects)

#morphing.average(task=task, subjects=subjects)

