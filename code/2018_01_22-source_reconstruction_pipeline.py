# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:30:25 2018

@author: oussama.abdoun
"""

#IMPORT PACKAGES, DEFINE ENVIRONMENT VARIABLES, CREATE PATHES
#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/MEG/code/")
from header import *
#==============================================================================


#PARAMETERS
#==============================================================================
task = 'SMEG'
states = ['RS','FA','OM']
subjects = get_subjlist(task)
#subjects =  ['028','030','032', '037', '040', '042']
subjects = subjects[subjects.index('089'):]
subjects.sort()
#==============================================================================


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

