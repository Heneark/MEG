# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 11:20:52 2017

@author: oussama.abdoun
"""




#==============================================================================
#==============================================================================
#==============================================================================
# # # A REVOIR OU SUPPRIMER ...
#==============================================================================
#==============================================================================
#==============================================================================






import mne.gui
import os.path as op
import locale, os
import matplotlib

import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from csv_io_alex import *

matplotlib.use('Qt4Agg')
matplotlib.interactive(True)


os.environ["FREESURFER_HOME"] = "/opt/freesurfer"
os.environ["SUBJECTS_DIR"] = "/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/"


blk = '01'
subj = '012'
subj_dir = '/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/'
data_path = '/dycog/meditation/ERC/Raw data/SMEG/meg/'

ex_path = op.join('/home/oussama.abdoun/Documents/donneesMEG/', get_rawpath(subj) + blk + '_raw.fif')

#locale.setlocale(locale.LC_ALL, "en_US.utf8")
        
#file name
#raw_fname = data_path + get_rawpath(subj) + blk + '.ds'
        
mne.gui.coregistration(inst=ex_path, subject=subj)