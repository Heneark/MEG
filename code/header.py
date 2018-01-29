# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 18:02:35 2017

@author: oussama.abdoun
"""

import sys, os
import os.path as op
import glob
import joblib
import locale
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from numpy.linalg import inv
import pandas as pd
import scipy
import sklearn
from termcolor import colored
import time

#bugfix with mne.io.read_ctf
locale.setlocale(locale.LC_ALL, "en_US.utf8")


import mne
from mne.time_frequency import psd_welch
from mne.preprocessing import ICA
from mne.preprocessing import create_eog_epochs, create_ecg_epochs
from mne.event import make_fixed_length_events
from mne.minimum_norm import (make_inverse_operator, apply_inverse, write_inverse_operator, read_inverse_operator)


sys.path.append("/dycog/meditation/ERC/Analyses/MEG/code/")
from io_alex import *
from csv_io_alex import *


#ENVIRONMENT VARIABLES
#==============================================================================
os.environ["ETS_TOOLKIT"] = "wx"
os.environ["FREESURFER_HOME"] = "/opt/freesurfer"
os.environ["SUBJECTS_DIR"] = "/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/"
#==============================================================================


#PATH BASES
#==============================================================================
pathBase = '/dycog/meditation/ERC'
Analysis_path = op.join(pathBase, 'Analyses')
Raw_data_path = op.join(pathBase, 'Raw data')
#==============================================================================
