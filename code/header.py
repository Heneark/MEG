# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 18:02:35 2017

@author: oussama.abdoun
"""

import autoreject
from autoreject import get_rejection_threshold
import sys, os
import os.path as op
import glob
import inspect
import joblib# allows for parallel processing
import locale
import mpl_toolkits.axes_grid1
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from numpy.linalg import inv
import pandas as pd
import scipy
import sklearn# improves performance for e.g. anatomy.py
from termcolor import colored
import time
from tqdm import tqdm
import warnings
import xarray as xr

#bugfix with mne.io.read_ctf
try:
    locale.setlocale(locale.LC_ALL, "en_US.utf8")
except:
    pass


import mne
from mne.time_frequency import psd_welch
from mne.preprocessing import ICA
from mne.preprocessing import create_eog_epochs, create_ecg_epochs, read_ica
from mne.event import make_fixed_length_events
from mne.minimum_norm import make_inverse_operator, apply_inverse, write_inverse_operator, read_inverse_operator
from mne.report import Report


# ENVIRONMENT VARIABLES
#==============================================================================
os.environ["ETS_TOOLKIT"] = "wx"
os.environ["FREESURFER_HOME"] = "/opt/freesurfer"
os.environ["SUBJECTS_DIR"] = "/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/"
#==============================================================================


# CUSTOM SCRIPTS WITH PATH BASES
#==============================================================================
#from io_alex import *
from csv_io import *
#==============================================================================
