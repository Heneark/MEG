# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:49:56 2017

@author: oussama.abdoun
"""




#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# # # # # # # # # OUT OF DATE
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================




import os, sys
import os.path as op
os.environ["FREESURFER_HOME"] = "/opt/freesurfer"
os.environ["SUBJECTS_DIR"] = "/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/"
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from io_alex import *

import mne

from mne.minimum_norm import (make_inverse_operator, apply_inverse, write_inverse_operator, read_inverse_operator)
global subj
subj = '012'
blk = '01'
state = ''
task = 'WIM'
tag = 'testWIM'
vol = 0

blocks = []
for key in get_blocks(subj, task=task):
    blocks.append(key)

# PATHS 
#==============================================================================
pathBase = '/dycog/meditation/ERC'

Analysis_path = op.join(pathBase, 'Analyses')
data_path = op.join(pathBase, 'Analyses', task, 'meg')
sr_path = op.join(data_path, 'Source_Rec', get_id(subj))
cv_path = op.join(data_path, 'Covariance', get_id(subj))
evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', get_id(subj))

#==============================================================================


# READ FORWARD SOLUTION
#==============================================================================
fname_fwd = op.join(sr_path, subj + '_' + task + ('_vol' if vol else '') + '_fwd_solution_' + blk + '-fwd.fif')
fwd = mne.read_forward_solution(fname_fwd)

# Restrict forward solution as necessary for MEG
fwd = mne.pick_types_forward(fwd, meg=True, eeg=False)
#==============================================================================


# LOAD DATA 
#==============================================================================
#blocks = ['01','02','03','04','05','06','07','08','09','10','11']
#evokeds = erp_load(subj, data_path, blocks, tag=tag)
#evoked_list = []
#for key in evokeds:
#    evoked_list.append(evokeds[key])
#
#evoked = mne.combine_evoked(evoked_list, 'nave')

evoked = mne.read_evokeds(op.join(evoked_path, get_id(subj) + '_' + tag + '_evoked-ave.fif'))[0]

#mne.Evoked.set_eeg_reference(evoked)
#==============================================================================

# LOAD NOISE COVARIANCE FROM EPOCHS
#==============================================================================
#noise_covs = cov_load(subj, data_path, blocks, tag=tag)
#noise_cov = noise_covs[blk]

noise_cov = mne.read_cov('/dycog/meditation/ERC/Analyses/SMEG/meg/Covariance/012/012_SMEG_raw_covariance-cov.fif')
#==============================================================================

# Make an MEG inverse operator
#==============================================================================
info = evoked.info
inverse_operator = make_inverse_operator(info, fwd, noise_cov, verbose='DEBUG')

write_inverse_operator(op.join(sr_path, subj + '_' + task + '_inverse_operator_mne-inv.fif'), inverse_operator)

#inverse_operator = read_inverse_operator(op.join(sr_path, subj + '_' + task + '_inverse_operator_mne-inv.fif'))
src = inverse_operator['src']
#==============================================================================



##==============================================================================
#stc = lcmv(evoked, fwd, noise_cov, noise_cov, reg=0.05, pick_ori=None)
#
## Save result in stc files
#stc.save('lcmv-vol')
#
#stc.crop(0.0, 0.2)
#
## Save result in a 4D nifti file
#img = mne.save_stc_as_volume('lcmv_inverse.nii.gz', stc,
#                             forward['src'], mri_resolution=False)
#
#t1_fname = data_path + '/subjects/sample/mri/T1.mgz'
#
## Plotting with nilearn ######################################################
#plot_stat_map(index_img(img, 61), t1_fname, threshold=0.8,
#              title='LCMV (t=%.1f s.)' % stc.times[61])
#
## plot source time courses with the maximum peak amplitudes
#plt.figure()
#plt.plot(stc.times, stc.data[np.argsort(np.max(stc.data, axis=1))[-40:]].T)
#plt.xlabel('Time (ms)')
#plt.ylabel('LCMV value')
#plt.show()
##==============================================================================





# Compute inverse solution
#==============================================================================
method = "dSPM"
snr = 3.
lambda2 = 1. / snr ** 2
stc = apply_inverse(evoked, inverse_operator, lambda2, method=method, pick_ori=None)
stc.save(op.join(sr_path, subj + '_' + task + '_source_estimate'))
stc.crop(0.325,0.7)
stc.save(op.join(sr_path, subj + '_' + task + '_source_estimate_cropped'))

#==============================================================================


#plt.plot(1e3 * stc.times, stc.data[::100, :].T)
#plt.xlabel('time (ms)')
#plt.ylabel('%s value' % method)
#plt.show()
#
#
#vertno_max, time_max = stc.get_peak(hemi='rh')
#
#subjects_dir = "/dycog/meditation/ERC/Analyses/ANAT/T1/FreeSurfer/"
#brain = stc.plot(surface='inflated', hemi='rh', subjects_dir=subjects_dir, time_unit='s')
#brain.add_foci(vertno_max, coords_as_verts=True, hemi='rh', color='blue', scale_factor=0.6)
#brain.show_view('lateral')