#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 19:04:54 2018

@author: benjamin.ador
"""

import matplotlib
matplotlib.use('Qt5Agg')
from header import *
import visbrain
from visbrain import Brain
from visbrain.objects import BrainObj

#%% PARAMETERS
task = 'SMEG'
names = ['T_ECG_included', 'T_ECG_excluded']#'R_ECG_included', 'R_ECG_excluded', 
name = names[1]
noise_cov = 'empty_room_cov'
surface = 'ico4'
sfreq = 200
window = 'baseline'

stat_path = op.join(Analysis_path, task, 'meg', 'Stats', window)
stat_file = op.join(stat_path, '{}-{}-surface_{}-{}Hz-{}.nc'.format(name, noise_cov, surface, sfreq, window))
test_key = 'RS'
p_thresh = 0.05


#%% DATA
stats = xr.open_dataarray(stat_file, test_key)
stats.load()

t_clu_size = np.unique(np.where(stats.loc['p_val'].values < p_thresh)[1], return_counts=True)
t = t_clu_size[0][t_clu_size[1].argmax()] if t_clu_size[0].size else 0
t=19

time = stats.time.values[t]
data = np.where(stats.loc['p_val', :, time].values < p_thresh, stats.loc['T_stat', :, time].values, -np.inf)
cmax = np.max(np.abs(stats.loc['T_stat', :, time].values))
clim = (-cmax, cmax)


#%% PLOT HEMISPHERES SEPARATELY
brains = dict()
#os.makedirs(op.splitext(stat_file)[0], exist_ok=True)
#for h,hemi in enumerate(stats.hemisphere.values):
#    b_obj = BrainObj('inflated', translucent=False, hemisphere=hemi, sulcus=True)
#    b_obj.add_activation(data=data[h], vertices=stats.src.values, hemisphere=hemi, clim=clim, hide_under=clim[0], cmap='cool', vmin=clim[0], vmax=clim[1], under=None)
#    brains[hemi] = Brain(brain_obj=b_obj, bgcolor=None)
    
#    views = {'R': (110,0), 'L': (250,0)}#, 'F': (180,-45), 'B': (0,45)}
#    for view, rotation in views.items():
#        brains[hemi].rotate(custom=rotation)
#        vb.screenshot(op.join(op.splitext(stat_file)[0], '{}_{}-{}_{}.png'.format(test_key, round(time, 3), hemi, view)), print_size=(5,3), dpi=360, autocrop=True, transparent=True)
#
#brains[hemi].cbar_select('brain')
#brains[hemi].cbar_control('brain', cblabel='T statistic', bgcolor=None, txtcolor='black')
#vb.screenshot(op.join(op.splitext(stat_file)[0], '{}_{}-colorbar.png'.format(test_key, round(time, 3))), canvas='colorbar', dpi=360, transparent=True)#


#%% SHOW WHOLE BRAIN
b_obj = BrainObj('inflated', translucent=False, hemisphere='both', sulcus=True)
for h,hemi in enumerate(stats.hemisphere.values):
    b_obj.add_activation(data=data[h], vertices=stats.src.values, hemisphere=hemi, clim=clim, hide_under=clim[0], cmap='cool', vmin=clim[0], vmax=clim[1], under=None)
#    b_obj.add_activation(data=stats.loc['T_stat', :, time].values[h], vertices=stats.src.values, hemisphere=hemi, clim=clim, hide_under=clim[0], cmap='cool', vmin=clim[0], vmax=clim[1], under=None)
hemi = 'both'
brains[hemi] = Brain(brain_obj=b_obj, bgcolor=None)

#views = {'F': (180,-45), 'B': (0,45)}
#for view, rotation in views.items():
#    brains[hemi].rotate(custom=rotation)
#    vb.screenshot(op.join(op.splitext(stat_file)[0], '{}_{}-both_{}.png'.format(test_key, round(time, 3), view)), print_size=(5,3), dpi=360, autocrop=True, transparent=True)

brains[hemi].cbar_select('brain')
brains[hemi].cbar_control('brain', cblabel='T statistic', bgcolor=None, txtcolor='black')
brains[hemi].show()
