# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 16:26:17 2017

@author: oussama.abdoun
"""

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from header import *
#==============================================================================

#==============================================================================
#remove_ECG_artifacts = True
REM = [True, False]
# To remove all of them, set to 'None'
#number_of_eog_components_to_remove = None
NUM = [1, None]
#==============================================================================
task = 'SMEG'

#subjects = get_subjlist(task = task)
subjects = ['002', '004', '007', '010', '012', '014', '016', '018', '019', '021', '028', '030', '037', '032', '040', '042']

picformat = 'svg'

for subj in subjects:
    
    evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', get_id(subj))
    plots_path = op.join(Analysis_path, task, 'meg', 'Evoked', 'plots')
    ICA_path = op.join(Analysis_path, task, 'meg', 'ICA', get_id(subj))
    fnames = glob.glob(op.join(evoked_path, '*1-40_blk*-ave.fif'))
    
    for fname in fnames:
        
        for remove_ECG_artifacts in REM:
            for number_of_eog_components_to_remove in NUM:
        
        
                blk = op.split(fname)[1].strip('-ave.fif')[-2:]
                evoked = mne.read_evokeds(fname, baseline=[-.400, -.300])[0]
                
        #==============================================================================
                #load ICA file, and apply ICA
                ica_fname = sorted(glob.glob(op.join(ICA_path, get_id(subj) + '_ICA_ncomp0.99*_blk' + blk + '-ica.fif')), key = op.getmtime)[-1]
                ica = mne.preprocessing.read_ica(ica_fname)
                
                n_comp = ica.get_components().shape[1]
                exclude_inds = []
                
                t = '_withECGartifacts'
                
                if remove_ECG_artifacts == True:
                    t = ''
                    exclude_inds.extend(ica.labels_['ecg'])
                    
                if number_of_eog_components_to_remove:
                    s = '_{}EOGcomp'.format(number_of_eog_components_to_remove)
                    exclude_inds.extend(ica.labels_['eog'][:number_of_eog_components_to_remove])
                else:
                    s = '_{}EOGcomp'.format(len(ica.labels_['eog']))
                    exclude_inds.extend(ica.labels_['eog'])
                                    
                include = list(set(range(n_comp)) - set(exclude_inds))
                
                ica.apply(evoked, include = include)
        #==============================================================================        
                
                fig = evoked.plot_joint(times=(-.180,-.010,.010,.230,.380), title = op.split(fname)[1].strip('-ave.fif') + (t if t else '') + s,  ts_args=dict(ylim=dict(mag=[-150, 150])), topomap_args=dict(vmin=-50,vmax=50))
                fig.savefig(op.join(plots_path, op.split(fname)[1].strip('-ave.fif') + (t if t else '') + s + '.' + picformat), format = picformat)
                plt.close(fig)