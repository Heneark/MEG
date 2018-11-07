# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 16:54:53 2017

@author: oussama.abdoun
"""

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from header import *
#==============================================================================


#==============================================================================
#==============================================================================
#==============================================================================
# # # USE WITH "%gui wx"
#==============================================================================
#==============================================================================
#==============================================================================


task = 'SMEG'
states = ['RS','OM','FA']
subjects = ['002', '007', '012', '016', '028', '030', '032', '037', '040', '042']
#subjects = ['002', '007', '012', '016', '028', '032', '037', '042']

#==============================================================================
tag = 'Cardiac'

ecgchoice = 'ECGinclude'
#ecgchoice = 'ECGexclude'

#Crop the stcs
tmin = .350
tmax = .500
#==============================================================================

tstep = None
full_data_array = None
S = 0

for subj in subjects:
    
    ST = 0    
    
    sr_path = op.join(Analysis_path, task, 'meg', 'Source_Rec', get_id(subj))
    stc_path = op.join(sr_path, 'STC')
    stc_paths = glob.glob(op.join(stc_path, '*' + ecgchoice + '*fsaverage-rh.stc'))
    stc_blocks = []

    if full_data_array is None:
        full_data_array = np.empty((len(subjects), len(states), mne.read_source_estimate(stc_paths[0]).crop(tmin, tmax).data.shape[1], mne.read_source_estimate(stc_paths[0]).data.shape[0]))
    
    if tstep is None:
        tstep = mne.read_source_estimate(stc_paths[0]).crop(tmin, tmax).tstep
    
    for path in stc_paths:
        #Careful with file names : classic [4:], withECGartifacts [7:]
        # Need to find a better way to do it ...
        stc_blocks.append(op.split(path)[1].strip('_fsaverage-rh.stc').split('blk_')[1].split('_')) 

    
    for state in states:
        
        data_list = []        
        
        blocks = []
        
        for key in get_blocks(subj, task = task, state = state):
            blocks.append(key)
        
        for i in range(len(stc_blocks)):
            b = list(set(stc_blocks[i]) & set(blocks))
            b.sort()
            
            if b:
                data_list.append(np.expand_dims(mne.read_source_estimate(stc_paths[i]).crop(tmin, tmax).data.transpose(), axis=2))
        #PONDERATION A REVOIR           
        stc_state_data = np.concatenate(data_list, axis=2).mean(axis=2)         
                   
                   
                   
        full_data_array[S, ST, :, :] = stc_state_data
                   
        ST += 1
    S += 1
    

X = full_data_array[:,1,:,:] - full_data_array[:,0,:,:]

#X = full_data_array[:,1,:,:]

print('Computing connectivity.')
connectivity = mne.spatial_tris_connectivity(mne.grade_to_tris(5))

n_permutations = 10000

p_threshold = 0.005
n_subjects = full_data_array.shape[0]
t_threshold = -scipy.stats.distributions.t.ppf(p_threshold / 2., n_subjects - 1)

print('Clustering.')
T_obs, clusters, cluster_p_values, H0 = clu = mne.stats.spatio_temporal_cluster_1samp_test(X, connectivity=connectivity, n_jobs=6, threshold=t_threshold, n_permutations = n_permutations)
#    Now select the clusters that are sig. at p < 0.05 (note that this value
#    is multiple-comparisons corrected).
good_cluster_inds = np.where(cluster_p_values < 0.05)[0]



print('Visualizing clusters.')

#    Now let's build a convenient representation of each cluster, where each
#    cluster becomes a "time point" in the SourceEstimate
stc_all_cluster_vis = mne.stats.summarize_clusters_stc(clu, p_thresh=0.05, tstep=tstep, vertices= [np.arange(10242), np.arange(10242)], subject='fsaverage')

#    Let's actually plot the first "time point" in the SourceEstimate, which
#    shows all the clusters, weighted by duration

# blue blobs are for condition A < condition B, red for A > B
brain = stc_all_cluster_vis.plot(hemi='split', views='lateral', time_label='Duration significant (ms)', colormap = 'mne')

#brain.save_image(op.join('/dycog/meditation/ERC/Analyses/SMEG/meg/Source_Rec/stats','spatio_temporal_cluster_1samp_10000_0.001_OM-RS.png'))
#brain.save_image(op.join('/dycog/meditation/ERC/Analyses/SMEG/meg/Source_Rec/stats','spatio_temporal_cluster_1samp_10000_0.001_OM-RS.png'))







                   
                   
                
                