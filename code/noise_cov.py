# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:43:30 2017

@author: oussama.abdoun
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:43:30 2017

@author: oussama.abdoun
"""


#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# # # # # # # # # # # # # # #  OUT OF DATE  # # # # # # # # # # # # # # ## # # 
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================


import os
import os.path as op
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/SMEG/code/")
from csv_io_alex import *
from io_alex import *
import mne
from termcolor import colored


mne.set_log_level(verbose = 'WARNING')


task = 'SMEG'

subjects = get_subjlist(task=task)

for subject in subjects:

    data_path = '/dycog/meditation/ERC/Raw data/' + get_task(subject) + '/meg/'
    pathBase = '/dycog/meditation/ERC/Analyses/' + get_task(subject) + '/data/' + get_id(subject) 
    cov_path = pathBase + '/Covariance/'
    
    d = get_date(subject)
    subjlist = get_subjlist(date=d)
    subjlist.remove(subject)
    
    
    if not op.exists(cov_path):
        os.makedirs(cov_path)
        
         
    try:
        raw_empty_room = mne.io.read_raw_ctf(op.join(data_path, get_rawpath(subject, noise=1)))
    except:
        print ('Subject {} :'.format(subject))
        print colored('Searching for Noise file...', 'green')
        
        try:
            raw_empty_room = mne.io.read_raw_ctf(op.join(data_path, get_rawpath(subject, noise=2)))
    
        except:
            print colored('No Noise file found in this subject\'s directory...', 'red')
            d = get_date(subject)
            print colored('Trying to find Noise file in other subjects directories for date = {}/{}/{}...'.format(d[6:8], d[4:6], d[0:4]), 'green')
            
            if subjlist:
                print colored('Found {} other subject(s)'.format(len(subjlist)), 'green')
                
                for subj in subjlist:
                    print colored('Looking for Noise file in {}\'s directory'.format(get_name(subj)), 'green')
                    
                    try:
                        raw_empty_room = mne.io.read_raw_ctf(op.join(data_path, get_rawpath(subj, noise=1)))
                        
                    except:
                        
                        try:
                            raw_empty_room = mne.io.read_raw_ctf(op.join(data_path, get_rawpath(subj, noise=2)))
                        except:
                            print colored('No Noise file for {}.'.format(get_name(subj)), 'red')
                        else:
                            blocks = []
                    
                            for key in get_blocks(subject):
                                blocks.append(key)    
                    
                            for blk in blocks:
                                raw_fname = op.join(data_path, get_rawpath(subject) + blk + '.ds')
                                raw = mne.io.read_raw_ctf(raw_fname)
                        
                                noise_cov = mne.compute_raw_covariance(raw_empty_room, tmin=0, tmax=None)
                        
                                fname = os.path.join(pathBase + '/Covariance/', get_id(subject) + '_covariance' + '_blk' + blk + '-cov.fif')
                                mne.write_cov(fname, noise_cov)
                        
                            print ('DONE')
                            
                            break
                        
                    else:
                        blocks = []
                    
                        for key in get_blocks(subject):
                            blocks.append(key)    
                            
                        for blk in blocks:
                            raw_fname = op.join(data_path, get_rawpath(subject) + blk + '.ds')
                            raw = mne.io.read_raw_ctf(raw_fname)
                        
                            noise_cov = mne.compute_raw_covariance(raw_empty_room, tmin=0, tmax=None)
                        
                            fname = os.path.join(pathBase + '/Covariance/', get_id(subject) + '_covariance' + '_blk' + blk + '-cov.fif')
                            mne.write_cov(fname, noise_cov)
                            
                        print ('DONE')
                            
                        break
                
                        
            else:
                print colored('No other subjects found for this date.', 'red')        
            
        else:
            blocks = []
        
            for key in get_blocks(subject):
                blocks.append(key)    
        
            for blk in blocks:
                raw_fname = op.join(data_path, get_rawpath(subject) + blk + '.ds')
                raw = mne.io.read_raw_ctf(raw_fname)
                
                noise_cov = mne.compute_raw_covariance(raw_empty_room, tmin=0, tmax=None)
                
                fname = os.path.join(pathBase + '/Covariance/', get_id(subject) + '_covariance' + '_blk' + blk + '-cov.fif')
                mne.write_cov(fname, noise_cov)
                
            print ('DONE')
        
    else:
        blocks = []
        
        for key in get_blocks(subject):
            blocks.append(key)    
        
        for blk in blocks:
            raw_fname = op.join(data_path, get_rawpath(subject) + blk + '.ds')
            raw = mne.io.read_raw_ctf(raw_fname)
            
            noise_cov = mne.compute_raw_covariance(raw_empty_room, tmin=0, tmax=None)
            
            fname = os.path.join(pathBase + '/Covariance/', get_id(subject) + '_covariance' + '_blk' + blk + '-cov.fif')
            mne.write_cov(fname, noise_cov)
        
        print ('Subject {} :'.format(subject) + ' DONE')

            