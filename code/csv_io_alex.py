# -*- coding: utf-8 -*-

import pandas as pd
import collections
import os.path as op
from termcolor import colored

#==============================================================================
#SUBJECT INFOS USING CSV FILE
#==============================================================================

#PATH BASES
#==============================================================================
code_path = os.getcwd()
Analysis_path = op.dirname(op.dirname(code_path))
pathBase = op.dirname(Analysis_path)

#pathBase = '/dycog/meditation/ERC'
#Analysis_path = op.join(pathBase, 'Analyses')
#code_path = op.join(Analysis_path, 'MEG', 'code')

Raw_data_path = op.join(pathBase, 'Raw data')
#==============================================================================


#GET SUBJECT LIST
#==============================================================================
def get_subjlist(task=None, date=None, pathBase = Analysis_path):
     
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)
     
     ld = ld[ld.include == '1']
     if task:
         ld = ld[ld.task == task]
 
     if date:
         ld = ld[ld.date == date]

     ld = ld.id
     try:
         
         subjlist = [ld.iloc[0]]
         
         for k in range(1, ld.shape[0]):
             z = ld.iloc[k]
             if z != ld.iloc[k-1]:
                     subjlist.append(z)
         
         return subjlist
         
     except:
         
         print(colored('No subjects found for this task ...', 'green'))
#==============================================================================


#GET SUBJECT ID
#==============================================================================
def get_id(subj, pathBase = Analysis_path):
     
     if subj.isdigit() == True:
         return subj
         
     if subj.isalpha() == True:
         ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)
             
         ld = ld[ld.name == subj]
         subj_id = ld['id'].iloc[0]
             
         return subj_id
#==============================================================================
    

#GET SUBJECT PROTOCOL
#==============================================================================
def get_protocol(subj, pathBase = Analysis_path):
         
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)
     
     ld = ld[(ld.id == subj) | (ld.name == subj)]    
     
     protocol = ld.iloc[0]['protocol']
         
     return protocol
#==============================================================================








#GET SUBJECT TASK
#==============================================================================
def get_task(subj, block = None, protocol = None, pathBase = Analysis_path):
         
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)
     
     ld = ld[(ld.id == subj) | (ld.name == subj)]
     if protocol:
         ld = ld[ld.protocol == protocol]
     if block:
         ld = ld[ld.block == block]
     
     task = ld.iloc[0]['task']
         
     return task
#==============================================================================
    

#GET SUBJECT DATE
#==============================================================================
def get_date(subj, pathBase = Analysis_path):
          
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)
     
     ld = ld[(ld.id == subj) | (ld.name == subj)]    
     
     date = ld.iloc[0]['date']
         
     return date 
#==============================================================================
 


#GET SUBJECT'S NAME
#==============================================================================
def get_name(subj, pathBase = Analysis_path):
          
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)
     
     ld = ld[(ld.id == subj) | (ld.name == subj)]    
     
     name = ld.iloc[0]['name']
         
     return name
#==============================================================================
     
     
     
     
#GET RAW DATA PATH
#==============================================================================
def get_rawpath(subj, task=None, noise=0, NoiseBase = op.join(Raw_data_path,'MEG_noise'), pathBase = Analysis_path):
        
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)

     if subj.isdigit() == True:   
         ldsub = ld.loc[ld['id'] == subj]
         subj_id = subj
         subj_name = ldsub.iloc[0]['name']
     
     if subj.isalpha() == True:
         ldsub = ld.loc[ld['name'] == subj]
         subj_name = subj
         subj_id = ldsub.iloc[0]['id']
         
     if task:
         ldsub = ldsub[ld.task == task]
         
     protocol = ldsub.iloc[0]['protocol']
     date = ldsub.iloc[0]['date']
     
     if noise==1:
         NoisePath = op.join(NoiseBase, subj_id, subj_name + '_' + protocol + '_' + date + '.misc', 'lyon_Noise_' + date)
         for k in range(10):
             if op.exists(NoisePath + '_0' + str(k) + '.ds'):
                 filename = NoisePath + '_0' + str(k) + '.ds'
 
     if noise==0:
         filename =  subj_name + '_' + protocol + '_' + date + '_'
         
     return subj_id, filename
#==============================================================================


#GET AN ORDERED DICT OF A SUBJECT'S BLOCKS AND CORRESPONDING STATES
#==============================================================================
def get_blocks(subj, state=None, task=None, protocol=None, pathBase = Analysis_path):
     
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)

     ld = ld[(ld.id == subj) | (ld.name == subj)]
     
     if state:
         ld = ld[(ld.state == state) | (ld.state == state + '1') | (ld.state == state + '2')]
         
     if task:
         ld = ld[ld.task == task]
         
     if protocol:
         ld = ld[ld.protocol == protocol]
         
     blockstate = dict()
     
     for k in range(ld.shape[0]):
         
         blockstate[ld['block'].iloc[k]] = ld['state'].iloc[k]
         
     odblockstate = collections.OrderedDict(sorted(blockstate.items()))
     
     return odblockstate
#==============================================================================
    

#COUNT THE NUMBER OF SUBJ BLOCKS FOR A GIVEN STATE
#==============================================================================
def state_count(state, pathBase = Analysis_path):
     
     ld = pd.read_table(op.join(pathBase, 'MEG', 'meta', 'listdata.tsv'), dtype=str)
 
     ld = ld[((ld.state == state) | (ld.state == state + '1') | (ld.state == state + '2') | (ld.state == state + '_postFA') | (ld.state == state + '_postOM')) & (ld.include == '1')]
     
     return ld.shape[0]
#==============================================================================
