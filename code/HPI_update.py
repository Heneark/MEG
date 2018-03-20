# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 17:50:20 2017

@author: oussama.abdoun
"""
#==============================================================================
# UPDATE FIF FILE WITH NEW HPI COORDINATES FROM A .hc file
#==============================================================================

#==============================================================================
import sys
sys.path.append("/dycog/meditation/ERC/Analyses/MEG/code/")
from header import *
#==============================================================================
warnings.filterwarnings("ignore",category=DeprecationWarning)

def update(task, subject, state, block, precision='0.5cm', opt='start'):
    """Load *-epo.fif Epochs from Analysis/<task>/meg/Epochs/<subject>/.
    Coordinates are updated according to .hc files in Analysis/<task>/meg/Coregistration/<subject>/.
    Corresponding bad segments due to head movement are rejected.
    The average *-ave.fif Evoked response is then saved in Analysis/<task>/meg/Evoked/<subject>/.
    Parameters:
        opt: argument to choose between multiple .hc files
    """
    
    #==============================================================================
    data_path = op.join(Analysis_path, task, 'meg', 'Epochs', subject)
    hc_path = op.join(Analysis_path, task, 'meg', 'Coregistration', subject)
    #==============================================================================
    
    # Find the .hc file corresponding to this block
    hc_files = glob.glob(op.join(hc_path, '*{}*_{}*.hc'.format(precision,block)))
    hc_file = ''
    if len(hc_files) > 1:
        for hc in hc_files:
            if not opt and not '+' in hc:
                hc_file = hc
                continue
            if opt in hc:
                hc_file = hc
                continue
        if not hc_file:
            hc_file = hc_files[0]
            warnings.warn("Multiple .hc files found, none matching the desired option. Using {}.".format(hc_file))
    else:
        hc_file = hc_files[0]
    
    bs_fname = hc_file[:-3] + '-bad.segments' + block
    bads = False
    if os.stat(bs_fname).st_size > 0:
        bad_segments = pd.read_table(bs_fname, header=None)
        bads = True
    
    # Read and load the content of the .hc file    
#==============================================================================
    with open(hc_file) as f:
        content = f.readlines()

    content = [x.strip() for x in content] 
    
    Nasion_ctf = [0,0,0,1]
    Left_ctf = [0,0,0,1]
    Right_ctf = [0,0,0,1]
    
    Nasion_ctf[0] = float(content[25].strip('x = '))
    Nasion_ctf[1] = float(content[26].strip('y = '))
    Nasion_ctf[2] = float(content[27].strip('z = '))
    
    Left_ctf[0] = float(content[29].strip('x = '))
    Left_ctf[1] = float(content[30].strip('y = '))
    Left_ctf[2] = float(content[31].strip('z = '))
    
    Right_ctf[0] = float(content[33].strip('x = '))
    Right_ctf[1] = float(content[34].strip('y = '))
    Right_ctf[2] = float(content[35].strip('z = '))
#==============================================================================
    
    files = glob.glob(op.join(data_path, '*{}-epo.fif'.format(block)))
        
        
    for data_fname in files:
        
        print(('Updating [' + op.split(data_fname)[1] + '] with [' + op.split(hc_file)[1] + ']'))        
        
        data = mne.read_epochs(data_fname)
        
        info = data.info
        
        # CTF -> Neuromag Coordinate transformation matrix
        #==============================================================================
        ctf_head_trans = np.asmatrix(list(info['ctf_head_t'].items())[1][1])
        
        # New Neuromag coordinates
        
        Nasion_head = np.asarray(np.dot(ctf_head_trans, Nasion_ctf)*0.01)[0][:3]
        
        Left_head = np.asarray(np.dot(ctf_head_trans, Left_ctf)*0.01)[0][:3]
        
        Right_head = np.asarray(np.dot(ctf_head_trans, Right_ctf)*0.01)[0][:3]
        #==============================================================================
        
        
        
        # DEV <-> CTF Coordinate transformation matrix
        #==============================================================================
        dev_ctf_trans = np.asmatrix(list(info['dev_ctf_t'].items())[1][1])
        ctf_dev_trans = inv(dev_ctf_trans)
        
        # New Dev coordinates                
        
        Nasion_dev = np.asarray(np.dot(ctf_dev_trans, Nasion_ctf)*0.01)[0][:3]
        
        Left_dev = np.asarray(np.dot(ctf_dev_trans, Left_ctf)*0.01)[0][:3]
        
        Right_dev = np.asarray(np.dot(ctf_dev_trans, Right_ctf)*0.01)[0][:3]
        #==============================================================================
        
        
        
        # UPDATE THE INFO 'hpi_results' DATA WITH NEW DEV COORDINATES
        #==============================================================================
        for i in range(len(data.info['hpi_results'][0]['dig_points'])):
            
            #Just checking it's the right kind of digitizer ...
            if data.info['hpi_results'][0]['dig_points'][i]['kind'] == 1:
                
                #Left, Nasion, Right
                if data.info['hpi_results'][0]['dig_points'][i]['ident'] == 1:
                    
                    data.info['hpi_results'][0]['dig_points'][i]['r'] = Left_dev
                    #data.info['hpi_results'][0]['dig_points'][i]['coord_frame'] = 1
                    
                if data.info['hpi_results'][0]['dig_points'][i]['ident'] == 2:
                    
                    data.info['hpi_results'][0]['dig_points'][i]['r'] = Nasion_dev
                    #data.info['hpi_results'][0]['dig_points'][i]['coord_frame'] = 1
                    
                if data.info['hpi_results'][0]['dig_points'][i]['ident'] == 3:
                    
                    data.info['hpi_results'][0]['dig_points'][i]['r'] = Right_dev
                    #data.info['hpi_results'][0]['dig_points'][i]['coord_frame'] = 1
        
        #==============================================================================            
                    
                    
        
        # UPDATE THE INFO 'dig' DATA WITH NEUROMAG HEAD COORDINATES
        #==============================================================================
        for i in range(len(data.info['dig'])):
            
            #Just checking it's the right kind of digitizer ...
            if data.info['dig'][i]['kind'] == 1:
                
                #Left, Nasion, Right
                if data.info['dig'][i]['ident'] == 1:
                    
                    data.info['dig'][i]['r'] = Left_head
                    data.info['dig'][i]['coord_frame'] = 4
                    
                if data.info['dig'][i]['ident'] == 2:
                    
                    data.info['dig'][i]['r'] = Nasion_head
                    data.info['dig'][i]['coord_frame'] = 4
                    
                if data.info['dig'][i]['ident'] == 3:
                    
                    data.info['dig'][i]['r'] = Right_head
                    data.info['dig'][i]['coord_frame'] = 4
                    
        #==============================================================================   
        if bads:
            for i in range(len(bad_segments)):
                events = data.events[:,0]/data.info['sfreq'] #extract event timing
                data.drop(np.logical_and(bad_segments.iloc[i][0] < events + data.times[-1], events + data.times[0] < bad_segments.iloc[i][1]), reason = 'head_movement')
                #reject the epoch if it ends after the beginning of the bad segment, and starts before the end of the bad segment
        
        #==============================================================================
        # AVERAGE AND SAVE TO FILE   
        evoked_path = op.join(Analysis_path, task, 'meg', 'Evoked', subject)
        if not op.isdir(evoked_path):
            os.makedirs(evoked_path)
        
        evoked = data.average()
        evoked.save(op.join(evoked_path, data_fname.split('/')[-1][:-8] + '-ave.fif'))
        
        print(colored ('Done', 'green'))
        #==============================================================================



