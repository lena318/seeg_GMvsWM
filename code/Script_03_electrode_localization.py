"""
2020.06.10
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose: 
    1. This is a wrapper script: Runs through meta-data to automatically calculate for all data
        Meta-data: data_raw/iEEG_times/EEG_times.xlsx
    2. Get electrode localization: Find which region each electrode is in for each atlas
    3. Calls electrode_localization.by_atlas in paper001/code/tools
        See this function for more detail

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Logic of code:
    1. Import appropriate tools 
    2. Get appropriate input and output paths and file names
    3. Setting appropriate parameters and preprocessing of data
    4. Get electrode localization
        1. Get electrode localization file for each patient
        2. For standard atlases:
            1. Get atlas image path
            2. Call electrode_localization.by_atlas
        3. For random atlases:
            1. For each random atlas permutation
                1. Get atlas image path
                2. Call electrode_localization.by_atlas
                

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input:
    None. This is a wrapper scipt that automatically runs based on meta-data file
    Meta-data: data_raw/iEEG_times/EEG_times.xlsx
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output:
    Saves electrode localization by atlas for each patinet's electrode localization file 
    in appropriate directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:

    python3.6 Script_03_electrode_localization.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

#%%
path = "/media/arevell/sharedSSD1/linux/papers/paper005" #Parent directory of project
import sys
import os
import pandas as pd
import numpy as np
from os.path import join as ospj
sys.path.append(ospj(path, "seeg_GMvsWM", "code", "tools"))
import electrode_localization

#%% Input/Output Paths and File names
ifname_EEG_times = ospj( path, "data/data_raw/iEEG_times/EEG_times.xlsx")
ifpath_electrode_localization = ospj( path, "data/data_raw/electrode_localization")
ifpath_segmentations = ospj( path, "data/data_processed/GM_WM_segmentations")
ofpath_electrode_localization = ospj( path, "data/data_processed/electrode_localization")

if not (os.path.isdir(ofpath_electrode_localization)): os.makedirs(ofpath_electrode_localization, exist_ok=True)
#%% Load Study Meta Data
data = pd.read_excel(ifname_EEG_times)    

#%% Processing Meta Data: extracting sub-IDs

sub_IDs_unique = np.unique(data.RID)


#%% Get Electrode Localization. Find which region each electrode is in for each atlas
for i in range(len(sub_IDs_unique)):
    #parsing data DataFrame to get iEEG information
    sub_ID = sub_IDs_unique[i]
    print("\nSubject: {0}".format(sub_ID))
   
    
    #getting electrode localization file
    ifpath_electrode_localization_sub_ID = ospj(ifpath_electrode_localization, "sub-{0}".format(sub_ID))
    if not os.path.exists(ifpath_electrode_localization_sub_ID): print("Path does not exist: {0}".format(ifpath_electrode_localization_sub_ID))
    ifname_electrode_localization_sub_ID = ospj(ifpath_electrode_localization_sub_ID, "sub-{0}_electrodenames_coordinates_native_and_T1.csv".format(sub_ID))
    
    #getting atlas/segmentation files
    ifpath_seg_sub_ID = ospj(ifpath_segmentations, "sub-{0}".format(sub_ID))
    if not os.path.exists(ifpath_seg_sub_ID): print("Path does not exist: {0}".format(ifpath_seg_sub_ID))
    ifname_seg_sub_ID = ospj(ifpath_seg_sub_ID, "sub-{0}_preop3T_to_T00_std_GM_WM_CSF.nii.gz".format(sub_ID))
    

    #Output electrode localization
    ofpath_electrode_localization_sub_ID = ospj(ofpath_electrode_localization, "sub-{0}".format(sub_ID))
    if not (os.path.isdir(ofpath_electrode_localization_sub_ID)): os.mkdir(ofpath_electrode_localization_sub_ID)#if the path doesn't exists, then make the directory
    ofpath_localization_files = ospj(ofpath_electrode_localization_sub_ID, "individual_files")
    if not (os.path.isdir(ofpath_localization_files)): os.mkdir(ofpath_localization_files)#if the path doesn't exists, then make the directory

    ofname_electrode_localization_concatenated = ospj(ofpath_electrode_localization_sub_ID, "sub-{0}_electrode_localization.csv".format(sub_ID))

    #localization
    electrode_localization.by_region(ifname_electrode_localization_sub_ID, ifname_seg_sub_ID, ospj(ofpath_localization_files, "sub-{0}_GM_WM_CSF.csv".format(sub_ID)))
    electrode_localization.distance_from_label(ifname_electrode_localization_sub_ID, ifname_seg_sub_ID, 2, ospj(ofpath_localization_files, "sub-{0}_WM_distance.csv".format(sub_ID)))

    #Concatenate files into one
    files = [f for f in sorted(os.listdir(ofpath_localization_files))]
    for f in range(len(files)):
        if f == 0:
            data = pd.read_csv(ospj(ofpath_localization_files, files[f]), sep=",", header=0)
        else:
            data= pd.concat([data, pd.read_csv(ospj(ofpath_localization_files, files[f]), sep=",", header=0) .iloc[:,4] ]  , axis = 1)
    pd.DataFrame.to_csv(data, ofname_electrode_localization_concatenated, header=True, index=False)
            
    
#%%












