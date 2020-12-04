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
import combine_FIRST_and_FAST_images as combo
import segmentation_FIRST_and_FAST_images as seg

#%% Input/Output Paths and File names
ifname_EEG_times = ospj(path,"data/data_raw/iEEG_times/EEG_times.xlsx")
ifpath_electrode_localization = ospj( path, "data/data_raw/electrode_localization")
ifpath_imaging = ospj( path, "data/data_raw/imaging")
ofpath_segmentation = ospj( path, "data/data_processed/GM_WM_segmentations")

if not (os.path.isdir(ofpath_segmentation)): os.mkdir(ofpath_segmentation)

#%% Load Study Meta Data
data = pd.read_excel(ifname_EEG_times)    

#%% Processing Meta Data: extracting sub-IDs

sub_IDs_unique = np.unique(data.RID)


#%% Get Electrode Localization. Find which region each electrode is in for each atlas
for i in range(len(sub_IDs_unique)):
    #parsing data DataFrame to get iEEG information
    sub_ID = sub_IDs_unique[i]
    print("\nSubject: {0}".format(sub_ID))
    #getting preop3T and T00 images
    
    ifname_preop3T = ospj(ifpath_imaging, "sub-{0}".format(sub_ID), "ses-preop3T", "anat", "sub-{0}_ses-preop3T_acq-3D_T1w".format(sub_ID))
    ifname_T00 = ospj(ifpath_electrode_localization, "sub-{0}".format(sub_ID), "sub-{0}_T00_mprage".format(sub_ID))
    ofpath_segmentation_sub_ID = ospj(ofpath_segmentation, "sub-{0}".format(sub_ID))
    
    if not (os.path.isdir(ofpath_segmentation_sub_ID)): os.mkdir(ofpath_segmentation_sub_ID)
    
    ofbase_flirt = ospj(ofpath_segmentation_sub_ID, "sub-{0}_preop3T_to_T00_std_linear".format(sub_ID))
    ofbase_fnirt = ospj(ofpath_segmentation_sub_ID, "sub-{0}_preop3T_to_T00_std_nonlinear".format(sub_ID))

    ifname_first =  ospj(ofpath_segmentation_sub_ID, "sub-{0}_ses-preop3T_acq-3D_T1w_std_bet_subcort_all_fast_firstseg.nii.gz".format(sub_ID))
    ifname_fast =  ospj(ofpath_segmentation_sub_ID,"sub-{0}_ses-preop3T_acq-3D_T1w_std_bet_seg.nii.gz".format(sub_ID))

    ofname_FIRST_FAST_COMBINED = ospj(ofpath_segmentation_sub_ID,"sub-{0}_ses-preop3T_acq-3D_T1w_std_bet_GM_WM_CSF.nii.gz".format(sub_ID))
    ofname_FIRST_FAST_COMBINED_to_T00 = ospj(ofpath_segmentation_sub_ID,"sub-{0}_preop3T_to_T00_std_GM_WM_CSF.nii.gz".format(sub_ID))
    
    #First and Fast segmentation (first is subcortical structures, fast is cortical gray matter)
    if not (os.path.exists(ifname_fast)):
        seg.first_and_fast_segmentation(ifname_preop3T, ifname_T00, ofpath_segmentation_sub_ID)
        
    #Combine First and Fast images
    if not (os.path.exists(ofname_FIRST_FAST_COMBINED)):
        combo.combine_first_and_fast(ifname_first, ifname_fast, ofname_FIRST_FAST_COMBINED)


    if not (os.path.exists( "{0}.nii.gz".format(ofbase_fnirt)  )):
        seg.register_preop3T_to_T00(ifname_preop3T, ifname_T00, ofpath_segmentation_sub_ID, ofbase_flirt, ofbase_fnirt)
        
    if not (os.path.exists(ofname_FIRST_FAST_COMBINED_to_T00)):
        seg.applywarp_to_combined_first_fast(ofname_FIRST_FAST_COMBINED, ifname_T00, ofbase_fnirt, ofpath_segmentation_sub_ID, ofname_FIRST_FAST_COMBINED_to_T00)
        



#%%












