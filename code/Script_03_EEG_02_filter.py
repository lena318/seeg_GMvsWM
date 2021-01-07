"""
2020.08.01
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose:
 filter raw eeg data in data_raw folder in a batch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Logic of code:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:

    #preictal
    inputfile =
    outputfile =

    #ictal
    inputfile =
    outputfile =

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
path = "/media/arevell/sharedSSD1/linux/papers/paper005" 
import sys
import os
from os.path import join as ospj
sys.path.append(ospj(path, "seeg_GMvsWM", "code", "tools"))
sys.path.append(ospj(path, "seeg_GMvsWM", "code" ,"tools", "ieegpy"))
from filter_eeg_data import filter_eeg_data
import pandas as pd


ifname_EEG_times = ospj(path,"data/data_raw/iEEG_times/EEG_times.xlsx")
ifpath_EEG = ospj(path,"data/data_raw/EEG")
ofpath_filtered_eeg = ospj(path,"data/data_processed/eeg/montage/referential/filtered")

if not (os.path.isdir(ofpath_filtered_eeg)): os.makedirs(ofpath_filtered_eeg, exist_ok=True)

#%% Load Study Meta Data
data = pd.read_excel(ifname_EEG_times)    

##%


for i in range(len(data)):
    #parsing data DataFrame to get iEEG information
    sub_ID = data.iloc[i].RID
    iEEG_filename = data.iloc[i].file
    ignore_electrodes = data.iloc[i].ignore_electrodes.split(",")
    start_time_usec = int(data.iloc[i].connectivity_start_time_seconds*1e6)
    stop_time_usec = int(data.iloc[i].connectivity_end_time_seconds*1e6)
    descriptor = data.iloc[i].descriptor
    #input filename EEG
    ifpath_EEG_sub_ID = ospj(ifpath_EEG, "sub-{0}".format(sub_ID))
    #Output filtered EEG
    ofpath_EEG_sub_ID = ospj(ofpath_filtered_eeg, "sub-{0}".format(sub_ID))
    if not (os.path.isdir(ofpath_EEG_sub_ID)): os.mkdir(ofpath_EEG_sub_ID)#if the path doesn't exists, then make the directory
    ifname_EEG = "{0}/sub-{1}_{2}_{3}_{4}_EEG.pickle".format(ifpath_EEG_sub_ID, sub_ID, iEEG_filename, start_time_usec, stop_time_usec)
    ofname_EEG_filtered = "{0}/sub-{1}_{2}_{3}_{4}_EEG_filtered.pickle.pickle".format(ofpath_EEG_sub_ID, sub_ID, iEEG_filename, start_time_usec, stop_time_usec)
    print("\n\n\nID: {0}\nDescriptor: {1}".format(sub_ID, descriptor))
    if (os.path.exists(ofname_EEG_filtered)):
        print("File already exists: {0}".format(ofname_EEG_filtered))
    else:#if file already exists, don't run below
        filter_eeg_data(ifname_EEG, ofname_EEG_filtered)


    

