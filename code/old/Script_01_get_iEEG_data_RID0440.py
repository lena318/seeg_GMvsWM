"""
2020.06.10
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose: script to get iEEG data in batches

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Logic of code:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input:
  username: first argument. Your iEEG.org username
  password: second argument. Your iEEG.org password


Please edit all other inputs below based on needs
  iEEG_filename: The file name on iEEG.org you want to download from
  start_time_usec: the start time in the iEEG_filename. In microseconds
  stop_time_usec: the stop time in the iEEG_filename. In microseconds.
    iEEG.org needs a duration input: this is calculated by stop_time_usec - start_time_usec
  ignore_electrodes: the electrode/channel names you want to exclude. EXACT MATCH on iEEG.org. Caution: some may be LA08 or LA8
  outputfile: the path and filename you want to save.
    PLEASE INCLUDE EXTENSION .pickle.
    Use this output naming convention: sub-RIDXXXX_iEEGFILENAME_STARTTIME_STOPTIME_EEG.pickle
    example: 'sub-RID0278_HUP138_phaseII_415723190000_248525740000_EEG.pickle'

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output:
Saves EEG timeseries in specified output directors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:

python3.6 Script_01_get_iEEG_data_RID0440.py "username" "password"

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


import sys
sys.path.append("..")
from get_iEEG_data import get_iEEG_data
import os
import pickle


#assumes current working directory is paper005/paper005/pipelines/scripts
paper_directory_to_save_eeg = "../../../data_raw/eeg"
username=sys.argv[1]
password=sys.argv[2]

sub_ID="RID0440"
iEEG_filename="HUP172_phaseII"
ignore_electrodes=["EKG1","EKG2","LB12","RB08","RE01","RE08","RE09","RE10"]

start_times_array=[
356850680000,
402651841658,
402704260829,
402756680000]
stop_times_array=[
356903099171,
402704260829,
402756680000,
402936680000]

for i in range(len(start_times_array)):
  print("\n\n\nID: {0}\nStart: {1}\nStop: {2}".format(sub_ID,start_times_array[i],stop_times_array[i]))
  start_time_usec=start_times_array[i]
  stop_time_usec=stop_times_array[i]

  EEG_outputfile="{0}/sub-{1}/sub-{1}_{2}_{3}_{4}_EEG.pickle".format(paper_directory_to_save_eeg,sub_ID,iEEG_filename,start_time_usec,stop_time_usec)
  #print(os.path.isdir("{0}/sub-{1}".format(paper_directory_to_save_eeg,sub_ID)))
  #Get iEEG data
  get_iEEG_data(username,password,iEEG_filename, start_time_usec, stop_time_usec, ignore_electrodes, EEG_outputfile)
