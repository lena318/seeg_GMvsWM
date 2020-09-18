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

import sys
sys.path.append("..")
from filter_eeg_data import filter_eeg_data
import os



#assumes current working directory is paper005/paper005/pipelines/scripts/
raw_eeg_file_path = "../../../data_raw/eeg"
filtered_eeg_file_path = "../../../data_processed/eeg/montage/referential/filtered"

RID = [f for f in sorted(os.listdir(raw_eeg_file_path))]

for sub_ID in RID:
    sub_ID_eeg_file_path = os.path.join(raw_eeg_file_path, sub_ID)
    raw_eeg_files = [f for f in sorted(os.listdir(sub_ID_eeg_file_path))]
    inputfile_path = os.path.join(raw_eeg_file_path, sub_ID)
    outputfile_path = os.path.join(filtered_eeg_file_path, sub_ID)
    for raw_eeg_file in raw_eeg_files:
        outputfile_name = raw_eeg_file.replace('.pickle', '_filtered.pickle')
        if not(os.path.isdir(outputfile_path)): os.mkdir(outputfile_path)
        inputfile = os.path.join(inputfile_path, raw_eeg_file)
        outputfile = os.path.join(outputfile_path,outputfile_name)
        filter_eeg_data(inputfile, outputfile)



