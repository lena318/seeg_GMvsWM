"""
2020.08.01
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose:
    Averaging PSD vs Distance arrays of all patients

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Logic of code:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:

~~~~~~~
"""

import pickle
import numpy as np
import os
import sys
sys.path.append("..")
import os
import matplotlib.pyplot as plt
import pandas as pd

paper_path = "../../.."
#assumes current working directory is paper005/paper005/pipelines/scripts/
input_path = os.path.join(paper_path, "data_processed/PSDvsDistance/montage/referential/filtered/original")
output_path = os.path.join(paper_path, "data_processed/PSDvsDistance/montage/referential/filtered/averaged")
RID = [f for f in sorted(os.listdir(input_path))]

PSDvsDistance_avg = np.zeros(shape = (201,813, 4))
PSDvsDistance = np.empty(shape = (201,906, len(RID), 4))
PSDvsDistance[:] = np.NaN

cnt = 0
for sub_ID in RID:
    sub_ID_input_path = os.path.join(input_path, sub_ID)
    input_files = [f for f in sorted(os.listdir(sub_ID_input_path))]
    count = 0#assumes 0 = interictal, 1= preictal, 2= ictal, 3 = postictal
    for input_file in input_files:
        inputfile = os.path.join(sub_ID_input_path, input_file)
        print("reading file {0}".format(input_file))
        with open(inputfile, 'rb') as f: PSD_df, PSD_interplation_df = pickle.load(f)
        index = np.where(PSD_interplation_df.index == 50)[0][0]
        len_index = len(PSD_interplation_df.index)
        len_col = len(PSD_interplation_df.columns)
        print("{0}".format(len_index))
        print("{0}\n\n".format(len_col))
        if len_index > 50:
            end = len_index
            PSD_interplation_df_new = PSD_interplation_df.drop(PSD_interplation_df.index[range(index+1, end)])
        PSDvsDistance[:,np.array(range(0,len_col )),cnt, count] = PSD_interplation_df_new
        count = count + 1
    cnt = cnt + 1

output_file_name = "PSDvsDistance_All_patients.pickle"
outputfile = os.path.join(output_path,output_file_name )

interictal = PSDvsDistance[:,:,:,0]
preictal = PSDvsDistance[:,:,:,1]
ictal = PSDvsDistance[:,:,:,2]
postictal = PSDvsDistance[:,:,:,3]
print("saving file {0}\n\n".format(outputfile))
with open(outputfile, 'wb') as f: pickle.dump([interictal, preictal, ictal ,postictal], f)


