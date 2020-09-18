"""
2020.08.01
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose:
    Compute SNR

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Logic of code:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:

interictal_file = "./data_processed/eeg/montage/referential/filtered/sub-RID0278/sub-RID0278_HUP138_phaseII_394423190000_394512890000_EEG_filtered.pickle"
preictal_file = "./data_processed/eeg/montage/referential/filtered/sub-RID0278/sub-RID0278_HUP138_phaseII_415933490000_416023190000_EEG_filtered.pickle"
ictal_file = "./data_processed/eeg/montage/referential/filtered/sub-RID0278/sub-RID0278_HUP138_phaseII_416023190000_416112890000_EEG_filtered.pickle"
postictal_file = "./data_processed/eeg/montage/referential/filtered/sub-RID0278/sub-RID0278_HUP138_phaseII_416112890000_416292890000_EEG_filtered.pickle"

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
from scipy.interpolate import interp1d
import seaborn as sns
from scipy import signal
np.seterr(divide = 'ignore')

paper_path = "../../.."
#assumes current working directory is paper005/paper005/pipelines/scripts/
inputfile_path_SNR = os.path.join(paper_path, "data_processed/SNR/montage/referential/filtered/original")
outputfile_path_SNR_avg = os.path.join(paper_path, "data_processed/SNR/montage/referential/filtered/average")
RID = [f for f in sorted(os.listdir(inputfile_path_SNR))]

count_patient = 0
breaks = [0,3,6]
interpolation_length = 200
SNR_all_patients = np.zeros(shape=(interpolation_length, len(breaks)+1,4, len(RID)))#4 = interictal, preictal, ictal, postictal
for sub_ID in RID:
    inputfile_path_SNR_sub_ID = os.path.join(inputfile_path_SNR, sub_ID)
    inputfile_name_SNR_files = [f for f in sorted(os.listdir(inputfile_path_SNR_sub_ID))]
    #outputfile_path_SNR_avg_sub_ID = os.path.join(outputfile_path_SNR_avg, sub_ID)
    #if not (os.path.isdir(outputfile_path_SNR_avg_sub_ID)): os.mkdir(outputfile_path_SNR_avg_sub_ID)
    count_ictal = 0
    for SNR_file in inputfile_name_SNR_files:
        outputfile_name = SNR_file.replace('SNR.pickle', 'SNR_average.pickle')
        inputfile = os.path.join(inputfile_path_SNR_sub_ID, SNR_file)
        #outputfile = os.path.join(outputfile_path_SNR_avg_sub_ID,outputfile_name)
        print("reading file {0}".format(inputfile))
        with open(inputfile, 'rb') as f: SNR, SNR_interpolation, SNR_avg_period, bins,bins_new, columns, dist_vec = pickle.load(f)
        dist_breaks = []
        for d in range(len(breaks)):
            if d == 0:
                dist_breaks.append(np.where(dist_vec[:,0] == breaks[d]))
            if d > 0:
                dist_breaks.append(np.where(np.logical_and(dist_vec[:,0]>breaks[d-1], dist_vec[:,0]<=breaks[d])))
        dist_breaks.append(np.where(dist_vec[:,0] > breaks[len(breaks)-1]))
        SNR_mean = np.zeros(shape=(len(bins), len(breaks)+1))
        for d in range(len(breaks)+1):
            if len(SNR[dist_breaks[d], :][0]) > 0:
                SNR_mean[:,d]= np.nanmean(SNR[dist_breaks[d], :][0], axis=0)
            if len(SNR[dist_breaks[d], :][0]) == 0:
                SNR_mean[:, d] = np.nan
        interp = interp1d(bins, SNR_mean, kind='linear', axis=0)
        bins_new = np.linspace(bins[0], bins[-1], num=interpolation_length, endpoint=True)
        SNR_mean_interplation = interp(bins_new)
        SNR_all_patients[:,:,count_ictal, count_patient] = SNR_mean_interplation
        count_ictal = count_ictal + 1
    count_patient = count_patient + 1


SNR_all_patients_mean = np.nanmean(SNR_all_patients, axis=3)
outputfile = os.path.join(outputfile_path_SNR_avg, "SNR_all_patients_mean.pickle")
print("saving file {0}".format(outputfile))
with open(outputfile, 'wb') as f: pickle.dump([SNR_all_patients_mean, SNR_all_patients, bins_new, RID, breaks], f)

