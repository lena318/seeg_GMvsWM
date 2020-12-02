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
sys.path.append(".")
import os
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import seaborn as sns
from scipy import signal
np.seterr(divide = 'ignore')

paper_path = "."
#assumes current working directory is paper005/paper005/pipelines/scripts/
inputfile_path_eeg = os.path.join(paper_path, "data_processed/eeg/montage/referential/filtered")
outputfile_path_SNR = os.path.join(paper_path, "data_processed/SNR/montage/referential/filtered/original")
RID = [f for f in sorted(os.listdir(inputfile_path_eeg))]
inputfile_path_electrode_localization = os.path.join(paper_path, "data_raw/electrode_localization")

for sub_ID in RID:
    inputfile_path_eeg_sub_ID = os.path.join(inputfile_path_eeg, sub_ID)
    inputfile_name_eeg_files = [f for f in sorted(os.listdir(inputfile_path_eeg_sub_ID))]
    inputfile_path_electrode_localization_sub_ID = os.path.join(inputfile_path_electrode_localization, sub_ID)
    outputfile_path_SNR_sub_ID = os.path.join(outputfile_path_SNR, sub_ID)
    if not (os.path.isdir(outputfile_path_SNR_sub_ID)): os.mkdir(outputfile_path_SNR_sub_ID)
    count_ictal = 0
    for eeg_file in inputfile_name_eeg_files:
        outputfile_name = eeg_file.replace('EEG_filtered.pickle', 'SNR.pickle')
        inputfile = os.path.join(inputfile_path_eeg_sub_ID, eeg_file)
        outputfile = os.path.join(outputfile_path_SNR_sub_ID,outputfile_name)
        inputfile_name_electrode_localization = [f for f in sorted(os.listdir(inputfile_path_electrode_localization_sub_ID))]
        electrode_localization = pd.read_csv(os.path.join(inputfile_path_electrode_localization, sub_ID, inputfile_name_electrode_localization[0]))
        print("reading file {0}".format(inputfile))
        with open(inputfile, 'rb') as f: data, fs = pickle.load(f)
        data_array = np.array(data)
        win = 4 * fs
        NFFT = int(fs*4)
        noverlap = NFFT/2
        columns = data.columns
        dist_vec = np.zeros(shape=(len(columns), 1))
        for e in range(len(columns)):
            if any(np.array(electrode_localization.iloc[:,0]) == columns[e]):  # if any analyzed signals are not in the electrode localization file, dont' compute.
                loc = np.where(np.array(electrode_localization.iloc[:, 0]) == columns[e])[0][0]
                dist_GM = electrode_localization.iloc[loc, 4]
                dist_vec[e,0] = dist_GM
                #Freqs, Pxx = signal.welch(data_array[:, loc], fs, nperseg=win)
                #power = np.trapz(Pxx[range(2, 201)], dx=0.25)
                #sns.scatterplot(x = Freqs[range(2, 201)], y= Pxx[range(2, 201)])
                Pxx, Freqs, bins, im = plt.specgram(x=data_array[:, e], Fs=fs, NFFT=NFFT, scale_by_freq=True, cmap='rainbow', noverlap=noverlap)
                Pxx_mean = np.mean(Pxx, axis=1)
                power = np.trapz(Pxx_mean[range(2, np.where(Freqs==50)[0][0]+1  )], dx=Freqs[1]-Freqs[0])
                #sns.scatterplot(x=Freqs2[range(2, 201)], y=Pxx3[range(2, 201)])
                interpolation_length = 200
                if count_ictal == 0:
                    if e ==0:
                        power_interictal = np.zeros(shape=(len(columns), 1))
                    power_interictal[e, 0] = power
                if e == 0:
                    #initialize
                    SNR_time = np.zeros(shape=(len(columns),len(bins)))
                    SNR = np.zeros(shape=(len(columns), len(bins)))
                    SNR_interpolation = np.zeros(shape=(len(columns), interpolation_length))
                    SNR_avg_period = np.zeros(shape=(len(columns),1))
                SNR_time[e,:] = np.trapz(Pxx[range(2, np.where(Freqs==50)[0][0]+1  ),:], dx=Freqs[1]-Freqs[0], axis=0)/power_interictal[e,0]#/len(bins)
                interp = interp1d(bins, SNR_time, kind='nearest')
                bins_new = np.linspace(bins[0], bins[-1], num=interpolation_length, endpoint=True)
                SNR_time_interplation = interp(bins_new)
                SNR_avg_period[e,0] = power/power_interictal[e,0]
                #sns.scatterplot(x = bins2 , y = power3)
        SNR[:,:] = SNR_time #need to restructure data since ictal time != postictal time.
        SNR_interpolation[:,:] = SNR_time_interplation
        count_ictal = count_ictal + 1
        print("saving file {0}\n\n".format(outputfile))
        with open(outputfile, 'wb') as f: pickle.dump([SNR, SNR_interpolation, SNR_avg_period, bins,bins_new, columns, dist_vec], f)


fig1 = sns.scatterplot(x = np.ndarray.flatten(dist_vec), y= np.ndarray.flatten(SNR_avg_period)  )
fig1.set(xlabel='Dist from GM (mm)', ylabel='SNR', title = "RID0309 Ictal SNR")

plt.show()