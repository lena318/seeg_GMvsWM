"""
2020.08.01
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose:
    Compute spectrogram and interpolate each segment to a length of 200

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
from scipy import signal
np.seterr(divide = 'ignore')

paper_path = "../../.."
#assumes current working directory is paper005/paper005/pipelines/scripts/
filtered_eeg_file_path = os.path.join(paper_path, "data_processed/eeg/montage/referential/filtered")
data_output_path = os.path.join(paper_path, "data_processed/PSDvsDistance/montage/referential/filtered/original")
RID = [f for f in sorted(os.listdir(filtered_eeg_file_path))]
electrode_localization_file_path = os.path.join(paper_path, "data_raw/electrode_localization")

for sub_ID in RID:
    sub_ID_eeg_file_path = os.path.join(filtered_eeg_file_path, sub_ID)
    eeg_files = [f for f in sorted(os.listdir(sub_ID_eeg_file_path))]
    inputfile_path = os.path.join(filtered_eeg_file_path, sub_ID)
    sub_ID_electrode_localization_file_path = os.path.join(electrode_localization_file_path, sub_ID)
    electrode_localization_file_name = [f for f in sorted(os.listdir(sub_ID_electrode_localization_file_path))]
    electrode_localization = pd.read_csv(os.path.join(electrode_localization_file_path, sub_ID, electrode_localization_file_name[0]))
    outputfile_path = os.path.join(data_output_path, sub_ID)
    if not (os.path.isdir(outputfile_path)): os.mkdir(outputfile_path)
    for eeg_file in eeg_files:
        outputfile_name = eeg_file.replace('EEG_filtered.pickle', 'PSDvsDistance_filtered.pickle')
        inputfile = os.path.join(inputfile_path, eeg_file)
        outputfile = os.path.join(outputfile_path,outputfile_name)
        print("reading file {0}".format(inputfile))
        with open(inputfile, 'rb') as f: data, fs = pickle.load(f)
        data_array = np.array(data)
        columns = data.columns
        win = 4 * fs
        count = 0
        freq_shape = np.shape(signal.welch(data_array[:,1], fs, nperseg=win))[1]
        PSD = np.zeros(shape=(freq_shape,1))
        Distance = np.zeros(shape=(1,1))
        for e in columns:
            if any(np.array(electrode_localization.iloc[:, 0]) == e):#if any analyzed signals are not in the electrode localization file, dont' compute.
                loc = np.where(np.array(electrode_localization.iloc[:, 0]) == e)[0][0]
                dist = np.reshape(electrode_localization.iloc[loc, 4], newshape=(1,1))
                freqs, Pxx = signal.welch(data_array[:,count], fs, nperseg=win)
                Pxx = np.reshape(Pxx, newshape=(freq_shape, 1))
                PSD = np.append(PSD, Pxx, axis=1)
                Distance = np.append(Distance, dist, axis=1)
            count = count + 1
        PSD = np.delete(PSD, [0], axis=1)
        Distance = np.delete(Distance, [0], axis=1)
        Distance = Distance.flatten()
        #interpolation
        max = np.max(Distance)
        interp = interp1d(Distance, PSD, kind='nearest', axis=1)
        Distance_interpolated = np.arange(start=0, stop=max, step = 0.01)
        PSD_interplation = interp(Distance_interpolated)
        I = pd.Index(freqs, name="frequency")
        C = pd.Index(Distance, name="distance")
        PSD_df = pd.DataFrame(PSD, index=I, columns=C)
        I = pd.Index(freqs, name="frequency")
        C = pd.Index(np.round(Distance_interpolated, 2), name="distance")
        PSD_interplation_df = pd.DataFrame(PSD_interplation, index=I, columns=C)
        #plot = np.delete(np.log10(np.array(PSD_interplation_df)), range(50, np.shape(PSD_interplation)[0]), axis=0)
        #plot = np.delete(plot, 88, axis=1)
        #sns.heatmap(plot, cmap=sns.color_palette("Spectral_r", n_colors=20, desat=0.6))
        #plt.show()
        print("saving file {0}\n\n".format(outputfile))
        with open(outputfile, 'wb') as f: pickle.dump([PSD_df, PSD_interplation_df], f)


