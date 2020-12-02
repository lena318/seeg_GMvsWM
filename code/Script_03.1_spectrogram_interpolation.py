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
np.seterr(divide = 'ignore')

paper_path = "../../.."
#assumes current working directory is paper005/paper005/pipelines/scripts/
filtered_eeg_file_path = os.path.join(paper_path, "data_processed/eeg/montage/referential/filtered")
spectrogram_file_path = os.path.join(paper_path, "data_processed/spectrogram/montage/referential/filtered/original")
spectrogram_file_path_interpolated = os.path.join(paper_path, "data_processed/spectrogram/montage/referential/filtered/interpolated")
RID = [f for f in sorted(os.listdir(filtered_eeg_file_path))]

for sub_ID in RID:
    sub_ID_eeg_file_path = os.path.join(filtered_eeg_file_path, sub_ID)
    eeg_files = [f for f in sorted(os.listdir(sub_ID_eeg_file_path))]
    inputfile_path = os.path.join(filtered_eeg_file_path, sub_ID)
    outputfile_path = os.path.join(spectrogram_file_path, sub_ID)
    outputfile_path_interpolated = os.path.join(spectrogram_file_path_interpolated, sub_ID)
    if not (os.path.isdir(outputfile_path)): os.mkdir(outputfile_path)
    if not (os.path.isdir(outputfile_path_interpolated)): os.mkdir(outputfile_path_interpolated)
    for eeg_file in eeg_files:
        outputfile_name = eeg_file.replace('EEG_filtered.pickle', 'spectrogram_filtered.pickle')
        outputfile_name_interpolation = eeg_file.replace('EEG_filtered.pickle', 'spectrogram_filtered_interpolated.pickle')
        inputfile = os.path.join(inputfile_path, eeg_file)
        outputfile = os.path.join(outputfile_path,outputfile_name)
        outputfile_interpolation = os.path.join(outputfile_path_interpolated, outputfile_name_interpolation)
        print("reading file {0}".format(inputfile))
        with open(inputfile, 'rb') as f: data, fs = pickle.load(f)
        data_array = np.array(data)
        NFFT = int(256)
        noverlap = 128
        spec_shape = plt.specgram(x=data_array[:,0], Fs=fs, NFFT=NFFT, scale_by_freq=True, noverlap=noverlap)[0].shape
        interpolation_length = 200
        spectrogram = np.zeros(shape = (spec_shape[0], spec_shape[1], data_array.shape[1]) )
        spectrogram_interpolation = np.zeros(shape=(spec_shape[0], interpolation_length, data_array.shape[1]))
        for i in range(np.shape(data_array)[1]):
            x = data_array[:,i]
            Pxx, freqs, bins, im = plt.specgram(x=x, Fs=fs, NFFT=NFFT, scale_by_freq=True, cmap='rainbow', noverlap=noverlap)
            spectrogram[:,:,i] = Pxx
            #interpolation
            interp = interp1d(bins, Pxx, kind='nearest')
            xnew = np.linspace(bins[0], bins[-1], num=interpolation_length, endpoint=True)
            Pxx_interplation = interp(xnew)
            spectrogram_interpolation[:, :, i] = Pxx_interplation
        print("saving file {0}".format(outputfile_name))
        with open(outputfile, 'wb') as f: pickle.dump([spectrogram, freqs, bins, data.columns[:]], f)
        print("saving file {0}\n\n".format(outputfile_name_interpolation))
        with open(outputfile_interpolation, 'wb') as f: pickle.dump([spectrogram_interpolation, freqs, xnew, data.columns[:]], f)
