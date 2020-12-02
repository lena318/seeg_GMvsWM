"""
2020.08.01
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose:
    Averaging interpolated spectrogram of all patients

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
spectrogram_file_path_interpolated = os.path.join(paper_path, "data_processed/spectrogram/montage/referential/filtered/interpolated")
spectrogram_file_path_interpolated_avg = os.path.join(paper_path, "data_processed/spectrogram/montage/referential/filtered/interpolated_avg")
RID = [f for f in sorted(os.listdir(spectrogram_file_path_interpolated))]
electrode_localization_file_path = os.path.join(paper_path, "data_raw/electrode_localization")
total_num_electrodes_GM = np.zeros(shape=[1, len(RID)])
total_num_electrodes_WM = np.zeros(shape=[1, len(RID)])
cnt = 0
distance_considered_WM = 3
for sub_ID in RID:
    sub_ID_spectrogram_file_path = os.path.join(spectrogram_file_path_interpolated, sub_ID)
    spectrogram_files = [f for f in sorted(os.listdir(sub_ID_spectrogram_file_path))]
    inputfile_path = os.path.join(spectrogram_file_path_interpolated, sub_ID)
    sub_ID_electrode_localization_file_path = os.path.join(electrode_localization_file_path, sub_ID)
    electrode_localization_file_name = [f for f in sorted(os.listdir(sub_ID_electrode_localization_file_path))]
    electrode_localization = pd.read_csv(os.path.join(electrode_localization_file_path, sub_ID, electrode_localization_file_name[0]))
    outputfile_path = os.path.join(spectrogram_file_path_interpolated_avg, sub_ID)
    if not (os.path.isdir(outputfile_path)): os.mkdir(outputfile_path)
    for spectrogram_file in spectrogram_files:
        outputfile_name = spectrogram_file.replace('.pickle', '_avg.pickle')
        inputfile = os.path.join(inputfile_path, spectrogram_file)
        outputfile = os.path.join(outputfile_path,outputfile_name)
        print("reading file {0}".format(inputfile))
        with open(inputfile, 'rb') as f: spectrogram_interpolation, freqs, xnew, columns = pickle.load(f)
        num_GM = len(np.where(electrode_localization.iloc[:, 4] <= distance_considered_WM)[0])
        num_WM = len(np.where(electrode_localization.iloc[:, 4] > distance_considered_WM)[0])
        spectrogram_GM = np.zeros(
            shape=(spectrogram_interpolation.shape[0], spectrogram_interpolation.shape[1], num_GM))
        spectrogram_WM = np.zeros(
            shape=(spectrogram_interpolation.shape[0], spectrogram_interpolation.shape[1], num_WM))
        count = 0;
        count_GM = 0;
        count_WM = 0
        for e in columns:
            if any(np.array(electrode_localization.iloc[:, 0]) == e):#if any analyzed signals are not in the electrode localization file, dont' compute.
                loc = np.where(np.array(electrode_localization.iloc[:, 0]) == e)[0][0]
                dist_GM = electrode_localization.iloc[loc, 4]
                if dist_GM <= distance_considered_WM and dist_GM >=0:
                    spectrogram_GM[:, :, count_GM] = spectrogram_interpolation[:, :, count]
                    count_GM = count_GM + 1
                if dist_GM > distance_considered_WM:
                    spectrogram_WM[:, :, count_WM] = spectrogram_interpolation[:, :, count]
                    count_WM = count_WM + 1
            count = count + 1
        #some electrodes are have locations, but are not recorded on iEEG.org, so need to delete the part of the arrays
        spectrogram_GM = np.delete(spectrogram_GM, range(count_GM, num_GM), axis=2)
        spectrogram_WM = np.delete(spectrogram_WM, range(count_WM, num_WM), axis=2)
        spectrogram_GM_mean = np.mean(spectrogram_GM, axis=2)
        spectrogram_WM_mean = np.mean(spectrogram_WM, axis=2)
        print("saving file {0}\n\n".format(outputfile))
        with open(outputfile, 'wb') as f: pickle.dump([spectrogram_GM_mean, spectrogram_WM_mean, freqs, xnew, columns], f)
    total_num_electrodes_GM[0,cnt] = spectrogram_GM.shape[2]
    total_num_electrodes_WM[0,cnt] =  spectrogram_WM.shape[2]
    cnt = cnt + 1

np.sum(total_num_electrodes_GM)
np.sum(total_num_electrodes_WM)
#averaging by patient across each peri-ictal time
spectrogram_file_path_interpolated_avg_all = os.path.join(paper_path, "data_processed/spectrogram/montage/referential/filtered/interpolated_all_patients_combined")
ALL_data_GM = np.zeros(shape=(129, 200, 4 ))
ALL_data_WM = np.zeros(shape=(129, 200, 4 ))
for i in range(0,4):
    spectrogram_ALL_GM = np.zeros(shape=(129, 200, len(RID) ))
    spectrogram_ALL_WM = np.zeros(shape=(129, 200, len(RID)))
    count_sub = 0
    for sub_ID in RID:
        inputfile_path = os.path.join(spectrogram_file_path_interpolated_avg, sub_ID)
        spectrogram_file = [f for f in sorted(os.listdir(inputfile_path))][i]#assumes 0 = interictal, 1 = preictal, 2= ictal, 3 = postictal
        inputfile = os.path.join(inputfile_path, spectrogram_file)
        outputfile_name = "ALL_patient_spectrograms_averaged.pickle"
        outputfile = os.path.join(spectrogram_file_path_interpolated_avg_all, outputfile_name)
        if not (os.path.isdir(spectrogram_file_path_interpolated_avg_all)): os.mkdir(spectrogram_file_path_interpolated_avg_all)
        print("reading file {0}".format(inputfile))
        with open(inputfile, 'rb') as f: spectrogram_GM_mean, spectrogram_WM_mean, freqs, xnew, columns = pickle.load(f)
        spectrogram_ALL_GM[:, :, count_sub] = np.log10(spectrogram_GM_mean)
        spectrogram_ALL_WM[:, :, count_sub] = np.log10(spectrogram_WM_mean)
        count_sub = count_sub + 1
    print("{0}\n\n".format(i))
    spectrogram_ALL_GM_mean = np.mean(spectrogram_ALL_GM, axis=2)
    spectrogram_ALL_WM_mean = np.mean(spectrogram_ALL_WM, axis=2)
    ALL_data_GM[:,:,i] = spectrogram_ALL_GM_mean
    ALL_data_WM[:,:,i] = spectrogram_ALL_WM_mean

print("saving file {0}\n\n".format(outputfile))
with open(outputfile, 'wb') as f: pickle.dump([ALL_data_GM, ALL_data_WM, freqs, xnew], f)

np.sum(total_num_electrodes_GM)
np.sum(total_num_electrodes_WM)
