# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 19:21:54 2020

@author: andyr
"""

import pickle
import numpy as np
import os
import sys
from os.path import join as ospj
path = ospj("/media","arevell","sharedSSD","linux","papers","paper005") #Parent directory of project
#path = ospj("E:\\","linux","pastates","paper005") #Parent directory of project
sys.path.append(ospj(path, "seeg_GMvsWM", "code", "tools"))
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import signal
from scipy.stats import mannwhitneyu
import numpy as np
import copy
import json, codecs
import pingouin as pg
import statsmodels.api as sm
import bct
import networkx as nx
import matplotlib.colors
from colormap import rgb2hex #  pip install colormap;  pip install easydev
from scipy.stats import pearsonr, spearmanr
from mpl_toolkits.mplot3d import Axes3D
np.seterr(divide = 'ignore')
from sklearn.linear_model import LinearRegression
from statsmodels.gam.api import GLMGam, BSplines
import pygam #  pip install pygam
from scipy import stats
from sklearn.linear_model import TweedieRegressor
from pygam import LinearGAM
from imagingToolsRevell import printProgressBar


#%% Input/Output Paths and File names
fname_EEG_times = ospj( path, "data","data_raw","iEEG_times","EEG_times.xlsx")
fpath_electrode_localization = ospj( path, "data","data_processed","electrode_localization")
fname_atlases_csv = ospj( path, "data/data_raw/atlases/atlas_names.csv")
fpath_FC = ospj(path, "data","data_processed","connectivity_matrices","function")
fpath_SC = ospj(path, "data","data_processed","connectivity_matrices","structure")



fpath_figure = ospj(path, "seeg_GMvsWM","figures","sfc")


if not (os.path.isdir(fpath_figure)): os.makedirs(fpath_figure, exist_ok=True)

#% Load Study Meta Data
data = pd.read_excel(fname_EEG_times)    
atlases = pd.read_csv(fname_atlases_csv)    
#% Processing Meta Data: extracting sub-IDs

sub_IDs_unique =  np.unique(data.RID)[np.argsort( np.unique(data.RID, return_index=True)[1])]
#%%

#%%

DATA = [None] * len(sub_IDs_unique)
for pt in range(len(sub_IDs_unique)):
    DATA[pt] = [None] 
    for atlas in range(len(DATA[pt])):
        DATA[pt][atlas] = [None] * 3
        DATA[pt][atlas][2] = [None] * 4

descriptors = ["interictal","preictal","ictal","postictal"]

count_patient = 0
reference = ["CAR"]
ref = 0
for i in range(len(data)):
    #parsing data DataFrame to get iEEG information
    sub_ID = data.iloc[i].RID
    sub_RID = "sub-{0}".format(sub_ID)
    iEEG_filename = data.iloc[i].file
    ignore_electrodes = data.iloc[i].ignore_electrodes.split(",")
    start_time_usec = int(data.iloc[i].connectivity_start_time_seconds*1e6)
    stop_time_usec = int(data.iloc[i].connectivity_end_time_seconds*1e6)
    descriptor = data.iloc[i].descriptor
    print( "\n\n{0}: {1}".format(sub_ID,descriptor) )
    if descriptor == descriptors[0]: state = 0
    if descriptor == descriptors[1]: state = 1
    if descriptor == descriptors[2]: state = 2
    if descriptor == descriptors[3]: state = 3
    
    #Inputs and Outputs
    #input filename EEG
    fpath_FC_sub_ID = ospj(fpath_FC, sub_RID)
    fpath_electrode_localization_sub_ID =  ospj(fpath_electrode_localization, sub_RID)
    fname_FC = ospj(fpath_FC_sub_ID,reference[ref],  "sub-{0}_{1}_{2}_{3}_crossCorrelation.pickle".format(sub_ID, iEEG_filename, start_time_usec, stop_time_usec))
    fname_FC_metadata = ospj(fpath_FC_sub_ID,reference[ref],  "sub-{0}_{1}_{2}_{3}_metadata.json".format(sub_ID, iEEG_filename, start_time_usec, stop_time_usec))
    fpath_SC_sub_ID = ospj(fpath_SC, sub_RID)

    fname_electrode_localization = ospj(fpath_electrode_localization_sub_ID, "sub-{0}_electrode_localization.csv".format(sub_ID))

    
    
    #GET DATA
    #get localization and FC files
    with open(fname_FC, 'rb') as f: FC = pickle.load(f)
    electrode_localization = pd.read_csv(fname_electrode_localization)
    metadata = json.loads(codecs.open(fname_FC_metadata, 'r', encoding='utf-8').read())
    
    #get electrode names in FC and localization files
    electrode_names_FC = np.asarray(metadata["channels"])
    electrode_names_localization = np.array(electrode_localization["electrode_name"])
    
    
    #Preprocessing files
    #find electrodes in both files (i.e there are electrodes in localization files but not in FC files, and there are electrodes in FC files but not localized)
    electrode_names_intersect, electrode_names_FC_ind, electrode_names_localization_ind =  np.intersect1d(electrode_names_FC, electrode_names_localization, return_indices = True)

    #Equalizing discrepancies in localization files and FC files
    #removing electrodes in localization files not in FC
    electrode_localization_intersect = electrode_localization.iloc[electrode_names_localization_ind].reset_index(drop=True)

    
    #removing electrodes in FC not in localization files
    FC_intersect = copy.deepcopy(FC)
    FC_intersect = FC_intersect[:,electrode_names_FC_ind[:,None], electrode_names_FC_ind[None,:], :] 
    electrode_names_FC_intersect = electrode_names_FC[electrode_names_FC_ind]

    
        
    #Only GM and WM tissues considered (not CSF=1, or outside brain = 0)
    labels = np.array(electrode_localization_intersect["Tissue_segmentation_region_number"])
    labels_gmwm_ind = np.where(labels >= 2)[0] #Only GM and WM tissues considered
    #removing electrode localization electrodes not in GM or WM
    electrode_localization_intersect_GMWM = electrode_localization_intersect.iloc[labels_gmwm_ind].reset_index(drop=True)
    #removing FC electrodes not in GM or WM
    FC_intersect_GMWM = copy.deepcopy(FC_intersect)
    FC_intersect_GMWM = FC_intersect_GMWM[:,labels_gmwm_ind[:,None], labels_gmwm_ind[None,:], :] 
    electrode_names_FC_intersect_GMWM = np.array(electrode_localization_intersect_GMWM["electrode_name"])
    
    #averaging FC
    FC_intersect_GMWM_mean = [None] * len(FC_intersect_GMWM)
    for f in range(len(FC_intersect_GMWM_mean)):
        FC_intersect_GMWM_mean[f] = np.nanmean(FC_intersect_GMWM[f], axis=2)    
    np.allclose(FC_intersect_GMWM_mean[f], FC_intersect_GMWM_mean[f].T, rtol=1e-05, atol=1e-08)#check if symmetric
    



    #calculating distances matrices to see if GM-GM distances are different from WM-WM distances
    #takes a while to compute, so only need to do it once
    if descriptor == descriptors[0]:
        
        electrode_localization_intersect_GMWM
        distance_matrix = np.zeros(shape = (len(electrode_localization_intersect_GMWM),len(electrode_localization_intersect_GMWM) ))
        for e1 in range(len(electrode_localization_intersect_GMWM)):
            for e2 in range(len(electrode_localization_intersect_GMWM)):
                p1 = [electrode_localization_intersect_GMWM.iloc[e1]["x_coordinate"], electrode_localization_intersect_GMWM.iloc[e1]["y_coordinate"], electrode_localization_intersect_GMWM.iloc[e1]["z_coordinate"]]
                p2 = [electrode_localization_intersect_GMWM.iloc[e2]["x_coordinate"], electrode_localization_intersect_GMWM.iloc[e2]["y_coordinate"], electrode_localization_intersect_GMWM.iloc[e2]["z_coordinate"]]
                distance_matrix[e1,e2] = np.sqrt(  np.sum((np.array(p1)-np.array(p2))**2, axis=0)    )

    #Structure
    if (os.path.exists(fpath_SC_sub_ID)): #If structure exists
        for a in range(len(atlases)):
            atlas_fname = os.path.splitext(os.path.splitext(  atlases.iloc[a]["atlas_filename"] )[0])[0]
            atlas_name = os.path.splitext(os.path.splitext(  atlases.iloc[a]["atlas_name"] )[0])[0]
            print(f"{sub_ID}; {atlas_name}")
            basename = f"sub-{sub_ID}_preop3T_{atlas_fname}"
            fname_SC =  ospj( fpath_SC_sub_ID, f"sub-{sub_ID}.{basename}.count.pass.connectogram.txt") 
            
            
            SC = pd.read_table(fname_SC, header=None, dtype=object)
            #cleaning up structural data 
            SC = SC.drop([0,1], axis=1)
            SC = SC.drop([0], axis=0)
            SC = SC.iloc[:, :-1]
            SC_regionNames = np.array([e[len(basename)+1:] for e in  np.array(SC.iloc[0])])
            SC = np.array(SC.iloc[1:, :]).astype('float64')  #finally turn into numpy array
            
            
            #reorder structure
            #remove electrodes not in FC 
    
            #find electrodes in both files (i.e there are electrodes in localization files but not in FC files, and there are electrodes in FC files but not localized)
            electrode_names_intersect_SC, electrode_names_FC_ind_SC, electrode_names_localization_ind_SC =  np.intersect1d(electrode_names_FC_intersect_GMWM, structure_electrode_names, return_indices = True)
            #Equalizing discrepancies in SC files and FC files
            #removing electrodes in SC, but not in FC files
            SC = copy.deepcopy(structure)
            SC = SC[electrode_names_localization_ind_SC[:,None], electrode_names_localization_ind_SC[None,:]] 
            
            electrode_names_SC = structure_electrode_names[electrode_names_localization_ind_SC]
            
            print(electrode_names_SC == electrode_names_FC_intersect_GMWM)
            

 
            DATA[count_patient][a][0] = [electrode_localization_intersect_GMWM, distance_matrix ]
            DATA[count_patient][a][1] = [SC, distance_matrix   ] 
            DATA[count_patient][a][2][state] = [FC_intersect_GMWM, FC_intersect_GMWM_mean]

    if descriptor == descriptors[3]: count_patient = count_patient + 1




frequencies = np.hstack(np.asarray(order_of_matrices_in_pickle_file))

#%%
    
    
    
    
    