#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 13:54:02 2021

@author: arevell
"""

#%% Path
path = "/media/arevell/sharedSSD/linux/papers/paper005" 

# Imports

import sys
import os
import pandas as pd
import numpy as np
from os.path import join as ospj
import seaborn as sns
import matplotlib.pyplot as plt
import time
from os.path import splitext as split
from matplotlib import colors
from colormap import rgb2hex 
from scipy.stats import mannwhitneyu
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
#% Custom Imports

sys.path.append(ospj(path, "seeg_GMvsWM", "code", "tools"))
import network_measures

#% Input/Output Files and paths

fname_diffusion_imaging = ospj( path, "data/data_raw/iEEG_times/diffusion_imaging.csv")
fname_atlases_csv = ospj( path, "data/data_raw/atlases/atlas_names.csv")
fpath_connectivity = ospj( path, "data/data_processed/connectivity_matrices/structure")
fpath_networkMeasures = ospj( path, "data/data_processed/network_measures/structure")
fname_network_measures = ospj( fpath_networkMeasures, "network_measures_structure.csv")
fname_meanVols = ospj(path, "data/data_processed/atlas_morphology/volumes_and_sphericity.csv")
fpath_figure = ospj(path, "figures", "network")
if not (os.path.isdir(fpath_figure)): os.makedirs(fpath_figure, exist_ok=True)

#% Load Data

data = pd.read_csv(fname_diffusion_imaging) 
atlases = pd.read_csv(fname_atlases_csv)    
measures = pd.read_csv(fname_network_measures)    
vols = pd.read_csv(fname_meanVols)    

#renmae atlases in measures
for i in range(len(measures)):
    measures["atlas"].loc[ i ] = np.array(atlases["atlas_name"][measures["atlas"].loc[ i ]  == np.array([split(file_name)[0] for file_name in [split(file_name)[0] for file_name in atlases.atlas_filename]])     ])[0]                                     
    

#change name of all random atlases
for i in range(len(measures)):
    if measures["atlas"][i][0:11] == "RandomAtlas":
        measures["atlas"].loc[ i ] = measures["atlas"].loc[i][0:18]

for i in range(len(vols)):
    if vols["atlas_name"][i][0:11] == "RandomAtlas":
        vols["atlas_name"].loc[ i ] = vols["atlas_name"].loc[i][0:18]

vols =  vols.groupby(['atlas_name'], as_index=False).mean()
    
sub_IDs_unique = np.unique(data.RID)[np.argsort( np.unique(data.RID, return_index=True)[1])]

#% Add in whether sub_ID is a control or note
sub_IDs = np.array(measures.RID)
patient_types = []

for i in range(len(sub_IDs)):
    patient_types.append( np.array(data.type[np.where(sub_IDs[i] == data.RID )[0]])[0]    ) 
    if patient_types[i] == "SEEG" or patient_types[i] == "ECoG" or patient_types[i] == "EcoG":
        patient_types[i] = "patient"
measures["type"] =    patient_types  

#% add in volumes
atlases_IDs = np.array(measures.atlas)
vols_IDs = []

for i in range(len(atlases_IDs)):
    vols_IDs.append( np.log10(  np.array(   vols.volume_voxels[atlases_IDs[i] == vols.atlas_name ])[0])       ) 

measures["vols"] =    vols_IDs  



#% averages
df_all =  measures.groupby(["atlas"], as_index=False).mean()
df_bysubject =  measures.groupby(["type", "atlas"], as_index=False).mean()
 
#%% Plotting parameters
subset = ["AAL2", "AAL600", "Craddock_200", "Craddock_400", "DKT31_OASIS", "MMP", "Talairach", "Yeo_17_liberal", "Yeo_7_liberal",
          "Schaefer_17_100","Schaefer_17_200","Schaefer_17_300","Schaefer_17_400","Schaefer_17_500","Schaefer_17_600","Schaefer_17_700","Schaefer_17_800","Schaefer_17_900","Schaefer_17_1000",
          "RandomAtlas0000010","RandomAtlas0000030","RandomAtlas0000050","RandomAtlas0000075","RandomAtlas0000100","RandomAtlas0000200","RandomAtlas0000300","RandomAtlas0000400",
          "RandomAtlas0000500","RandomAtlas0000750","RandomAtlas0001000","RandomAtlas0002000","RandomAtlas0005000","RandomAtlas0010000"]
subset_legend = ["AAL2", "AAL600", "Craddock 400", "Craddock 200", "DKT31 OASIS", "MMP","Talairach", "Yeo 17 liberal", "Yeo 7 liberal", "Schaefer", "Random Atlas"]
legend_order1 = [6, 5, 1, 2, 3, 4, 0, 7, 8, 9 , 10 ]

legend_order = [6, 0, 5, 7, 1, 8, 2, 9, 3, 10, 4]
#% Assign subsets
df_all_subset = df_all.loc[df_all['atlas'].isin(subset)].reset_index(drop = True)

df_bysubject_ctrl = df_bysubject.loc[df_bysubject['type'].isin(["control"])].reset_index(drop = True)
df_bysubject_ctrl_subset = df_bysubject_ctrl.loc[df_bysubject_ctrl['atlas'].isin(subset)].reset_index(drop = True)
df_bysubject_pt = df_bysubject.loc[df_bysubject['type'].isin(["patient"])].reset_index(drop = True)
df_bysubject_pt_subset = df_bysubject_pt.loc[df_bysubject_pt['atlas'].isin(subset)].reset_index(drop = True)


subset_ctrl_vs_pt = ["AAL2"]
df_ctrl_vs_pt = measures.loc[measures['atlas'].isin(subset_ctrl_vs_pt)].reset_index(drop = True)


#%
random_color = "#000000cc"
schaefer_color = "#33333399"
colors_list = ["#8256b3", "#b3b156", "#6ba069","#6ba069","#5687b3","#b38256","#b35659","#b156b3","#b156b3",
          schaefer_color,schaefer_color,schaefer_color,schaefer_color,schaefer_color,schaefer_color,schaefer_color,schaefer_color,schaefer_color,schaefer_color,
          random_color, random_color, random_color, random_color, random_color, random_color,random_color, random_color,random_color, random_color,random_color, random_color,random_color, random_color]
colors_legend = ["#8256b3", "#b3b156", "#6ba069","#6ba069","#5687b3","#b38256","#b35659","#b156b3","#b156b3",schaefer_color, random_color]

markers = { "AAL2": "v", 
           "AAL600": "v", 
           "Craddock_200": ",", 
           "Craddock_400": "o", 
           "DKT31_OASIS": "v", 
           "MMP": "v",
           "Talairach": "v", 
           "Yeo_17_liberal": "v", 
           "Yeo_7_liberal": "^",
          "Schaefer_17_100": "o","Schaefer_17_200": "o","Schaefer_17_300": "o","Schaefer_17_400": "o","Schaefer_17_500": "o","Schaefer_17_600": "o","Schaefer_17_700": "o","Schaefer_17_800": "o","Schaefer_17_900": "o","Schaefer_17_1000": "o",
          "RandomAtlas0000010": "o","RandomAtlas0000030": "o","RandomAtlas0000050": "o","RandomAtlas0000075": "o","RandomAtlas0000100": "o","RandomAtlas0000200": "o","RandomAtlas0000300": "o","RandomAtlas0000400": "o",
          "RandomAtlas0000500": "o","RandomAtlas0000750": "o","RandomAtlas0001000": "o","RandomAtlas0002000": "o","RandomAtlas0005000": "o","RandomAtlas0010000": "o"  }
markers_legend = ["v", "v" , "o", "s", "v", "v","v","v","^","o","o"]


norm = 100
small = 40
size = [norm,norm,norm,norm,norm,norm,norm,norm,norm,
        small,small,small,small,small,small,small,small,small,small,
        small,small,small,small,small,small,small,small,small,small,small,small,small,small]


#Assign colors
order =  np.intersect1d(df_all_subset["atlas"],  subset ,return_indices=True)[2] 
colors_plot = np.array(colors_list)[  order  ]
size =  np.array(size)[  order  ]



#%% plot figure for paper

fig = plt.figure(constrained_layout=False, figsize=(12.5,7) , dpi=300)
left = 0.03; right = 0.99
bottom1 = 0.20; bottom2 = 0.00; 
top1 = 0.95; top2 = 0.12
gs1 = fig.add_gridspec(nrows = 2, ncols = 3, left = left, right = right, bottom = bottom1, top = top1)
gs2 = fig.add_gridspec(nrows = 1, ncols = 1, left = 0, right = right, bottom = bottom2, top = top2)
gs = [gs1, gs2]
ax = [None] * 2
ax[0] = [ fig.add_subplot(gs[0][0, 0]), fig.add_subplot(gs[0][0, 1]),fig.add_subplot(gs[0][0, 2]),fig.add_subplot(gs[0][1, 0]),fig.add_subplot(gs[0][1, 1]),fig.add_subplot(gs[0][1, 2]) ]
ax[1] = [ fig.add_subplot(gs[1][0, 0])]
network_names = ["Density", "degree_mean", "clustering_coefficient_mean", "characteristic_path_length", "small_worldness"]
network_labels = ["Density", "Mean Degree", "Mean Clustering Coefficient", "Characteristic Path Length", "Small Worldness"]
for i in range(5):
    sns.scatterplot(x = "vols", y = network_names[i], data = df_all_subset, ax = ax[0][i], 
                    hue = "atlas", palette = colors_plot.tolist(), style = "atlas", s= size,  legend = False, markers = markers)
    ax[0][i].set_xlabel('')
    ax[0][i].set_ylabel("")
    ax[0][i].set_title("")

    
ax[0][0].set_ylim([0.16, 1.05]);ax[0][1].set_ylim([-100, 1700]);ax[0][2].set_ylim([0.58, 1.02]);ax[0][3].set_ylim([0.95, 2.3]);ax[0][4].set_ylim([0.9, 3.3]);

sns.violinplot(  data = df_ctrl_vs_pt, x = "type", y= "Density", ax = ax[0][5] , palette = ["#6899c7", "#c0d2e3"],linewidth=1, inner="quartile")
sns.swarmplot(data = df_ctrl_vs_pt, x = "type", y= "Density", ax = ax[0][5] , palette = ["#000000", "#333333"])
ax[0][5].set_ylabel(network_labels[0])
ax[0][5].set_xlabel("")
ax[0][5].set_ylabel("")

#fig.suptitle('N = 41 (13 controls, 28 patients)', fontsize = 13)
fig.text(0.5, bottom2 + (bottom1-bottom2)*0.6, 'Volume ($log_{10}$  $mm^3$)', ha='center', size=13)
#write p value
pval = mannwhitneyu(     df_ctrl_vs_pt.loc[df_ctrl_vs_pt['type'].isin(["control"])]["Density"] , df_ctrl_vs_pt.loc[df_ctrl_vs_pt['type'].isin(["patient"])]["Density"] )[1]

#legend
legend_elements = []

for i in range(len(subset_legend)):
    legend_elements.append( Line2D(   [0], [0], marker=  np.array(markers_legend)[legend_order][i], color=np.array(colors_legend)[legend_order][i], label=np.array(subset_legend)[legend_order][i], markersize=8, linestyle = 'None')   )


ax[1][0].legend(handles=legend_elements, loc='center', ncol=6 , frameon=False,fontsize=13)
ax[1][0].axis('off')

text_size = 15
ax[0][0].text(0.05, 0.95, network_labels[0], horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][0].transAxes)
ax[0][1].text(0.95, 0.95, network_labels[1], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][1].transAxes)
ax[0][2].text(0.05, 0.95, "Mean Clustering\nCoefficient", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][2].transAxes)
ax[0][3].text(0.50, 0.95, "Characteristic \nPath Length", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][3].transAxes)
ax[0][4].text(0.95, 0.95, network_labels[4], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][4].transAxes)
ax[0][5].text(0.5, 0.1, f"Atlas: {subset_ctrl_vs_pt[0]}", horizontalalignment='center', transform=ax[0][5].transAxes, size = 15 )
ax[0][5].text(0.5, 0.4, "p > 0.05", horizontalalignment='center', transform=ax[0][5].transAxes, size = 15 )

ax[0][1].text(0.5, 1.05, "N = 41 (13 controls, 28 patients)", horizontalalignment='center', transform=ax[0][1].transAxes, size = text_size )
fname_fig = ospj(fpath_figure,"network_measures_all.pdf")
plt.savefig(fname_fig  )


#%%Supplemental: controls and patients separated
dat = [df_bysubject_ctrl_subset, df_bysubject_pt_subset]
titles = ['CONTROLS (N = 13)', 'PATIENTS (N = 28)']
fnames_supp = ["network_measures_controls.pdf", "network_measures_patients.pdf"]
for p in range(2):
    fig = plt.figure(constrained_layout=False, figsize=(12.5,7) , dpi=300)
    left = 0.03; right = 0.99
    bottom1 = 0.09
    top1 = 0.95
    gs1 = fig.add_gridspec(nrows = 2, ncols = 3, left = left, right = right, bottom = bottom1, top = top1)
    gs = [gs1]
    ax = [None] 
    ax[0] = [ fig.add_subplot(gs[0][0, 0]), fig.add_subplot(gs[0][0, 1]),fig.add_subplot(gs[0][0, 2]),fig.add_subplot(gs[0][1, 0]),fig.add_subplot(gs[0][1, 1]),fig.add_subplot(gs[0][1, 2]) ]
    for i in range(5):
        sns.scatterplot(x = "vols", y = network_names[i], data = dat[p], ax = ax[0][i], 
                        hue = "atlas", palette = colors_plot.tolist(), style = "atlas", s= size,  legend = False, markers = markers)
        ax[0][i].set_xlabel('')
        ax[0][i].set_ylabel("")
        ax[0][i].set_title("")
        
    legend_elements = []
    for i in range(len(subset_legend)):
        legend_elements.append( Line2D(   [0], [0], marker=  np.array(markers_legend)[legend_order1][i], color=np.array(colors_legend)[legend_order1][i], label=np.array(subset_legend)[legend_order1][i], markersize=10, linestyle = 'None')   )
        
        
    ax[0][0].set_ylim([0.16, 1.05]);ax[0][1].set_ylim([-100, 1700]);ax[0][2].set_ylim([0.58, 1.02]);ax[0][3].set_ylim([0.95, 2.3]);ax[0][4].set_ylim([0.9, 3.3]);
    ax[0][5].legend(handles=legend_elements, loc='center', ncol=2 , frameon=False)
    ax[0][5].axis('off')
    
    #fig.suptitle(titles[p])
    fig.text(0.5,0.02, 'Volume ($log_{10}$  $mm^3$)', ha='center')

    ax[0][0].text(0.05, 0.95, network_labels[0], horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][0].transAxes)
    ax[0][1].text(0.95, 0.95, network_labels[1], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][1].transAxes)
    ax[0][2].text(0.05, 0.95, "Mean Clustering\nCoefficient", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][2].transAxes)
    ax[0][3].text(0.50, 0.95, "Characteristic \nPath Length", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][3].transAxes)
    ax[0][4].text(0.95, 0.95, network_labels[4], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][4].transAxes)
    ax[0][1].text(0.5, 1.05, titles[p], horizontalalignment='center', transform=ax[0][1].transAxes, size = text_size )
    fname_fig = ospj(fpath_figure,fnames_supp[p])
    plt.savefig(fname_fig  )

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
#%% Supplemental: Other Atlases FIRST
colormap = colors.LinearSegmentedColormap.from_list("", ["#c43c3c","#d6a677","#d6d677","#a7d677","#77a7d6","#7778d6", "#d677a7"]) #rainbow ROYGBIV
N = 10
#Colorbar - Make one unique color for each electrode stick base on the rainbow above
colors_rainbow = []
for c in np.linspace(start=0, stop=255, num=N).astype(int):
    colors_rainbow.append( rgb2hex(int(colormap(c)[0]*255), int(colormap(c)[1]*255), int(colormap(c)[2]*255))  )
random_color = "#000000cc"
schaefer_color = "#33333399"
colors_rainbow.append(schaefer_color)
colors_rainbow.append(random_color)
rainbow = np.array(colors_rainbow)
#%
subset = ["AAL1", "AAL3", "AAL_JHU_combined", "AICHA", "BN", "Brodmann", "Harvard_Oxford_cortical_and_subcortical", "Gordon_Petersen", "Hammersmith", "JHU",
          "Schaefer_17_100","Schaefer_17_200","Schaefer_17_300","Schaefer_17_400","Schaefer_17_500","Schaefer_17_600","Schaefer_17_700","Schaefer_17_800","Schaefer_17_900","Schaefer_17_1000",
          "RandomAtlas0000010","RandomAtlas0000030","RandomAtlas0000050","RandomAtlas0000075","RandomAtlas0000100","RandomAtlas0000200","RandomAtlas0000300","RandomAtlas0000400",
          "RandomAtlas0000500","RandomAtlas0000750","RandomAtlas0001000","RandomAtlas0002000","RandomAtlas0005000","RandomAtlas0010000"]
subset_legend = ["AAL1", "AAL3", "AAL-JHU", "AICHA", "BN", "Brodmann", "HO cort + sub", "Gordon_Petersen", "Hammersmith", "JHU", "Schaefer", "Random Atlas"]
legend_order = [7, 3, 4, 6, 1, 2, 0, 8, 5, 9 , 10, 11 ]
#legend_order = [6, 0, 5, 7, 1, 8, 2, 9, 3, 10, 4, 11]
#% Assign subsets
df_all_subset = df_all.loc[df_all['atlas'].isin(subset)].reset_index(drop = True)

df_bysubject_ctrl = df_bysubject.loc[df_bysubject['type'].isin(["control"])].reset_index(drop = True)
df_bysubject_ctrl_subset = df_bysubject_ctrl.loc[df_bysubject_ctrl['atlas'].isin(subset)].reset_index(drop = True)
df_bysubject_pt = df_bysubject.loc[df_bysubject['type'].isin(["patient"])].reset_index(drop = True)
df_bysubject_pt_subset = df_bysubject_pt.loc[df_bysubject_pt['atlas'].isin(subset)].reset_index(drop = True)


#%
colors_list = [rainbow[6], rainbow[4], rainbow[5], rainbow[1], rainbow[2], rainbow[8], rainbow[3], rainbow[0], rainbow[7], rainbow[9], 
          schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color,
          random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color]
colors_legend = [rainbow[6], rainbow[4], rainbow[5], rainbow[1], rainbow[2], rainbow[8], rainbow[3], rainbow[0], rainbow[7], rainbow[9],  schaefer_color, random_color]

markers = { "AAL1": "v", 
           "AAL3": "v", 
           "AAL_JHU_combined": "v", 
           "AICHA": "v", 
           "BN": "v", 
           "Brodmann": "v",
           "Harvard_Oxford_cortical_and_subcortical": "v", 
           "Gordon_Petersen": "v", 
           "Hammersmith": "v",
           "JHU": "v",
          "Schaefer_17_100": "o","Schaefer_17_200": "o","Schaefer_17_300": "o","Schaefer_17_400": "o","Schaefer_17_500": "o","Schaefer_17_600": "o","Schaefer_17_700": "o","Schaefer_17_800": "o","Schaefer_17_900": "o","Schaefer_17_1000": "o",
          "RandomAtlas0000010": "o","RandomAtlas0000030": "o","RandomAtlas0000050": "o","RandomAtlas0000075": "o","RandomAtlas0000100": "o","RandomAtlas0000200": "o","RandomAtlas0000300": "o","RandomAtlas0000400": "o",
          "RandomAtlas0000500": "o","RandomAtlas0000750": "o","RandomAtlas0001000": "o","RandomAtlas0002000": "o","RandomAtlas0005000": "o","RandomAtlas0010000": "o"  }
markers_legend = ["v", "v" ,"v", "v", "v", "v","v","v","v", "v", "o","o"]


norm = 100
small = 40
size = [norm,norm,norm,norm,norm,norm,norm,norm,norm,norm,
        small,small,small,small,small,small,small,small,small,small,
        small,small,small,small,small,small,small,small,small,small,small,small,small,small]


#Assign colors
order =  np.intersect1d(df_all_subset["atlas"],  subset ,return_indices=True)[2] 
colors_plot = np.array(colors_list)[  order  ]
size =  np.array(size)[  order  ]

#%
#plot
fig = plt.figure(constrained_layout=False, figsize=(12.5,7) , dpi=300)
left = 0.03; right = 0.99
bottom1 = 0.09
top1 = 0.95; top2 = 0.15
gs1 = fig.add_gridspec(nrows = 2, ncols = 3, left = left, right = right, bottom = bottom1, top = top1)
gs = [gs1]
ax = [None] 
ax[0] = [ fig.add_subplot(gs[0][0, 0]), fig.add_subplot(gs[0][0, 1]),fig.add_subplot(gs[0][0, 2]),fig.add_subplot(gs[0][1, 0]),fig.add_subplot(gs[0][1, 1]),fig.add_subplot(gs[0][1, 2]) ]

network_names = ["Density", "degree_mean", "clustering_coefficient_mean", "characteristic_path_length", "small_worldness"]
network_labels = ["Density", "Mean Degree", "Mean Clustering Coefficient", "Characteristic Path Length", "Small Worldness"]
for i in range(5):
    sns.scatterplot(x = "vols", y = network_names[i], data = df_all_subset, ax = ax[0][i], 
                    hue = "atlas", palette = colors_plot.tolist(), style = "atlas", s= size,  legend = False, markers = markers)
    ax[0][i].set_xlabel('')
    ax[0][i].set_ylabel("")
    ax[0][i].set_title("")
ax[0][0].set_ylim([0.16, 1.05]);ax[0][1].set_ylim([-100, 1700]);ax[0][2].set_ylim([0.58, 1.02]);ax[0][3].set_ylim([0.95, 2.3]);ax[0][4].set_ylim([0.9, 3.3]);


#fig.suptitle('Network Measures; N = 41 (13 controls, 28 patients)')
fig.text(0.5, 0.02, 'Volume ($log_{10}$  $mm^3$)', ha='center')
#write p value

#legend
legend_elements = []

for i in range(len(subset_legend)):
    legend_elements.append( Line2D(   [0], [0], marker=  np.array(markers_legend)[legend_order][i], color=np.array(colors_legend)[legend_order][i], label=np.array(subset_legend)[legend_order][i], markersize=8, linestyle = 'None')   )


ax[0][5].legend(handles=legend_elements, loc='center', ncol=2 , frameon=False)
ax[0][5].axis('off')

ax[0][0].text(0.05, 0.95, network_labels[0], horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][0].transAxes)
ax[0][1].text(0.95, 0.95, network_labels[1], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][1].transAxes)
ax[0][2].text(0.05, 0.95, "Mean Clustering\nCoefficient", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][2].transAxes)
ax[0][3].text(0.50, 0.95, "Characteristic \nPath Length", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][3].transAxes)
ax[0][4].text(0.95, 0.95, network_labels[4], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][4].transAxes)
ax[0][1].text(0.5, 1.05, "Remaining Atlases Plot 1", horizontalalignment='center', transform=ax[0][1].transAxes, size = text_size )


fname_fig = ospj(fpath_figure,"supplemental_all_atlases_FIRST.pdf")
plt.savefig(fname_fig, figsize=(12.5,7)  )

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
#%% Supplemental: Other Atlases SECOND

subset = ["CerebrA", "HO_cortical_nonsymmetric", "HO_subcortical_only", "Hammersmith", "Harvard_Oxford_cortical_only", "JHU_tracts", "Julich", "Yeo_17", "Yeo_7", "Harvard_Oxford_cortical_and_subcortical",
          "Schaefer_17_100","Schaefer_17_200","Schaefer_17_300","Schaefer_17_400","Schaefer_17_500","Schaefer_17_600","Schaefer_17_700","Schaefer_17_800","Schaefer_17_900","Schaefer_17_1000",
          "RandomAtlas0000010","RandomAtlas0000030","RandomAtlas0000050","RandomAtlas0000075","RandomAtlas0000100","RandomAtlas0000200","RandomAtlas0000300","RandomAtlas0000400",
          "RandomAtlas0000500","RandomAtlas0000750","RandomAtlas0001000","RandomAtlas0002000","RandomAtlas0005000","RandomAtlas0010000"]
subset_legend = ["CerebrA", "HO cort nonsymmetric", "HO sub only", "Hammersmith", "HO cort only", "JHU tracts", "Julich", "Yeo 17", "Yeo 7", "HO cort + sub", "Schaefer", "Random Atlas"]
legend_order = [5, 6, 2, 9, 1, 0, 3, 4, 7, 8 , 10, 11 ]
loa = np.array(legend_order)
#legend_order = [6, 0, 5, 7, 1, 8, 2, 9, 3, 10, 4, 11]
#% Assign subsets
df_all_subset = df_all.loc[df_all['atlas'].isin(subset)].reset_index(drop = True)

df_bysubject_ctrl = df_bysubject.loc[df_bysubject['type'].isin(["control"])].reset_index(drop = True)
df_bysubject_ctrl_subset = df_bysubject_ctrl.loc[df_bysubject_ctrl['atlas'].isin(subset)].reset_index(drop = True)
df_bysubject_pt = df_bysubject.loc[df_bysubject['type'].isin(["patient"])].reset_index(drop = True)
df_bysubject_pt_subset = df_bysubject_pt.loc[df_bysubject_pt['atlas'].isin(subset)].reset_index(drop = True)


#%
colors_list = [rainbow[np.where(loa == 0)[0][0]], rainbow[np.where(loa == 1)[0][0]], rainbow[np.where(loa == 2)[0][0]], rainbow[np.where(loa == 3)[0][0]], 
               rainbow[np.where(loa == 4)[0][0]], rainbow[np.where(loa == 5)[0][0]], rainbow[np.where(loa == 6)[0][0]], rainbow[np.where(loa == 7)[0][0]], 
               rainbow[np.where(loa == 8)[0][0]], rainbow[np.where(loa == 9)[0][0]], 
          schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color, schaefer_color,
          random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color, random_color]
colors_legend = [rainbow[np.where(loa == 0)[0][0]], rainbow[np.where(loa == 1)[0][0]], rainbow[np.where(loa == 2)[0][0]], rainbow[np.where(loa == 3)[0][0]], 
               rainbow[np.where(loa == 4)[0][0]], rainbow[np.where(loa == 5)[0][0]], rainbow[np.where(loa == 6)[0][0]], rainbow[np.where(loa == 7)[0][0]], 
               rainbow[np.where(loa == 8)[0][0]], rainbow[np.where(loa == 9)[0][0]],   schaefer_color, random_color]

markers = { "CerebrA": "v", 
           "HO_cortical_nonsymmetric": "v", 
           "HO_subcortical_only": "v", 
           "Hammersmith": "v", 
           "Harvard_Oxford_cortical_only": "v", 
           "JHU_tracts": "v",
           "Julich": "v", 
           "Yeo_17": "v", 
           "Yeo_7": "v",
           "Harvard_Oxford_cortical_and_subcortical": "v",
          "Schaefer_17_100": "o","Schaefer_17_200": "o","Schaefer_17_300": "o","Schaefer_17_400": "o","Schaefer_17_500": "o","Schaefer_17_600": "o","Schaefer_17_700": "o","Schaefer_17_800": "o","Schaefer_17_900": "o","Schaefer_17_1000": "o",
          "RandomAtlas0000010": "o","RandomAtlas0000030": "o","RandomAtlas0000050": "o","RandomAtlas0000075": "o","RandomAtlas0000100": "o","RandomAtlas0000200": "o","RandomAtlas0000300": "o","RandomAtlas0000400": "o",
          "RandomAtlas0000500": "o","RandomAtlas0000750": "o","RandomAtlas0001000": "o","RandomAtlas0002000": "o","RandomAtlas0005000": "o","RandomAtlas0010000": "o"  }
markers_legend = ["v", "v" ,"v", "v", "v", "v","v","v","v", "v", "o","o"]


norm = 100
small = 40
size = [norm,norm,norm,norm,norm,norm,norm,norm,norm,norm,
        small,small,small,small,small,small,small,small,small,small,
        small,small,small,small,small,small,small,small,small,small,small,small,small,small]


#Assign colors
order =  np.intersect1d(df_all_subset["atlas"],  subset ,return_indices=True)[2] 
colors_plot = np.array(colors_list)[  order  ]
size =  np.array(size)[  order  ]

#%
#plot
fig = plt.figure(constrained_layout=False, figsize=(12.5,7) , dpi=300)
left = 0.03; right = 0.99
bottom1 = 0.09
top1 = 0.95; top2 = 0.15
gs1 = fig.add_gridspec(nrows = 2, ncols = 3, left = left, right = right, bottom = bottom1, top = top1)
gs = [gs1]
ax = [None] 
ax[0] = [ fig.add_subplot(gs[0][0, 0]), fig.add_subplot(gs[0][0, 1]),fig.add_subplot(gs[0][0, 2]),fig.add_subplot(gs[0][1, 0]),fig.add_subplot(gs[0][1, 1]),fig.add_subplot(gs[0][1, 2]) ]

network_names = ["Density", "degree_mean", "clustering_coefficient_mean", "characteristic_path_length", "small_worldness"]
network_labels = ["Density", "Mean Degree", "Mean Clustering Coefficient", "Characteristic Path Length", "Small Worldness"]
for i in range(5):
    sns.scatterplot(x = "vols", y = network_names[i], data = df_all_subset, ax = ax[0][i], 
                    hue = "atlas", palette = colors_plot.tolist(), style = "atlas", s= size,  legend = False, markers = markers)
    ax[0][i].set_xlabel('')
    ax[0][i].set_ylabel("")
    ax[0][i].set_title("")
ax[0][0].set_ylim([0.16, 1.05]);ax[0][1].set_ylim([-100, 1700]);ax[0][2].set_ylim([0.58, 1.02]);ax[0][3].set_ylim([0.95, 2.3]);ax[0][4].set_ylim([0.9, 3.3]);


#fig.suptitle('Network Measures; N = 41 (13 controls, 28 patients)')
fig.text(0.5, 0.02, 'Volume ($log_{10}$  $mm^3$)', ha='center')
#write p value

#legend
legend_elements = []

for i in range(len(subset_legend)):
    legend_elements.append( Line2D(   [0], [0], marker=  np.array(markers_legend)[legend_order][i], color=np.array(colors_legend)[legend_order][i], label=np.array(subset_legend)[legend_order][i], markersize=8, linestyle = 'None')   )


ax[0][5].legend(handles=legend_elements, loc='center', ncol=2 , frameon=False)
ax[0][5].axis('off')



ax[0][0].text(0.05, 0.95, network_labels[0], horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][0].transAxes)
ax[0][1].text(0.95, 0.95, network_labels[1], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][1].transAxes)
ax[0][2].text(0.05, 0.95, "Mean Clustering\nCoefficient", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][2].transAxes)
ax[0][3].text(0.50, 0.95, "Characteristic \nPath Length", horizontalalignment='left', verticalalignment='top',size = text_size , transform=ax[0][3].transAxes)
ax[0][4].text(0.95, 0.95, network_labels[4], horizontalalignment='right', verticalalignment='top',size = text_size , transform=ax[0][4].transAxes)
ax[0][1].text(0.5, 1.05, "Remaining Atlases Plot 2", horizontalalignment='center', transform=ax[0][1].transAxes, size = text_size )

fname_fig = ospj(fpath_figure,"supplemental_all_atlases_SECOND.pdf")
plt.savefig(fname_fig, figsize=(12.5,7)  )


