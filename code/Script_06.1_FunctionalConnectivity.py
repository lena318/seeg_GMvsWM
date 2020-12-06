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

path = "/media/arevell/sharedSSD1/linux/papers/paper005" #Parent directory of project
import pickle
import numpy as np
import os
import sys
from os.path import join as ospj
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
import pingouin as pg
import statsmodels.api as sm
import bct
import networkx as nx
np.seterr(divide = 'ignore')



#%%
#from NetworkX:

def circular_layout(G, scale=1, center=None, dim=2):
    """Position nodes on a circle."""
    import numpy as np
    if dim < 2:
        raise ValueError("cannot handle dimensions < 2")
    paddims = max(0, (dim - 2))
    if len(G) == 0:
        pos = {}
    elif len(G) == 1:
        pos = {nx.utils.arbitrary_element(G): center}
    else:
        theta = np.linspace(1, 2, len(G) + 1)[:-1] * 2 * np.pi
        theta = theta.astype(np.float32)
        pos = np.column_stack(
            [-np.cos(theta), np.sin(theta), np.zeros((len(G), paddims))]
        )
        pos = dict(zip(G, pos))
    return pos




#%% Input/Output Paths and File names
ifname_EEG_times = ospj( path, "data/data_raw/iEEG_times/EEG_times.xlsx")
ifpath_electrode_localization = ospj( path, "data/data_processed/electrode_localization")

ifpath_FC = ospj(path, "data/data_processed/connectivity_matrices/function")
ofpath_FC_analysis = ospj(path, "data/data_processed/FC_analysis/montage/referential/filtered/original")
ofpath_figure = ospj(path, "seeg_GMvsWM/figures/networks")

if not (os.path.isdir(ofpath_FC_analysis)): os.makedirs(ofpath_FC_analysis, exist_ok=True)
if not (os.path.isdir(ofpath_figure)): os.makedirs(ofpath_figure, exist_ok=True)
#%% Load Study Meta Data
data = pd.read_excel(ifname_EEG_times)    

#%% Processing Meta Data: extracting sub-IDs

sub_IDs_unique =  np.unique(data.RID)[np.argsort( np.unique(data.RID, return_index=True)[1])]
#%%

diff_FC_GMWM = [np.zeros(1)] *4

FC_GMGM_all_ii = np.zeros(1)
FC_WMWM_all_ii = np.zeros(1)
FC_GMGM_all_pi = np.zeros(1)
FC_WMWM_all_pi = np.zeros(1)
FC_GMGM_all_ic = np.zeros(1)
FC_WMWM_all_ic = np.zeros(1)
FC_GMGM_all_po = np.zeros(1)
FC_WMWM_all_po = np.zeros(1)

calculate_distances_boolean = True #calculating distances matrices to see if GM-GM distances are different from WM-WM distances. takes a while to compute, so only need to do it once
if (calculate_distances_boolean):
    GM_distances = np.zeros(1)
    WM_distances = np.zeros(1)
descriptors = ["interictal","preictal","ictal","postictal"]

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
    
    #Inputs and OUtputs
    #input filename EEG
    ifpath_FC_sub_ID = ospj(ifpath_FC, sub_RID)
    ifpath_electrode_localization_sub_ID =  ospj(ifpath_electrode_localization, sub_RID)
    ifname_FC_filtered = ospj(ifpath_FC_sub_ID, "sub-{0}_{1}_{2}_{3}_functionalConnectivity.pickle".format(sub_ID, iEEG_filename, start_time_usec, stop_time_usec))

    ifname_electrode_localization = ospj(ifpath_electrode_localization_sub_ID, "sub-{0}_electrode_localization.csv".format(sub_ID))
    ofpath_FC_sub_ID = ospj(ofpath_FC_analysis, sub_RID)
    ofname_FC_sub_ID = ospj(ofpath_FC_sub_ID, "sub-{0}_{1}_{2}_{3}_PSDvsDistance_filtered.pickle".format(sub_ID, iEEG_filename, start_time_usec, stop_time_usec))
    if not (os.path.isdir(ofpath_FC_sub_ID)): os.mkdir(ofpath_FC_sub_ID)
    
    
    #GET DATA
    #get localization and FC files
    with open(ifname_FC_filtered, 'rb') as f: broadband, alphatheta, beta, lowgamma, highgamma, electrode_row_and_column_names, order_of_matrices_in_pickle_file = pickle.load(f)
    electrode_localization = pd.read_csv(ifname_electrode_localization)
    FC = [broadband, alphatheta, beta, lowgamma, highgamma]
    
    #remname to standard 4 character electrode name
    for e in range(len(electrode_row_and_column_names)):
        electrode_name = electrode_row_and_column_names[e]
        if (len(electrode_name) == 3): electrode_name = "{0}{1}{2}".format(electrode_name[0:2], 0, electrode_name[2])
        electrode_row_and_column_names[e] = electrode_name
    #get electrode names in FC and localization files
    electrode_names_FC = electrode_row_and_column_names
    electrode_names_localization = np.array(electrode_localization["electrode_name"])
    
    
    
    #Preprocessing files
    #find electrodes in both files (i.e there are electrodes in localization files but not in FC files, and there are electrodes in FC files but not localized)
    electrode_names_intersect, electrode_names_FC_ind, electrode_names_localization_ind =  np.intersect1d(electrode_names_FC, electrode_names_localization, return_indices = True)

    #Equalizing discrepancies in localization files and FC files
    #removing electrodes in localization files not in FC
    electrode_localization_intersect = electrode_localization.iloc[electrode_names_localization_ind]
    electrode_names_localization_intersect = electrode_names_localization[electrode_names_localization_ind]
    
    #removing electrodes in FC not in localization files
    FC_intersect = copy.deepcopy(FC)
    for f in range(len(FC)):
        FC_intersect[f] = FC_intersect[f][electrode_names_FC_ind[:,None], electrode_names_FC_ind[None,:], :] 
    electrode_names_FC_intersect = electrode_names_FC[electrode_names_FC_ind]
    print(np.equal(electrode_names_localization_intersect, electrode_names_FC_intersect))#check electrode localization order aligns with FC
    
        
    #Only GM and WM tissues considered (not CSF=1, or outside brain = 0)
    labels = np.array(electrode_localization_intersect["region_number"])
    labels_gm_ind = np.where(labels >= 2)[0] #Only GM and WM tissues considered
    #removing electrode localization electrodes not in GM or WM
    electrode_localization_intersect_GMWM = electrode_localization_intersect.iloc[labels_gm_ind]
    #removing FC electrodes not in GM or WM
    FC_intersect_GMWM = copy.deepcopy(FC_intersect)
    for f in range(len(FC_intersect_GMWM)):
        FC_intersect_GMWM[f] = FC_intersect_GMWM[f][labels_gm_ind[:,None], labels_gm_ind[None,:], :] 
    
    
    
    #averaging FC
    FC_intersect_GMWM_mean = [None] * len(FC_intersect_GMWM)
    for f in range(len(FC_intersect_GMWM_mean)):
        FC_intersect_GMWM_mean[f] = np.nanmean(FC_intersect[f], axis=2)    
    np.allclose(FC_intersect_GMWM_mean[f], FC_intersect_GMWM_mean[f].T, rtol=1e-05, atol=1e-08)#check if symmetric
    
    #plot distributions of FC
    f=0
    #plot = sns.displot( FC_intersect_GMWM_mean[f][  np.triu_indices( len(FC_intersect_GMWM_mean[f]), k = 1)   ], kind = "kde" ); plot.axes[0,0].set_ylim(0,10); plot.axes[0,0].set_xlim(0,1); plt.title("Broadband \n{1} \n{0}".format(sub_RID,descriptor ))
    
    #Get distances
    distances = np.array(electrode_localization_intersect_GMWM["distances_label_2"])
    distances_order_ind = np.argsort(distances)
    distances[distances_order_ind]
    electrode_localization_intersect_GMWM_distances_order = electrode_localization_intersect_GMWM.iloc[distances_order_ind]
    
    #Ordering FC nodes based on distances
    FC_intersect_GMWM_mean_distances_order = copy.deepcopy(FC_intersect_GMWM_mean)
    for f in range(len(FC)):
        FC_intersect_GMWM_mean_distances_order[f] = FC_intersect_GMWM_mean_distances_order[f][distances_order_ind[:,None], distances_order_ind[None,:]] 
    

    #plt.figure(figsize=(10, 10), dpi = 100); sns.heatmap(FC_intersect_GMWM_mean_distances_order[f], square=True, vmin = 0, vmax = 1); plt.title("Broadband \n{1} \n{0}".format(sub_RID,descriptor ))



    #Subset GM and WM FC
    distances = np.array(electrode_localization_intersect_GMWM_distances_order["distances_label_2"])
    distances_GM_ind = np.where(distances <= 0)[0]
    distances_WM_ind = np.where(distances > 0)[0]
    #Subset GM
    FC_intersect_GMWM_mean_distances_order_GM = copy.deepcopy(FC_intersect_GMWM_mean_distances_order)
    for f in range(len(FC)):
        FC_intersect_GMWM_mean_distances_order_GM[f] = FC_intersect_GMWM_mean_distances_order_GM[f][distances_GM_ind[:,None], distances_GM_ind[None,:]] 
    #Subset WM
    FC_intersect_GMWM_mean_distances_order_WM = copy.deepcopy(FC_intersect_GMWM_mean_distances_order)
    for f in range(len(FC)):
        FC_intersect_GMWM_mean_distances_order_WM[f] = FC_intersect_GMWM_mean_distances_order_WM[f][distances_WM_ind[:,None], distances_WM_ind[None,:]] 
    

    #plt.figure(figsize=(10, 10), dpi = 100); sns.heatmap(FC_intersect_GMWM_mean_distances_order_GM[f], square=True, vmin = 0, vmax = 1); plt.title("Broadband GM-GM \n{1} \n{0}".format(sub_RID,descriptor ))
    #plt.figure(figsize=(10, 10), dpi = 100); sns.heatmap(FC_intersect_GMWM_mean_distances_order_WM[f], square=True, vmin = 0, vmax = 1); plt.title("Broadband WM-WM \n{1} \n{0}".format(sub_RID,descriptor ))
    f = 0
    FC_GMGM = FC_intersect_GMWM_mean_distances_order_GM[f][np.triu_indices( len(FC_intersect_GMWM_mean_distances_order_GM[f]), k = 1) ] 
    FC_WMWM = FC_intersect_GMWM_mean_distances_order_WM[f][np.triu_indices( len(FC_intersect_GMWM_mean_distances_order_WM[f]), k = 1) ] 


    df_long_GMGM = pd.DataFrame({"Tissue": np.repeat("GM-GM", len(FC_GMGM)), "FC": FC_GMGM})
    df_long_WMWM = pd.DataFrame({"Tissue": np.repeat("WM-WM", len(FC_WMWM)), "FC": FC_WMWM})
    
    FC_df_long = pd.concat([df_long_GMGM, df_long_WMWM] )

    np.nanmean(FC_GMGM )
    np.nanmean(FC_WMWM  )
    np.nanmean(FC_WMWM  ) - np.nanmean(FC_GMGM ) #df_diff_FC_GMWM
    if descriptor == descriptors[0]: diff_FC_GMWM[0] = np.concatenate([diff_FC_GMWM[0] ,np.array([np.nanmean(FC_WMWM  ) - np.nanmean(FC_GMGM )] )  ])
    if descriptor == descriptors[1]: diff_FC_GMWM[1] = np.concatenate([diff_FC_GMWM[1] ,np.array([np.nanmean(FC_WMWM  ) - np.nanmean(FC_GMGM )] )  ])
    if descriptor == descriptors[2]: diff_FC_GMWM[2] = np.concatenate([diff_FC_GMWM[2] ,np.array([np.nanmean(FC_WMWM  ) - np.nanmean(FC_GMGM )] )  ])
    if descriptor == descriptors[3]: diff_FC_GMWM[3] = np.concatenate([diff_FC_GMWM[3] ,np.array([np.nanmean(FC_WMWM  ) - np.nanmean(FC_GMGM )] )  ])

    #plot = sns.displot( data =FC_df_long  , kind = "kde", x="FC" , hue="Tissue"); plot.axes[0,0].set_ylim(0,10); plot.axes[0,0].set_xlim(0,1); plt.title("Broadband \n{1} \n{0}".format(sub_RID,descriptor ))

    
    if descriptor == descriptors[0]:
        FC_GMGM_all_ii = np.concatenate([FC_GMGM_all_ii,FC_GMGM ])
        FC_WMWM_all_ii = np.concatenate([FC_WMWM_all_ii,FC_WMWM ])
        per = 0
    if descriptor == descriptors[1]:
        FC_GMGM_all_pi = np.concatenate([FC_GMGM_all_pi,FC_GMGM ])
        FC_WMWM_all_pi = np.concatenate([FC_WMWM_all_pi,FC_WMWM ])
        per = 1
    if descriptor == descriptors[2]:
        FC_GMGM_all_ic = np.concatenate([FC_GMGM_all_ic,FC_GMGM ])
        FC_WMWM_all_ic = np.concatenate([FC_WMWM_all_ic,FC_WMWM ])
        per = 2
    if descriptor == descriptors[3]:
        FC_GMGM_all_po = np.concatenate([FC_GMGM_all_po,FC_GMGM ])
        FC_WMWM_all_po = np.concatenate([FC_WMWM_all_po,FC_WMWM ])
        per = 3
        
        
    #Network Analysis
        
    #Subset GM
    FC_intersect_GMWM_mean_distances_order
    for f in range(len(FC)):
        FC_intersect_GMWM_mean_distances_order_GM[f] = FC_intersect_GMWM_mean_distances_order_GM[f][distances_GM_ind[:,None], distances_GM_ind[None,:]] 
    #Subset WM
    FC_intersect_GMWM_mean_distances_order_WM = copy.deepcopy(FC_intersect_GMWM_mean_distances_order)
    for f in range(len(FC)):
        FC_intersect_GMWM_mean_distances_order_WM[f] = FC_intersect_GMWM_mean_distances_order_WM[f][distances_WM_ind[:,None], distances_WM_ind[None,:]] 
        
    tmp = copy.deepcopy(FC_intersect_GMWM_mean_distances_order[f])
    tmp[np.where(tmp <0.3 ) ] = 0
    G = nx.from_numpy_array(tmp) 
    G.edges(data=True)
    pos = circular_layout(G)
    #plt.figure(figsize=(10, 10), dpi = 100); nx.draw_circular(G, with_labels=True, font_weight='bold')

    labels = {} 
    for idx, node in enumerate(G.nodes()): 
        labels[node] = np.array(electrode_localization_intersect_GMWM_distances_order["electrode_name"])[idx] 
    colors_node = [] 
    for idx, node in enumerate(G.nodes()): 
        region_id = np.array(electrode_localization_intersect_GMWM_distances_order["region_number"])[idx] 
        if region_id ==3: clr = "#b6d4ee"
        if region_id ==2: clr = "#c6b4a5"
        colors_node.append(clr)
    colors_edge = [] 
    for idx, edge in enumerate(G.edges()): 
        region_id_1 = np.array(electrode_localization_intersect_GMWM_distances_order["region_number"])[edge[0]] 
        region_id_2 = np.array(electrode_localization_intersect_GMWM_distances_order["region_number"])[edge[1]] 
        if (region_id_1 ==3 and region_id_2 == 3 ): clr = "#8dbce4"
        if (region_id_1 ==2 and region_id_2 == 2 ): clr = "#90735c"
        if (region_id_1 ==2 and region_id_2 == 3 ): clr = "#33333333"
        if (region_id_1 ==3 and region_id_2 == 2 ): clr = "#33333333"
        colors_edge.append(clr)
        
    threshold = 0.4

    edges = [(u,v,d) for (u,v,d) in G.edges(data=True) if G.has_edge(v,u)]
    colors_edge = [] 
    for idx, edge in enumerate(G.edges()): 
        if edges[idx][2]["weight"] > threshold:
            region_id_1 = np.array(electrode_localization_intersect_GMWM_distances_order["region_number"])[edge[0]] 
            region_id_2 = np.array(electrode_localization_intersect_GMWM_distances_order["region_number"])[edge[1]] 
            if (region_id_1 ==3 and region_id_2 == 3 ): clr = "#8dbce4"
            if (region_id_1 ==2 and region_id_2 == 2 ): clr = "#90735c"
            if (region_id_1 ==2 and region_id_2 == 3 ): clr = "#33333333"
            if (region_id_1 ==3 and region_id_2 == 2 ): clr = "#33333333"
            colors_edge.append(clr)
        
    plt.figure(figsize=(10, 10), dpi = 300)
    nx.draw_networkx_nodes(G, pos, node_size=100, node_color=colors_node)
    elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > threshold]
    nx.draw_networkx_edges(G, pos, edgelist= elarge, width=1, edge_color = colors_edge)
    nx.draw_networkx_labels(G, pos, font_size=5, font_family="sans-serif", labels = labels)
    plt.title("Broadband \n{1} \n{0}".format(sub_RID,descriptor ))
    
    ofpath_figure_circles = ospj(ofpath_figure, "circles")
    if not (os.path.isdir(ofpath_figure_circles)): os.makedirs(ofpath_figure_circles, exist_ok=True)
    plt.savefig(ospj(ofpath_figure_circles, "Broadband_{0}_{2}_{1}".format(sub_RID,descriptor, per )))
    

    #calculating distances matrices to see if GM-GM distances are different from WM-WM distances
    if (calculate_distances_boolean): #takes a while to compute, so only need to do it once
        if descriptor == descriptors[0]:
            
            electrode_localization_intersect_GMWM_distances_order
            distance_matrix = np.zeros(shape = (len(electrode_localization_intersect_GMWM_distances_order),len(electrode_localization_intersect_GMWM_distances_order) ))
            for e1 in range(len(electrode_localization_intersect_GMWM_distances_order)):
                for e2 in range(len(electrode_localization_intersect_GMWM_distances_order)):
                    p1 = [electrode_localization_intersect_GMWM_distances_order.iloc[e1]["x_coordinate"], electrode_localization_intersect_GMWM_distances_order.iloc[e1]["y_coordinate"], electrode_localization_intersect_GMWM_distances_order.iloc[e1]["z_coordinate"]]
                    p2 = [electrode_localization_intersect_GMWM_distances_order.iloc[e2]["x_coordinate"], electrode_localization_intersect_GMWM_distances_order.iloc[e2]["y_coordinate"], electrode_localization_intersect_GMWM_distances_order.iloc[e2]["z_coordinate"]]
                    distance_matrix[e1,e2] = np.sqrt(  np.sum((np.array(p1)-np.array(p2))**2, axis=0)    )
            #plt.figure(figsize=(10, 10), dpi = 100); sns.heatmap(distance_matrix, square=True); plt.title("Distance Matrix \n{0} ".format(sub_RID ))
        
            #Subset Distances by GM-GM, WM-WM pairs
            distance_matrix_GM = distance_matrix[distances_GM_ind[:,None], distances_GM_ind[None,:]] 
            distance_matrix_WM = distance_matrix[distances_WM_ind[:,None], distances_WM_ind[None,:]] 
            #plt.figure(figsize=(10, 10), dpi = 100); sns.heatmap(distance_matrix_GM, square=True); plt.title("Distance Matrix GM \n{0} ".format(sub_RID ))
            #plt.figure(figsize=(10, 10), dpi = 100); sns.heatmap(distance_matrix_WM, square=True); plt.title("Distance Matrix WM \n{0} ".format(sub_RID ))
            distance_GM = distance_matrix_GM[np.triu_indices( len(distance_matrix_GM), k = 1) ] 
            distance_WM = distance_matrix_WM[np.triu_indices( len(distance_matrix_WM), k = 1) ] 
            
            np.nanmean(distance_GM )
            np.nanmean(distance_WM  )
            GM_distances = np.concatenate([GM_distances, distance_GM]   )
            WM_distances = np.concatenate([WM_distances, distance_WM]    )
            
            mannwhitneyu(distance_GM ,distance_WM  )[1]

for l in range(len(diff_FC_GMWM)):
    diff_FC_GMWM[l] = np.delete(diff_FC_GMWM[l], 0)
FC_GMGM_all_ii = np.delete(FC_GMGM_all_ii, 0); FC_WMWM_all_ii = np.delete(FC_WMWM_all_ii, 0) #removing the initialized zero
FC_GMGM_all_pi = np.delete(FC_GMGM_all_pi, 0); FC_WMWM_all_pi = np.delete(FC_WMWM_all_pi, 0) #removing the initialized zero
FC_GMGM_all_ic = np.delete(FC_GMGM_all_ic, 0); FC_WMWM_all_ic = np.delete(FC_WMWM_all_ic, 0) #removing the initialized zero
FC_GMGM_all_po = np.delete(FC_GMGM_all_po, 0); FC_WMWM_all_po = np.delete(FC_WMWM_all_po, 0) #removing the initialized zero
#%%


df_long_GMGM_all_ii = pd.DataFrame({"Tissue": np.repeat("GM-GM", len(FC_GMGM_all_ii)), "FC": FC_GMGM_all_ii})
df_long_GMGM_all_pi = pd.DataFrame({"Tissue": np.repeat("GM-GM", len(FC_GMGM_all_pi)), "FC": FC_GMGM_all_pi})
df_long_GMGM_all_ic = pd.DataFrame({"Tissue": np.repeat("GM-GM", len(FC_GMGM_all_ic)), "FC": FC_GMGM_all_ic})
df_long_GMGM_all_po = pd.DataFrame({"Tissue": np.repeat("GM-GM", len(FC_GMGM_all_po)), "FC": FC_GMGM_all_po})

df_long_WMWM_all_ii = pd.DataFrame({"Tissue": np.repeat("WM-WM", len(FC_WMWM_all_ii)), "FC": FC_WMWM_all_ii})
df_long_WMWM_all_pi = pd.DataFrame({"Tissue": np.repeat("WM-WM", len(FC_WMWM_all_pi)), "FC": FC_WMWM_all_pi})
df_long_WMWM_all_ic = pd.DataFrame({"Tissue": np.repeat("WM-WM", len(FC_WMWM_all_ic)), "FC": FC_WMWM_all_ic})
df_long_WMWM_all_po = pd.DataFrame({"Tissue": np.repeat("WM-WM", len(FC_WMWM_all_po)), "FC": FC_WMWM_all_po})


FC_df_long_all_ii = pd.concat([df_long_GMGM_all_ii, df_long_WMWM_all_ii] )
FC_df_long_all_pi = pd.concat([df_long_GMGM_all_pi, df_long_WMWM_all_pi] )
FC_df_long_all_ic = pd.concat([df_long_GMGM_all_ic, df_long_WMWM_all_ic] )
FC_df_long_all_po = pd.concat([df_long_GMGM_all_po, df_long_WMWM_all_po] )
plot = sns.displot( data =FC_df_long_all_ii  , kind = "ecdf", x="FC" , hue="Tissue"); plot.axes[0,0].set_ylim(0,1); plot.axes[0,0].set_xlim(0,1); plt.title("{0} \nInterictal".format(order_of_matrices_in_pickle_file.iloc[f][0] ))
plot = sns.displot( data =FC_df_long_all_pi  , kind = "ecdf", x="FC" , hue="Tissue"); plot.axes[0,0].set_ylim(0,1); plot.axes[0,0].set_xlim(0,1); plt.title("{0} \nPreictal".format(order_of_matrices_in_pickle_file.iloc[f][0]  ))
plot = sns.displot( data =FC_df_long_all_ic  , kind = "ecdf", x="FC" , hue="Tissue"); plot.axes[0,0].set_ylim(0,1); plot.axes[0,0].set_xlim(0,1); plt.title("{0} \nIctal".format(order_of_matrices_in_pickle_file.iloc[f][0] ))
plot = sns.displot( data =FC_df_long_all_po  , kind = "ecdf", x="FC" , hue="Tissue"); plot.axes[0,0].set_ylim(0,1); plot.axes[0,0].set_xlim(0,1); plt.title("{0} \nPostictal".format(order_of_matrices_in_pickle_file.iloc[f][0]  ))



np.nanmean(FC_GMGM_all_ii )
np.nanmean(FC_WMWM_all_ii  )
print(np.nanmean(FC_WMWM_all_ii ) - np.nanmean(FC_GMGM_all_ii  ))

np.nanmean(FC_GMGM_all_pi )
np.nanmean(FC_WMWM_all_pi  )
print(np.nanmean(FC_WMWM_all_pi ) - np.nanmean(FC_GMGM_all_pi  ))

np.nanmean(FC_GMGM_all_ic )
np.nanmean(FC_WMWM_all_ic  )
print(np.nanmean(FC_WMWM_all_ic ) - np.nanmean(FC_GMGM_all_ic  ))

np.nanmean(FC_GMGM_all_po )
np.nanmean(FC_WMWM_all_po  )
print(np.nanmean(FC_WMWM_all_po ) - np.nanmean(FC_GMGM_all_po  ))


print(mannwhitneyu(FC_GMGM_all_ii, FC_WMWM_all_ii)[1])
print(mannwhitneyu(FC_GMGM_all_pi, FC_WMWM_all_pi)[1])
print(mannwhitneyu(FC_GMGM_all_ic, FC_WMWM_all_ic)[1])
print(mannwhitneyu(FC_GMGM_all_po, FC_WMWM_all_po)[1])

#%%

#calculating distances matrices to see if GM-GM distances are different from WM-WM distances
if (calculate_distances_boolean): 
    GM_distances = np.delete(GM_distances, 0); WM_distances = np.delete(WM_distances, 0) #removing the initialized zero
    df_long_distances_GM = pd.DataFrame({"Tissue": np.repeat("GM-GM", len(GM_distances)), "distances": GM_distances})
    df_long_distances_WM = pd.DataFrame({"Tissue": np.repeat("WM-WM", len(WM_distances)), "distances": WM_distances})
    df_long_distances = pd.concat([df_long_distances_GM, df_long_distances_WM] )
    np.nanmean(GM_distances )
    np.nanmean(WM_distances  )
    print(np.nanmean(GM_distances ) - np.nanmean(WM_distances  ))
    plot = sns.displot( data =df_long_distances  , kind = "hist", x="distances" , hue="Tissue", kde=True); plot.axes[0,0].set_ylim(0,); plot.axes[0,0].set_xlim(0,); plt.title("{0} \nPostictal".format(order_of_matrices_in_pickle_file.iloc[f][0]  ))
    
    print(mannwhitneyu(GM_distances ,WM_distances  )[1])






#%%
#Regressing out distances
FC_df_long_all = [FC_df_long_all_ii, FC_df_long_all_pi, FC_df_long_all_ic, FC_df_long_all_po]
per = range(4)
for p in range(len(per)):
    FC_vs_distance = pd.concat([FC_df_long_all[p], df_long_distances.drop(["Tissue"], axis=1)   ], axis=1)
    plt.figure(dpi = 100)
    plot =  sns.regplot(data =FC_vs_distance, x = "distances", y = "FC", marker='o', scatter_kws={'s':0.5});  plt.title("Distances vs FC {0}\n{1}".format(order_of_matrices_in_pickle_file.iloc[f][0], descriptors[p] ))
    
    
    
    linear_regression = pg.linear_regression(FC_vs_distance['distances'], FC_vs_distance['FC'])
    
    
    
    all(FC_df_long_all_ii["Tissue"] == df_long_distances["Tissue"])
    
    
    OLS_model = sm.OLS(FC_vs_distance['FC'],FC_vs_distance['distances']).fit()  # training the model
    predicted_values = OLS_model.predict()  # predicted values
    residual_values = OLS_model.resid # residual values
    residual_values = residual_values.rename("residuals")
    
    FC_vs_distance_residuals = pd.concat([FC_vs_distance, residual_values],  axis=1)
    
    
    plot = sns.displot( data =FC_vs_distance_residuals  , kind = "ecdf", x="residuals" , hue="Tissue"); plot.axes[0,0].set_ylim(0,); plt.title("{0} \n{1}".format(order_of_matrices_in_pickle_file.iloc[f][0], descriptors[p]  ))
    FC_vs_distance_residuals["Tissue"]
    

    print("Significance: {0}, {1}: {2}".format(order_of_matrices_in_pickle_file.iloc[f][0], descriptors[p] ,    mannwhitneyu(FC_vs_distance_residuals.loc[FC_vs_distance_residuals['Tissue'] == "GM-GM"]["residuals"] ,FC_vs_distance_residuals.loc[FC_vs_distance_residuals['Tissue'] == "WM-WM"]["residuals"]  )[1]))
  
    
#%%
# Set your custom color palette
colors = ["#c6b4a5", "#b6d4ee"]
sns.set_palette(sns.color_palette(colors))

diff_FC_GMWM



{"sub_ID": sub_IDs_unique, "interictal":diff_FC_GMWM[0] , "preictal":diff_FC_GMWM[1], "ictal":diff_FC_GMWM[2], "postictal":diff_FC_GMWM[3] }
df_diff_FC_GMWM = pd.DataFrame({"sub_ID": sub_IDs_unique, "interictal":diff_FC_GMWM[0] , "preictal":diff_FC_GMWM[1], "ictal":diff_FC_GMWM[2], "postictal":diff_FC_GMWM[3] })
#df_diff_FC_GMWM = df_diff_FC_GMWM.drop([4,6,7])
df_diff_FC_GMWM_long = pd.melt(df_diff_FC_GMWM, id_vars=['sub_ID'], value_vars=['interictal', 'preictal', 'ictal', 'postictal'])




plt.figure(dpi = 300); 
flierprops = dict(marker='o', markerfacecolor='#00000088', markersize=0.0, linestyle='none', markeredgecolor='#00000000')
ax = sns.boxplot(x="variable", y="value", data=df_diff_FC_GMWM_long,  linewidth=0.9, whis = 1, width=0.8, flierprops = flierprops, showfliers = False)
ax = sns.swarmplot(x="variable", y="value", data=df_diff_FC_GMWM_long, color=".25")


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.set_ylim([-0.015,0.2])
print(df_diff_FC_GMWM)
#%%

from scipy import stats

# estimate sample size via power analysis
from statsmodels.stats.power import TTestIndPower


def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)
effect = cohen_d(df_diff_FC_GMWM[ "ictal"], df_diff_FC_GMWM["preictal"])
alpha = 0.05/20
power = 0.8
# perform power analysis
analysis = TTestIndPower()
result = analysis.solve_power(effect, power=power, nobs1=None, ratio=1.0, alpha=alpha)
print('Sample Size: %.3f' % result)
effect = cohen_d(df_diff_FC_GMWM[ "ictal"], df_diff_FC_GMWM["interictal"])
result = analysis.solve_power(effect, power=power, nobs1=None, ratio=1.0, alpha=alpha)
print('Sample Size: %.3f' % result)
print(stats.ttest_rel(diff_FC_GMWM[1],diff_FC_GMWM[2])[1])



np.mean(diff_FC_GMWM[0])
np.mean(diff_FC_GMWM[1])
np.mean(diff_FC_GMWM[2])
np.mean(diff_FC_GMWM[3])


print(mannwhitneyu(diff_FC_GMWM[1] ,diff_FC_GMWM[2]  )[1])

df_diff_FC_GMWM.mean(axis=1)


plot = sns.displot( data =df_diff_FC_GMWM_long  , kind = "kde", x="value" , hue="sub_ID"); plot.axes[0,0].set_ylim(0,); plt.title("{0} \n{1}".format(order_of_matrices_in_pickle_file.iloc[f][0], descriptors[p]  ))
plot = sns.displot( data =df_diff_FC_GMWM.mean(axis=1)  , kind = "kde"); plot.axes[0,0].set_ylim(0,); plt.title("{0} \n{1}".format(order_of_matrices_in_pickle_file.iloc[f][0], descriptors[p]  ))

stats.ttest_1samp(diff_FC_GMWM[2]  , 0)[1]
stats.ttest_1samp(  np.concatenate([diff_FC_GMWM[0], diff_FC_GMWM[1],diff_FC_GMWM[2],diff_FC_GMWM[3] ])  , 0)[1]
stats.ttest_1samp( np.array(df_diff_FC_GMWM.mean(axis=1))  , 0)[1]
effect = np.mean(np.array(df_diff_FC_GMWM.mean(axis=1))) /np.std(np.array(df_diff_FC_GMWM.mean(axis=1)) )
result = analysis.solve_power(effect, power=power, nobs1=None, ratio=1.0, alpha=alpha)
import statsmodels
statsmodels.stats.power.tt_solve_power(effect_size=effect, nobs=None, alpha=0.05, power=0.8, alternative='two-sided')


FC_df_long_all_ii 
FC_df_long_all_pi 
FC_df_long_all_ic 
FC_df_long_all_po 

#%%

df_diff_FC_GMWM_by_elec = pd.DataFrame({"Tissue": FC_df_long_all_ii["Tissue"], "interictal":FC_df_long_all_ii["FC"] , "preictal":FC_df_long_all_pi["FC"] , "ictal":FC_df_long_all_ic["FC"] , "postictal":FC_df_long_all_po["FC"]  })

df_diff_FC_GMWM_by_elec_long = pd.melt(df_diff_FC_GMWM_by_elec, id_vars=['Tissue'], value_vars=['interictal', 'preictal', 'ictal', 'postictal'])




sns.violinplot(x="variable", y="value", hue="Tissue", data=df_diff_FC_GMWM_by_elec_long, palette="muted")


# Set your custom color palette
colors = ["#c6b4a5", "#b6d4ee"]
sns.set_palette(sns.color_palette(colors))

plt.figure(dpi = 300); 
flierprops = dict(marker='o', markerfacecolor='#00000000', markersize=0.0, linestyle='none', markeredgecolor='#00000033')
ax = sns.boxplot(x="variable", y="value", hue="Tissue", data=df_diff_FC_GMWM_by_elec_long,  linewidth=0.9, whis = 1, width=0.8, flierprops = flierprops)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([0.1,0.65])



#%%
#from NetworkX:

def circular_layout(G, scale=1, center=None, dim=2):
    # dim=2 only
    """Position nodes on a circle.

    Parameters
    ----------
    G : NetworkX graph or list of nodes
        A position will be assigned to every node in G.

    scale : number (default: 1)
        Scale factor for positions.

    center : array-like or None
        Coordinate pair around which to center the layout.

    dim : int
        Dimension of layout.
        If dim>2, the remaining dimensions are set to zero
        in the returned positions.
        If dim<2, a ValueError is raised.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node

    Raises
    -------
    ValueError
        If dim < 2

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> pos = nx.circular_layout(G)

    Notes
    -----
    This algorithm currently only works in two dimensions and does not
    try to minimize edge crossings.

    """
    import numpy as np

    if dim < 2:
        raise ValueError("cannot handle dimensions < 2")


    paddims = max(0, (dim - 2))

    if len(G) == 0:
        pos = {}
    elif len(G) == 1:
        pos = {nx.utils.arbitrary_element(G): center}
    else:
        # Discard the extra angle since it matches 0 radians.
        theta = np.linspace(1, 2, len(G) + 1)[:-1] * 2 * np.pi
        theta = theta.astype(np.float32)
        pos = np.column_stack(
            [-np.cos(theta), np.sin(theta), np.zeros((len(G), paddims))]
        )
        pos = dict(zip(G, pos))

    return pos

circular_layout(G)
