"""
2020.08.01
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Logic of code:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
path = "/media/arevell/sharedSSD/linux/papers/paper005" 
import sys
import os
from os.path import join as ospj
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#import custom
sys.path.append(ospj(path, "seeg_GMvsWM", "code", "tools"))
sys.path.append(ospj(path, "seeg_GMvsWM", "code" ,"tools", "ieegpy"))
from preprocessEEG import filter_eeg_data
from imagingToolsRevell import printProgressBar, show_eeg, plot_adj, plot_adj_allbands
import echobase

#%

ifname_EEG_times = ospj(path,"data/data_raw/iEEG_times/EEG_times.xlsx")
ifpath_EEG = ospj(path,"data/data_raw/EEG")
ifpath_electrodeLoc = ospj(path,"data/data_processed/electrode_localization")
ofpath_electrodesRemoved = ospj(path,"data/data_processed/EEG/electrodesRemoved/")
ofpath_filtered_eeg = ospj(path,"data/data_processed/EEG/filtered/")

if not (os.path.isdir(ofpath_filtered_eeg)): os.makedirs(ofpath_filtered_eeg, exist_ok=True)
if not (os.path.isdir(ofpath_electrodesRemoved)): os.makedirs(ofpath_electrodesRemoved, exist_ok=True)

#% Load Study Meta Data
eegTimes = pd.read_excel(ifname_EEG_times)    


#%%Remove electrodes outside brain

for i in range(len(eegTimes)):
    #parsing data DataFrame to get iEEG information
    sub_ID = eegTimes.iloc[i].RID
    print(sub_ID)
    iEEG_filename = eegTimes.iloc[i].file
    ignore_electrodes = eegTimes.iloc[i].ignore_electrodes.split(",")
    start_time_usec = int(eegTimes.iloc[i].connectivity_start_time_seconds*1e6)
    stop_time_usec = int(eegTimes.iloc[i].connectivity_end_time_seconds*1e6)
    descriptor = eegTimes.iloc[i].descriptor
    #input filename EEG
    ifpath_EEG_subID = ospj(ifpath_EEG, f"sub-{sub_ID}")
    ifname_EEG_subID = ospj(ifpath_EEG_subID, f"sub-{sub_ID}_{iEEG_filename}_{start_time_usec}_{stop_time_usec}_EEG.csv")
    ifname_EEG_subID_MD = ospj(ifpath_EEG_subID, f"sub-{sub_ID}_{iEEG_filename}_{start_time_usec}_{stop_time_usec}_EEG_metadata.csv")
    #input electrode localization
    ifname_electrodeLoc = ospj(ifpath_electrodeLoc, f"sub-{sub_ID}", f"sub-{sub_ID}_electrode_localization.csv")
    
    #read files
    eeg = pd.read_csv(ifname_EEG_subID)
    elecLoc = pd.read_csv(ifname_electrodeLoc)
    eeg_metadata = pd.read_csv(ifname_EEG_subID_MD)
    
    
    #remove eeg electrodes not in localization file
    name_loc = np.array(elecLoc["electrode_name"])
    name_eeg = np.array(eeg.columns) 
    
    name_intersect = np.intersect1d(name_loc, name_eeg, return_indices = True )
    eeg_ordered = eeg.iloc[:, name_intersect[2] ]
    name_eegOrdered = np.array(eeg_ordered.columns) 

    #remove eeg electrodes outside the brain
    ind_outside = np.where(np.array(elecLoc["Tissue_segmentation_region_number"]) == 0)
    name_outside = name_loc[ind_outside]
    
    name_outside_intersect = np.intersect1d(name_outside, name_eegOrdered, return_indices = True )
    if len(name_outside_intersect[2]) > 0:
        eeg_ordered_dropOutside = eeg_ordered.drop(eeg_ordered.columns[name_outside_intersect[2]], axis=1)
    else:
        eeg_ordered_dropOutside = eeg_ordered
    #Output electrodes removed EEG
    ofpath_EEG_subID = ospj(ofpath_electrodesRemoved, f"sub-{sub_ID}")
    ofname_EEG_subID = ospj(ofpath_EEG_subID, f"sub-{sub_ID}_{iEEG_filename}_{start_time_usec}_{stop_time_usec}_EEGremoved.csv")
    ofname_EEG_subID_MD = ospj(ofpath_EEG_subID, f"sub-{sub_ID}_{iEEG_filename}_{start_time_usec}_{stop_time_usec}_EEGremoved_metadata.csv")
    if not (os.path.isdir(ofpath_EEG_subID)): os.mkdir(ofpath_EEG_subID)#if the path doesn't exists, then make the directory
    pd.DataFrame.to_csv(eeg_ordered_dropOutside, ofname_EEG_subID, header=True, index=False)
    pd.DataFrame.to_csv(eeg_metadata, ofname_EEG_subID_MD, header=True, index=False)

#%%
for i in range(len(eegTimes)):
    #parsing data DataFrame to get iEEG information
    sub_ID = eegTimes.iloc[i].RID
    print(sub_ID)
    iEEG_filename = eegTimes.iloc[i].file
    ignore_electrodes = eegTimes.iloc[i].ignore_electrodes.split(",")
    start_time_usec = int(eegTimes.iloc[i].connectivity_start_time_seconds*1e6)
    stop_time_usec = int(eegTimes.iloc[i].connectivity_end_time_seconds*1e6)
    descriptor = eegTimes.iloc[i].descriptor
    #input filename EEG
    ifpath_EEG_subID = ospj(ofpath_electrodesRemoved, f"sub-{sub_ID}")
    ifname_EEG_subID = ospj(ifpath_EEG_subID, f"sub-{sub_ID}_{iEEG_filename}_{start_time_usec}_{stop_time_usec}_EEGremoved.csv")
    ifname_EEG_subID_MD = ospj(ifpath_EEG_subID, f"sub-{sub_ID}_{iEEG_filename}_{start_time_usec}_{stop_time_usec}_EEGremoved_metadata.csv")

    #read files
    eeg = pd.read_csv(ifname_EEG_subID)
    eeg_metadata = pd.read_csv(ifname_EEG_subID_MD)
    fs = int(eeg_metadata["fs"][0])
    data = np.array(eeg)
    
    #calculate parameters
    n_samp, n_chan = data.shape
    mw = 1 #Moving Window; in seconds
    ws = 2 #Window Size; in seconds
    times = np.floor(np.arange(0, n_samp, fs * mw))
    times = times[np.flatnonzero( n_samp - times > fs * ws)].astype(int) #only consider times where can fit into window size
    t = len(times)

    
 
    adj_xcorr = np.zeros(shape = (8, n_chan, n_chan , t))
    adj_pearson = np.zeros(shape = (8, n_chan, n_chan , t))
    adj_pearson_pval = np.zeros(shape = (8, n_chan, n_chan , t))
    adj_spearman = np.zeros(shape = (8, n_chan, n_chan , t))
    adj_spearman_pval = np.zeros(shape = (8, n_chan, n_chan , t))
    adj_coherence = np.zeros(shape = (8, n_chan, n_chan , t))
    adj_mi = np.zeros(shape = (8, n_chan, n_chan , t))

    

    count = 0
    for w in times:
        
        start = w
        stop = w + fs * ws

        data_win = data[start:stop,:]
        
        #show_eeg(data_win, fs, channel = 2)
        
        adj_xcorr[0, :,:,count], adj_xcorr[1, :,:,count], adj_xcorr[2, :,:,count], adj_xcorr[3, :,:,count], adj_xcorr[4, :,:,count], adj_xcorr[5, :,:,count], adj_xcorr[6, :,:,count], adj_xcorr[7, :,:,count] = echobase.crossCorrelation_wrapper(data_win, fs, avgref=True)
        #[adj_pearson[0, :,:,count], adj_pearson[1, :,:,count], adj_pearson[2, :,:,count], adj_pearson[3, :,:,count], adj_pearson[4, :,:,count], adj_pearson[5, :,:,count], adj_pearson[6, :,:,count], adj_pearson[7, :,:,count]], [adj_pearson_pval[0, :,:,count], adj_pearson_pval[1, :,:,count], adj_pearson_pval[2, :,:,count], adj_pearson_pval[3, :,:,count], adj_pearson_pval[4, :,:,count], adj_pearson_pval[5, :,:,count], adj_pearson_pval[6, :,:,count], adj_pearson_pval[7, :,:,count]]  = echobase.pearson_wrapper(data_win, fs, avgref=True)
        #[adj_spearman[0, :,:,count], adj_spearman[1, :,:,count], adj_spearman[2, :,:,count], adj_spearman[3, :,:,count], adj_spearman[4, :,:,count], adj_spearman[5, :,:,count], adj_spearman[6, :,:,count], adj_spearman[7, :,:,count]], [adj_spearman_pval[0, :,:,count], adj_spearman_pval[1, :,:,count], adj_spearman_pval[2, :,:,count], adj_spearman_pval[3, :,:,count], adj_spearman_pval[4, :,:,count], adj_spearman_pval[5, :,:,count], adj_spearman_pval[6, :,:,count], adj_spearman_pval[7, :,:,count]]  = echobase.spearman_wrapper(data_win, fs, avgref=True)
        #adj_coherence[0, :,:,count], adj_coherence[1, :,:,count], adj_coherence[2, :,:,count], adj_coherence[3, :,:,count], adj_coherence[4, :,:,count], adj_coherence[5, :,:,count], adj_coherence[6, :,:,count], adj_coherence[7, :,:,count] = echobase.coherence_wrapper(data_win, fs, avgref=True)
        #adj_mi[0, :,:,count], adj_mi[1, :,:,count], adj_mi[2, :,:,count], adj_mi[3, :,:,count], adj_mi[4, :,:,count], adj_mi[5, :,:,count], adj_mi[6, :,:,count], adj_mi[7, :,:,count] = echobase.mutualInformation_wrapper(data_win, fs, avgref=True)
    
        #plot_adj_allbands(adj_coherence[:, :,:,count], vmin= 0 , vmax = 1)
        #plt.show()
        printProgressBar(w, times[-1], prefix = "Progress:", suffix = f"{w}/{n_samp}. {np.round(t,2)} sec. Mov Win: {mw} sec." )
        count = count + 1

    def animation_adj( count=int):
        c=0
        sns.heatmap(adj_xcorr[c, :,:,count], square=True, ax = axes, vmin = 0, vmax = 1, cbar=False)
        axes.set_title(titles[c]+ f" {count}", size=10)
        c = c+1
    fig,axes = plt.subplots(1,1,figsize=(4,4), dpi = 100)
    titles = ["Broadband", "Delta", "Theta", "Alpha", "Beta", "Gamma - Low", "Gamma - Mid", "Gamma - High"] 
                
    animator = animation.FuncAnimation(fig, animation_adj, interval = 200, repeat=True)
    plt.show()
    
    
    
    
    

    
def data_gen(t=0):
    cnt = 0
    while cnt < 1000:
        cnt += 1
        t += 0.1
        yield t, np.sin(2*np.pi*t) * np.exp(-t/10.)
def init():
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlim(0, 10)
    del xdata[:]
    del ydata[:]
    line.set_data(xdata, ydata)
    return line,

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.grid()
xdata, ydata = [], []


def run(data):
    # update the data
    t, y = data
    xdata.append(t)
    ydata.append(y)
    xmin, xmax = ax.get_xlim()

    if t >= xmax:
        ax.set_xlim(xmin, 2*xmax)
        ax.figure.canvas.draw()
    line.set_data(xdata, ydata)

    return line,


ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10,
                          repeat=False, init_func=init)

%matplotlib qt
plt.show()   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#%%
for i in range(len(eegTimes)):
    #parsing data DataFrame to get iEEG information
    sub_ID = eegTimes.iloc[i].RID
    iEEG_filename = eegTimes.iloc[i].file
    ignore_electrodes = eegTimes.iloc[i].ignore_electrodes.split(",")
    start_time_usec = int(eegTimes.iloc[i].connectivity_start_time_seconds*1e6)
    stop_time_usec = int(eegTimes.iloc[i].connectivity_end_time_seconds*1e6)
    descriptor = eegTimes.iloc[i].descriptor
    #input filename EEG
    ifpath_EEG_sub_ID = ospj(ifpath_EEG, "sub-{0}".format(sub_ID))
    #Output filtered EEG
    ofpath_EEG_sub_ID = ospj(ofpath_filtered_eeg, "sub-{0}".format(sub_ID))
    if not (os.path.isdir(ofpath_EEG_sub_ID)): os.mkdir(ofpath_EEG_sub_ID)#if the path doesn't exists, then make the directory
    ifname_EEG = "{0}/sub-{1}_{2}_{3}_{4}_EEG.pickle".format(ifpath_EEG_sub_ID, sub_ID, iEEG_filename, start_time_usec, stop_time_usec)
    ofname_EEG_filtered = "{0}/sub-{1}_{2}_{3}_{4}_EEG_filtered.pickle.pickle".format(ofpath_EEG_sub_ID, sub_ID, iEEG_filename, start_time_usec, stop_time_usec)
    print("\n\n\nID: {0}\nDescriptor: {1}".format(sub_ID, descriptor))
    if (os.path.exists(ofname_EEG_filtered)):
        print("File already exists: {0}".format(ofname_EEG_filtered))
    else:#if file already exists, don't run below
        filter_eeg_data(ifname_EEG, ofname_EEG_filtered)


    

