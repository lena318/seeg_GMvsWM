"""
2020.08.01
Andy Revell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Purpose:
    Plot spectrogram

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
from scipy import signal
import pickle
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

paper_path = "."
#assumes current working directory is paper005/paper005/pipelines/scripts/
input_file_path = os.path.join(paper_path, "data_processed/PSDvsDistance/montage/referential/filtered/averaged")
inputfile_name = "PSDvsDistance_All_patients.pickle"
inputfile = os.path.join(input_file_path, inputfile_name)
with open(inputfile, 'rb') as f: interictal, preictal, ictal ,postictal = pickle.load(f)


output_path = os.path.join(paper_path, "paper005/figures")
outputfile = os.path.join(output_path, "PSDvsDistance.png")



tmp = np.nanmean(preictal, axis=2)
plot = np.delete(np.log10(np.array(tmp)), range(30, np.shape(tmp)[0]), axis=0)
plot = np.delete(plot, range(800, np.shape(tmp)[1]), axis=1)
sns.heatmap(plot, cmap=sns.color_palette("Spectral_r", n_colors=20, desat=0.6), vmin=2,vmax=4.5)
plt.show()

All_data = [interictal, preictal, ictal, postictal]

ylim = [51 800]
fontsize = 30
plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['font.family'] = "serif"
sns.set_context("talk")

gridspec = dict(wspace=0.00, width_ratios=[0.5,0.1, 1, 2, 1], hspace= 0.1)#create a blank space between interictal and preictal
fig, axs = plt.subplots(1, 5, sharex='col',  gridspec_kw=gridspec)
plt.xlabel("Distance (mm)")
plt.ylabel("Frequency (Hz)")
cbar_ax = fig.add_axes([.91, .3, .03, .4])
#(ax1, ax2, ax3, ax4, ax5), (ax6, ax7, ax8, ax9, ax10) = axs
for i in range(0,2):
    k = 0
    for j in range(0, 5):
        if i == 0:
            plot = ALL_data_GM[:, :, k]
        if i == 1:
            plot = ALL_data_WM[:, :, k]
        plot = np.delete(plot, range(ylim, np.shape(plot)[0]), axis=0)
        yticklabels = True
        if j > 1:
            yticklabels=True
        else:
            yticklabels = True
        if j == 1:#create a blank space between interictal and preictal
            axs[i][j].remove()
        if j != 1:
            sns.heatmap(plot,cmap=sns.color_palette("Spectral_r", n_colors=20, desat=0.6), vmin=-1,vmax=2.5, ax = axs[i][j],
                        cbar=i == 0, cbar_ax=None if i else cbar_ax, cbar_kws={'label': '$\mathregular{log_{10}}({V^2}/{Hz})$'})
            old_ticks = axs[i][j].get_xticks()
            new_ticks = np.linspace(np.min(old_ticks), np.max(old_ticks), 5)
            new_ticks = new_ticks[range(1,4)]
            axs[i][j].set_xticks(new_ticks)
            axs[i][j].set_xticklabels(["25%", "50%", "75%"], rotation=0)
            if j == 0:
                axs[i][j].set_xticks([new_ticks[1]])
                axs[i][j].set_xticklabels(["~6hrs before"], rotation=0)
            if j > 1:
                axs[i][j].set_yticks([])
                if j <4:
                    axs[i][j].axvline(x=np.shape(ALL_data_GM)[1], c = "black", linewidth=4, linestyle='--')
            k = k + 1
        axs[i][j].invert_yaxis()
fig.text(0.56, 0.04, 'Time (% of normalized seizure length)', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig.text(0.07, 0.5, 'Frequency (Hz)', va='center', rotation='vertical', fontdict = {'fontsize': fontsize*0.7})
fig.text(0.05, 0.7, 'Gray\nMatter', ha='center', fontdict = {'fontsize': fontsize*1.1})
fig.text(0.05, 0.225, 'White\nMatter', ha='center',  fontdict = {'fontsize': fontsize*1.1})
fig.text(0.165, 0.9, 'Interictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig.text(0.3, 0.9, 'Preictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig.text(0.56, 0.9, 'Ictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig.text(0.81, 0.9, 'Postictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig.text(0.5, 0.95, 'GM vs WM Spectrogram', ha='center', fontdict = {'fontsize': fontsize*1.1})
fig.text(0.9, 0.915, 'Patients: 5\nSeizures: 5\nGM electrodes: 346\nWM electrodes: 185', ha='left', fontdict = {'fontsize': 12})
#plt.show()
plt.savefig(outputfile, dpi = 300)

