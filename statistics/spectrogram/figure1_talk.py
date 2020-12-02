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


~~~~~~~
"""

from scipy import signal
import pickle
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd
from matplotlib.patches import Ellipse

paper_path = "."
#assumes current working directory is paper005/paper005/pipelines/scripts/

#Spectrogram
input_path = os.path.join(paper_path, "data_processed/spectrogram/montage/referential/filtered/interpolated_all_patients_combined")
inputfile = os.path.join(input_path, [f for f in sorted(os.listdir(input_path))][0])
with open(inputfile, 'rb') as f: ALL_data_GM, ALL_data_WM, freqs, xnew = pickle.load(f)

#PSD VS Distance
input_path = os.path.join(paper_path, "data_processed/PSDvsDistance/montage/referential/filtered/averaged")
inputfile_name = "PSDvsDistance_All_patients.pickle"
inputfile = os.path.join(input_path, inputfile_name)
with open(inputfile, 'rb') as f: interictal, preictal, ictal ,postictal = pickle.load(f)


#SNR
input_path = os.path.join(paper_path, "data_processed/SNR/montage/referential/filtered/average")
inputfile_name = "SNR_all_patients_mean.pickle"
inputfile = os.path.join(input_path, inputfile_name)
with open(inputfile, 'rb') as f: SNR_all_patients_mean, SNR_all_patients, bins_new, RID, breaks = pickle.load(f)

output_path = os.path.join(paper_path, "paper005/figures/figure2")
outputfile_spectrogram = os.path.join(output_path, "figure2_talk_spectrogram.png")
outputfile_PSDvsDistance = os.path.join(output_path, "figure2_talk_PSDvsDistance.png")
outputfile_SNR = os.path.join(output_path, "figure2_talk_SNR.png")




#Spectrogram
plt.rcParams['figure.figsize'] = (20.0, 10)
plt.rcParams['font.family'] = "serif"
sns.set_context("talk")
fontsize = 30



#making plot layout
fig1 = plt.figure(constrained_layout=False)
left = 0.07; right = 0.8
gs1 = fig1.add_gridspec(nrows=1, ncols=5, left=left, right=right, bottom = 0.525, top = 0.875 ,wspace=0.00, width_ratios=[1,0.1, 1, 1.2, 1])
gs2 = fig1.add_gridspec(nrows=1, ncols=5, left=left, right=right, bottom = 0.125, top = 0.475 ,wspace=0.00, width_ratios=[1,0.1, 1, 1.2, 1])

gs = [gs1, gs2]
axes = [None] * len(gs)
for a in range(len(gs)):
    axes[a] = [None] * gs[a]._ncols

for a in range(len(axes)):
    for b in range(len(axes[a])):
        axes[a][b] = fig1.add_subplot(gs[a][:, b])

cbar_ax1 = fig1.add_axes([0.865, 0.525, 0.02, 0.875-0.525])
cbar_ax2 = fig1.add_axes([0.865, 0.125, 0.02, 0.875-0.525])
ylim = 31
for i in range(0,2):
    k = 0
    for j in range(0, 5):
        if i == 0:
            plot = ALL_data_GM[:, :, k]
        if i == 1:
            plot = ALL_data_WM[:, :, k]
        plot = np.delete(plot, range(ylim, np.shape(plot)[0]), axis=0)
        plot = np.delete(plot, range(0, 1), axis=0)
        if j == 1:#create a blank space between interictal and preictal
            axes[i][j].remove()
        if j != 1:
            if i == 0:
                sns.heatmap(plot, cmap=sns.color_palette("coolwarm", n_colors=40, desat=0.8), vmin=-1,vmax=2.5, ax = axes[i][j], yticklabels=10, center = 0.5,
                            cbar=k <1, cbar_ax=None if k else cbar_ax1, cbar_kws={'label': '$\mathregular{log_{10}}({V^2}/{Hz})$'})
                print(k<1)
            if i == 1:
                sns.heatmap(plot, cmap=sns.color_palette("coolwarm", n_colors=40, desat=0.8), vmin=-1, vmax=2.5, ax = axes[i][j], yticklabels=10, center = 0.5,
                            cbar=k <1 , cbar_ax=None if k else cbar_ax2, cbar_kws={'label': '$\mathregular{log_{10}}({V^2}/{Hz})$'})
                print(k < 1)
            old_ticks = axes[i][j].get_xticks()
            new_ticks = np.linspace(np.min(old_ticks), np.max(old_ticks), 5)
            new_ticks = new_ticks[range(1,4)]
            axes[i][j].set_xticks(new_ticks)
            axes[i][j].set_xticklabels(["25%", "50%", "75%"], rotation=0)
            if j == 0:
                axes[i][j].set_xticks([new_ticks[1]])
                axes[i][j].set_xticklabels(["~6hrs before"], rotation=0)
                axes[i][j].set_yticklabels([0, 10,20,30,40,50], rotation=0)
            if j > 1:
                axes[i][j].set_yticks([])
                if j <4:
                    axes[i][j].axvline(x=np.shape(ALL_data_GM)[1], c = "black", linewidth=4, linestyle='--')
            k = k + 1
        axes[i][j].invert_yaxis()
        cbar_ax1.yaxis.set_ticks_position('left')
        cbar_ax1.yaxis.set_label_position('left')
        cbar_ax2.yaxis.set_ticks_position('left')
        cbar_ax2.yaxis.set_label_position('left')

fig1.text(0.02, 0.5, 'Frequency (Hz)', ha='center', va='center', rotation='vertical', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.5, 0.05, 'Time (% seizure length)', ha='center', va='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.145, 0.9, 'Interictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.31, 0.9, 'Preictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.515, 0.9, 'Ictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.72, 0.9, 'Postictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
subplot_x=0.94
fig1.text(subplot_x, 0.3, 'Gray\nMatter', ha='center', va='center', fontdict = {'fontsize': fontsize*1.1})
fig1.text(subplot_x, 0.7, 'White\nMatter', ha='center', va='center',  fontdict = {'fontsize': fontsize*1.1})
fig1.text(0.5, 0.98, 'Spectrogram', ha='center', va='top', fontdict = {'fontsize': fontsize*1.1})

plt.savefig(outputfile_spectrogram, dpi = 100)








#PSDvsDistance
plt.rcParams['figure.figsize'] = (20.0, 6)
plt.rcParams['font.family'] = "serif"
sns.set_context("talk")
fontsize = 30

#making plot layout
fig1 = plt.figure(constrained_layout=False)
left = 0.07; right = 0.9
gs = fig1.add_gridspec(nrows=1, ncols=4, left=left, right=right, bottom = 0.15, top = 0.80 ,wspace=0.1, width_ratios=[1, 1, 1, 1])
axes = [None] * gs._ncols
for a in range(len(axes)):
    axes[a] = fig1.add_subplot(gs[:, a])
PSDvsDist_limits = [51,800]
PSDvsDistance = [interictal, preictal, ictal, postictal]
cbar_ax3 = fig1.add_axes([0.965, 0.15, .02, 0.65])
k = 0
for j in range(0, 4):
    plot = PSDvsDistance[k]
    plot = np.nanmean(plot, axis=2)
    plot = np.delete(np.array(plot), range(30, np.shape(plot)[0]), axis=0)
    plot = np.delete(plot, range(800, np.shape(plot)[1]), axis=1)
    plot = np.log10(plot)
    sns.heatmap(plot,cmap=sns.color_palette("Spectral_r", n_colors=100, desat=0.6), vmin=2,vmax=4.5, ax = axes[j],
                cbar=k < 1, cbar_ax=None if k else cbar_ax3, cbar_kws={'label': '$\mathregular{log_{10}}({V^2}/{Hz})$'})
    old_ticks = axes[j].get_xticks()
    new_ticks = np.linspace(np.min(old_ticks), np.max(old_ticks), 5)
    new_ticks = new_ticks[range(1,4)]
    axes[j].set_xticks([0,200,400,600,800])
    axes[j].set_xticklabels(["0", "2", "4", "6", "8"], rotation=0)
    if j ==0:
        axes[j].set_yticks([0, 10, 20, 30])
        axes[j].set_yticklabels([0, 10, 20, 30], rotation=0)
    if j > 0:
        axes[j].set_yticks([])
    k = k + 1
    axes[j].invert_yaxis()
cbar_ax3.yaxis.set_ticks_position('left')
cbar_ax3.yaxis.set_label_position('left')

fig1.text(0.02, 0.5, 'Frequency (Hz)', ha='center', va='center', rotation='vertical', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.5, 0.025, 'Distance (mm)', ha='center', va='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.16, 0.845, 'Interictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.375, 0.845, 'Preictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.59, 0.845, 'Ictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.81, 0.845, 'Postictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.5, 0.98, 'Power vs Distance', ha='center', va='top', fontdict = {'fontsize': fontsize*1.1})

axes[0].add_patch(Ellipse((90, 7.1), width=160, height=12,edgecolor='#0000cccc',facecolor='none',linewidth=5))
axes[0].add_patch(Ellipse((700, 7.1), width=160, height=12,edgecolor='#0000cccc',facecolor='none',linewidth=5))

axes[2].add_patch(Ellipse((90, 7.1), width=160, height=12,edgecolor='#cc0000cc',facecolor='none',linewidth=5))
axes[2].add_patch(Ellipse((700, 7.1), width=160, height=12,edgecolor='#cc0000cc',facecolor='none',linewidth=5))

plt.savefig(outputfile_PSDvsDistance, dpi = 100)






#SNR

plt.rcParams['figure.figsize'] = (20.0, 9.5)
plt.rcParams['font.family'] = "serif"
sns.set_context("talk")
fontsize = 30

#making plot layout
fig1 = plt.figure(constrained_layout=False)
left = 0.07; right = 0.88
gs = fig1.add_gridspec(nrows=1, ncols=5, left=left, right=right, bottom = 0.10, top = 0.80 ,wspace=0.0, width_ratios=[1,0.1, 1, 1.2, 1])
axes = [None] * gs._ncols
for a in range(len(axes)):
    axes[a] = fig1.add_subplot(gs[:, a])

k = 0
#line colors
legend_labels = breaks[:]
names = []
for f in range(len(legend_labels)):
    if f == 0:
        names.append("= {0}".format(legend_labels[f]))
    if f > 0:
        names.append("({0}, {1}]".format(legend_labels[f - 1], legend_labels[f]))
names.append("≥ {0}]".format(legend_labels[len(legend_labels) - 1]))
cols = ["Time"]
for f in range(len(names)):
    cols.append(names[f])
line_colors = sns.color_palette("Blues", n_colors=len(names)-2, desat=0.6)
black = (0, 0, 0)
red = (0.7, 0, 0)
line_colors.append(red)
line_colors.insert(0, black)
for j in range(0, 5):
    time = np.linspace(0, 100, num=len(SNR_all_patients_mean), endpoint=True)
    time = time.reshape(time.shape[0],1)
    SNR_to_plot = SNR_all_patients_mean[:,:,k]
    SNR_to_plot = np.concatenate([time, SNR_to_plot], axis=1)
    SNR_mean_df_wide = pd.DataFrame(SNR_to_plot, columns = cols)
    SNR_mean_df_long = pd.melt(SNR_mean_df_wide, ['Time'])
    SNR_mean_df_long.columns = ['Time', 'Distance (mm)', 'Values']
    if j == 1:#create a blank space between interictal and preictal
        axes[j].remove()
    if j != 1:
        sns.lineplot(x='Time', y='Values', hue='Distance (mm)', data=SNR_mean_df_long, ax = axes[j],
                     palette=line_colors)
        axes[j].set_xlim(0,100)
        axes[j].set_ylim(0, 60)
        axes[j].set(xlabel='', ylabel='')
        if j < 4:
            axes[j].get_legend().set_visible(False)
        if j == 4:
            handles, labels = axes[j].get_legend_handles_labels()
            axes[j].legend(handles=handles[1:], labels=labels[1:], title="Distance (mm)",
                              bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        old_ticks = axes[j].get_xticks()
        new_ticks = np.linspace(np.min(old_ticks), np.max(old_ticks), 5)
        new_ticks = new_ticks[range(1, 4)]
        axes[j].set_xticks(new_ticks)
        axes[j].set_xticklabels(["25%", "50%", "75%"], rotation=0)
        if j == 0:
            axes[j].set_xticks([new_ticks[1]])
            axes[j].set_xticklabels(["~6hrs before"], rotation=0)
         #   axes[i][j].set_yticklabels([0, 10, 20, 30, 40, 50], rotation=0)
        if j > 1:
            axes[j].set_yticks([])
            if j < 4:
                axes[j].axvline(x=np.shape(ALL_data_GM)[1], c="black", linewidth=4, linestyle='--')
        k = k + 1



fig1.text(0.02, 0.5, 'SNR', ha='center', va='center', rotation='vertical', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.5, 0.025, 'Time (% seizure length)', ha='center', va='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.155, 0.845, 'Interictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.36, 0.845, 'Preictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.565, 0.845, 'Ictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.78, 0.845, 'Postictal', ha='center', fontdict = {'fontsize': fontsize*0.7})
fig1.text(0.5, 0.97, 'Signal to Noise Ratio', ha='center', va='top', fontdict = {'fontsize': fontsize*1.1})

fig1.text(0.08, 0.5, '$\mathregular{SNR}=(\\frac{Power_{s}}{Power_{i}})$', ha='left', fontdict = {'fontsize': fontsize*1.2})

plt.savefig(outputfile_SNR, dpi = 100)



