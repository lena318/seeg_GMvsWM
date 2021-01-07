#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:23:28 2020

@author: arevell
"""

#combine subcortical segmentation (FIST) and regular segmentation (FAST)

import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt



def show_slices(img_data, low = 0.33, middle = 0.5, high = 0.66):
    """ Function to display row of image slices """
    slices1 = [   img_data[:, :, int((img_data.shape[2]*low)) ] , img_data[:, :, int(img_data.shape[2]*middle)] , img_data[:, :, int(img_data.shape[2]*high)]   ]
    slices2 = [   img_data[:, int((img_data.shape[1]*low)), : ] , img_data[:, int(img_data.shape[1]*middle), :] , img_data[:, int(img_data.shape[1]*high), :]   ]
    slices3 = [   img_data[int((img_data.shape[0]*low)), :, : ] , img_data[int(img_data.shape[0]*middle), :, :] , img_data[int(img_data.shape[0]*high), :, :]   ]
    slices = [slices1, slices2, slices3]
    plt.style.use('dark_background')
    fig = plt.figure(constrained_layout=False, dpi=300, figsize=(5, 5))
    gs1 = fig.add_gridspec(nrows=3, ncols=3, left=0, right=1, bottom=0, top=1, wspace=0.00, hspace = 0.00)
    axes = []
    for r in range(3): #standard
        for c in range(3):
            axes.append(fig.add_subplot(gs1[r, c]))
    r = 0; c = 0
    for i in range(9):
        if (i%3 == 0 and i >0): r = r + 1; c = 0
        axes[i].imshow(slices[r][c].T, cmap="gray", origin="lower")
        c = c + 1
        axes[i].set_xticklabels([])
        axes[i].set_yticklabels([])
        axes[i].set_xticks([])
        axes[i].set_yticks([])
        axes[i].axis("off")







# Progress bar function
def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = "X", printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


