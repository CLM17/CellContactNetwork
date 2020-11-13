# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:09:36 2020

@author: LocalAdmin
"""

from skimage.filters.thresholding import threshold_triangle,threshold_yen,threshold_niblack, threshold_local, threshold_otsu, threshold_li, threshold_minimum, threshold_multiotsu, threshold_mean, threshold_isodata
from skimage.morphology import binary_closing, binary_opening, binary_dilation, binary_erosion
from skimage.draw import line, disk, circle_perimeter
from network_creator_helpers import get_well_locations, get_image

import numpy as np
from skimage import io
import os
import matplotlib.pyplot as plt

well = 'D04'
folder = 'Data/WKS023/2020-09-09'
well_folder = os.path.join(folder, 'well '+well)

xc_well, yc_well, diameter_well = get_well_locations(folder, well)
full_img = get_image(well_folder, well)

ch_DAPI = full_img[:,:,2] / np.max(full_img[:,:,2])
ch_MEM = full_img[:,:,0] / np.max(full_img[:,:,0])

threshold = threshold_li(ch_MEM)
binary = ch_MEM > threshold
io.imshow(binary)
io.imsave(os.path.join(well_folder, well+'_binary.tif'),binary)

'''
plt.figure()
io.imshow(ch_MEM)
plt.figure()
#io.imshow(ch_MEM > threshold)
#io.imsave(os.path.join(well_folder,'binary.tif'), ch_MEM>threshold)
'''