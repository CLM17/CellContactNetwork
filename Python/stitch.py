# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 08:52:26 2020

@author: LocalAdmin
"""

import os
import numpy as np
from skimage import io
from Helpers.network_creator_helpers import distance_transform_watershed
from skimage.filters.thresholding import threshold_otsu

root = 'Data/WKS024'
well = 'B02'
w = 8

well_folder = os.path.join(root, well)
tile_folder = os.path.join(well_folder, 'tiles')

file_folder = os.path.join(tile_folder, 'tile_36.tif')
img = io.imread(file_folder)

ch_dapi = img[:,:,2] / np.max(img[:,:,2])
ch_marker = img[:,:,0] / np.max(img[:,:,0])

#th = threshold_otsu(ch_marker)
#io.imshow(ch_marker>th, cmap='Greys')

th = threshold_otsu(ch_dapi)
th_dapi = ch_dapi > th

labels, markers = distance_transform_watershed(th_dapi, th_dapi, connectivity=np.ones((3,3)), sigma=1, min_distance=10)

io.imshow(labels,cmap='Greys')
output_path_labels = os.path.join(well_folder, 'dist_trans_wts.tif')
output_path_markers = os.path.join(well_folder, 'dist_trans_wts_markers.tif')
io.imsave(output_path, labels)
#distance_transform_watershed(th_dapi, th_marker, connectivity=np.ones((3,3)), sigma=1, min_distance=10)
