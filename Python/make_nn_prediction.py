# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 09:28:08 2020

@author: clmvandenheuve
"""

#------------------------------SPECIFY PARAMETERS-----------------------------

# Overlapping tile parameters
w = 512  # width / height of 1 tile
ol = 100 # overlap between tiles
root = 'M:\\tnw\\bn\\dm\\Shared\\Lukas\\BEP\\Experiments\\JJ005\\20x'
well_list = ['C02']

# Unet model parameters
nr_color_channels = 3
nr_classes = 2
output_title = '_boundaries.tif'
unet_folder = 'M:\\tnw\\bn\\dm\\Shared\\Lukas\\BEP\\CellContactNetwork\\Unet Master multiclass'
model_name = 'cos_boundaries_final.hdf5'

#----------------------------------START CODE---------------------------------

from sys import path
path.append(unet_folder)
import os
from skimage import io
from Helpers.make_nn_prediction_helpers import make_nn_prediction_from_fused
from skimage import img_as_ubyte

for well in well_list:
    print('Predicting well ' + well)
    
    well_folder = os.path.join(root, well)
    full_prediction = make_nn_prediction_from_fused(well_folder, well, w, ol, nr_color_channels, nr_classes, unet_folder, model_name)
    
    # save the result
    print('\nSaving the fused prediction...')
    io.imsave(os.path.join(well_folder, well+output_title), img_as_ubyte(full_prediction))
    print('\nBoundary prediction is successfully saved.')
