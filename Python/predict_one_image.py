# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 14:57:56 2020

@author: clmvandenheuve
"""

# Unet model parameters
nr_color_channels = 3
nr_classes = 2
output_title = '_predict.tif'
unet_folder = 'M:\\tnw\\bn\\dm\\Shared\\Lukas\\BEP\\CellContactNetwork\\Unet Master multiclass'
model_name = 'hela_boundaries_final.hdf5'

from sys import path
path.append(unet_folder)
import os
from skimage import io
from skimage import img_as_ubyte

import os
import numpy as np
from skimage import io
from model import *  # unet model in unet folder
from data import *   # unet data in unet folder

def predict_tiles(tiles, unet_folder, model_name, nr_classes):
        
    num_test_images = len(tiles)

    # predict the tiles
    model = unet(nb_classes = nr_classes,
                 pretrained_weights = os.path.join(unet_folder, 'models', model_name),
                 input_size = tiles[0].shape)
    
    print('\nMaking a prediction for all '+str(num_test_images)+' tiles...\n')
    testGene = predictGenerator(tiles, target_size = tiles[0].shape)
    results = model.predict_generator(testGene, num_test_images, verbose=1)
    
    return results

img = io.imread("M:\\tnw\\bn\\dm\\Shared\\Lukas\\BEP\\Experiments\\Test\\small\\D02\\D02_fused_RGB_512.tif")
predict = predict_tiles([img], unet_folder, model_name, nr_classes)
io.imsave("M:\\tnw\\bn\\dm\\Shared\\Lukas\\BEP\\Experiments\\Test\\small\\D02\\D02_boundaries.tif", predict[0])