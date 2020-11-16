# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 13:53:24 2020

@author: clmvandenheuve
"""

import numpy as np
import os
from skimage import io

mask_folder = 'data/hela/Masks'
img_folder = 'data/hela/train_test'
output_folder = 'data/hela/train/'

file_list = os.listdir(mask_folder)

for file in file_list:
    if not '.tif' in file:
        file_list.remove(file)
print(file_list)

c = 0
for file in file_list:
    mask = io.imread(os.path.join(mask_folder,file))
    img = io.imread(os.path.join(img_folder,file))
    ch_MEM = img[:,:,0]
    fname = str(c)+'.png'
    io.imsave(os.path.join(output_folder, 'image', fname), ch_MEM / np.max(ch_MEM))
    io.imsave(os.path.join(output_folder, 'label', fname), mask / np.max(mask))
    c = c+1