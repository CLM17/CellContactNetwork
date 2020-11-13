#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:48:21 2020

@author: lukasvdh

This script does a brute force stitching. Use it if stitching in FIJI does not work.
Brute force stitching = sticking the tiles next to each other without caluclating the overlap between images.

"""

import os
from skimage import io
from Helpers.brute_force_stitch_helpers import paste_tiles


root = 'M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/WKS025/10x'
well = 'D03'
w = 8

fused = paste_tiles(root, well, w)

output_path = os.path.join(root, well, well+'_fused.tif')

if os.path.isfile(output_path):
    a = input(well + "_fused.tif already exists. Do you want to overwite? (y/n)")
    if a != "y":
        raise ValueError("User terminated the script. Fused image is not saved.")

print('Saving the fused image...')
io.imsave(output_path, fused)
print('Fused image saved succesfully in well folder.')