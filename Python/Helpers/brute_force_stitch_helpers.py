#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:25:23 2020

@author: lukasvdh
"""
import numpy as np
import os
from skimage import io, img_as_ubyte
    
def get_tile_dimensions (tile_folder, file_list):
    
    '''
    Outputs an array with the dimensions of the first tile in tile_folder.
    '''
    
    # pick the first tile and measure its dimensions
    for file in sorted(file_list):
        if '.tif' in file:
            file_path = os.path.join(tile_folder, file)
            tile = io.imread(file_path)
            break
    
    return tile.shape

def paste_tiles(root, well, w):
    '''
    This function pastes tiles next to each other in the right order. (= "brute force stitching").
    '''
    
    print('Preparing brute-force stitching...')
    
    well_folder = os.path.join(root, well)
    tile_folder = os.path.join(well_folder, 'tiles')
    file_list = os.listdir(tile_folder)
    
    n,m,z = get_tile_dimensions (tile_folder, file_list)
    
    # initialize an empty fused image
    fused = np.zeros((w * n, w * m, z))
    # initialize counters to keep track of the (x,y) position
    xc, yc = [0,0]
    
    print('Stitching the tiles...')
    
    for file in sorted(file_list):
        if '.tif' in file:
            file_path = os.path.join(tile_folder, file)
            tile = io.imread(file_path)
            tile = tile / np.max(tile)
            
            fused[ n*yc:n*(yc+1), m*xc:m*(xc+1), : ] = tile
            
            yc = yc + 1
            if yc == w:
                yc = 0
                xc = xc + 1
    
    return img_as_ubyte(fused)