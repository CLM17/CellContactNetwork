# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 08:52:26 2020

@author: LocalAdmin
"""

import numpy as np
from skimage import io
from skimage import img_as_uint

import os

#%%

def make_spiral_grid(width, height):
    '''
    The high content microscope images a well in spiral shape.
    This function makes a spiral-shaped grid, that we can use to find the locations of the image.
    '''
    NORTH, S, W, E = (0, -1), (0, 1), (-1, 0), (1, 0) # directions
    turn_right = {NORTH: E, E: S, S: W, W: NORTH} # old -> new direction
    
    if width < 1 or height < 1:
        raise ValueError
    
    x, y = int(np.ceil(width/2)-1), int(np.ceil(height/2)-1) #- 1 # start near the center
    dx, dy = NORTH # initial direction
    matrix = [[None] * width for _ in range(height)]
    count = 0
    while True:
        matrix[y][x] = count # visit
        count += 1
        # try to turn right
        new_dx, new_dy = turn_right[dx,dy]
        new_x, new_y = x + new_dx, y + new_dy
        if (0 <= new_x < width and 0 <= new_y < height and
            matrix[new_y][new_x] is None): # can turn right
            x, y = new_x, new_y
            dx, dy = new_dx, new_dy
        else: # try to move straight
            x, y = x + dx, y + dy
            if not (0 <= x < width and 0 <= y < height):
                print(np.asarray(matrix))
                a = input('This is the spiral matrix. Do you want to continue? (y/n)')
                if a!='y':
                    raise ValueError('You aborted the stitching preparation')

                return np.asarray(matrix) # nowhere to go

#%%
def make_column_grid(width, height):
    '''
    The FIJI stitching algorithm can stitch images in a column-by-column grid.
    This function makes a column-shaped grid.
    '''
    matrix = np.arange(width * height).reshape(height,width).transpose()
    return matrix

#%%
def convert_array_to_string(a):
    '''
    This function takes as input a 2D array a.
    It outputs the same array, but integers are strings formatted as '00', '01', etc.
    '''
    if len(a.shape) == 1:
        a = [a]
    a_as_str = []
    for row in a:
        str_row = []
        for num in row:
            str_row.append("%02d" % num)
        a_as_str.append(str_row)        

    return a_as_str

#%%
def prepare_dataset_for_stitching(raw, root, well_list, w, n_channels):
    '''
    This function combines DAPI and MEMBRITE channels of acquired images and stores them under the correct filename,
    so they can be used for stitching in imageJ.
    
    Inputs:
        target_folder: path to folder of the experiment. It is named with the experiment data (e.g. 2020-09-01)
        img_folder: name of folder with the raw images (e.g. 'AcquireOnly.V3_09-01-20_04;16;23')
        well_rows: list of well rows that were imaged, e.g. ['B', 'D']
        well_cols: list of well columns that were imaged, e.g. [2, 3, 4]
        w = number of images in the stitch, e.g. 8
        n_channels = number of channels
    '''

    # Define a spiral-shaped grid
    spiral_grid = make_spiral_grid(w,w)
    column_grid = make_column_grid(w, w)
    column_grid = convert_array_to_string(column_grid)

    nrs = convert_array_to_string( np.arange(w ** 2) )[0]

    for well in well_list:
            
        raw_folder = os.path.join(raw, well)
        well_folder = os.path.join(root, well)
        tile_folder = os.path.join(well_folder, 'tiles')
        
        file_list = os.listdir( raw_folder )

        # get first part of file name (is the same for all files in folder)
        splitted_file = file_list[1].split('_')
        img_info = splitted_file[0] + '_' + splitted_file[1] + '_'

        # Create tile new directory for this well in root (if it didn't exist already)
        if os.path.exists( well_folder ):
            ans = input("A folder for well "+well+" already exists in root. Do you want to continue? (y/n)")
            if ans!="y":
                raise ValueError("You stopped the stitching preparation.")
        else:
            os.mkdir( well_folder )
            
        if not os.path.exists( tile_folder ):
            os.mkdir( tile_folder )

        print("Busy copying the images of well "+well+"...")

        count = 0
        for nr in nrs:
            
            # load the first channel
            fname0 = img_info + well + 'f' + nr + 'd0.tif'
            ch0 = io.imread( os.path.join(raw_folder, fname0) )[:,:,0]
            # initiate empty image with the right dimentsions
            N,M = ch0.shape
            img = np.zeros((n_channels,N,M))
            img[0,:,:] = ch0 / np.max(ch0)

            for ch_nr in range(1,n_channels):
                
                fname = img_info + well + 'f' + nr + 'd' + str(ch_nr) + '.tif'
                ch = io.imread( os.path.join(raw_folder, fname) )[:,:,0]
                img[ch_nr,:,:] = ch / np.max(ch)
            
            #fname_marker2 = img_info + well + 'f' + nr + 'd2.tif'
            #img_marker2 = io.imread( os.path.join(raw_folder, fname_marker2) )[:,:,0]

            r,c = np.where(spiral_grid == count)
            tile_index = column_grid[r[0]][c[0]]
            fname = 'tile_' + tile_index + '.tif'
            io.imsave(os.path.join(tile_folder, fname), img_as_uint(img))
                        
            count = count + 1
            
