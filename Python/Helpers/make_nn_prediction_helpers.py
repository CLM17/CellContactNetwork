# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 09:28:39 2020

@author: clmvandenheuve
"""

import os
import numpy as np
from skimage import io
from model import *  # unet model in unet folder
from data import *   # unet data in unet folder

def crop_img(P, w, ol):
    '''
    This function find the limits of an image dimension s.t. tiles fit inside it.
    The operation is done for one dimension only.
    P = number of pixels of the input image, in one dimension.
    w = width of tile (e.g. 512).
    ol = overlap of tiles (e.g. 100 pixels)
    '''
    
    numtiles = (P - w) / (w - ol)
    numtiles_rounded = np.floor(numtiles)
    
    num_included = numtiles_rounded * (w - ol) + w 
    num_cropped = P - num_included
    
    # Notify the user if more than 25 pixels are cropped on both sides:
    if num_cropped > 50:
        print("WARNING: You are losing more than 25 pixels on both sides of the image.\n Check your image afterwards to see if it still contains the whole well.\n")    
    
    # calculate margins (min and max)
    m_min = np.floor(num_cropped / 2)
    m_max = num_cropped - m_min
    
    min_ = int(m_min)
    max_ = int(P - m_max)
    
    return min_, max_, numtiles_rounded+1

#%%

def create_tiles(cropped_ch, xnumtiles, ynumtiles, w, ol):
    '''
    Create tiles with overlap.
    '''
    xtile_positions = np.arange(0, xnumtiles) * (w - ol)
    ytile_positions = np.arange(0, ynumtiles) * (w - ol)
    
    tiles = []
    for x in xtile_positions:
        for y in ytile_positions:
            x = int(x)
            y = int(y)
            tiles.append(cropped_ch[y:y+w, x:x+w, :])
            
    return tiles

#%%
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

#%%
def combine_overlapping_tiles(results, w, ol, xnumtiles, ynumtiles, cr_img):
    
    # Combine tiles into columns
    
    rg = np.linspace(0,1,ol).reshape(ol,1) # range
    vert_combiner_bottom = np.vstack( (np.matlib.repmat(rg,1,w), np.ones((w-ol, w))) )
    vert_combiner_top = np.vstack( ( np.ones((w-ol, w)), np.matlib.repmat(np.flipud(rg),1,w)) )
    vert_combiner_center = np.vstack( ( (np.matlib.repmat(rg,1,w), np.ones((w-2*ol, w)), np.matlib.repmat(np.flipud(rg),1,w)) ) )
    
    xtile_positions = np.arange(0, xnumtiles) * (w - ol)
    ytile_positions = np.arange(0, ynumtiles) * (w - ol)
    
    Ncr, Mcr, Ccr = cr_img.shape
    col_list = []
    c = 0
    for x in xtile_positions:
        col = np.zeros((Ncr, w))
        for y in ytile_positions:
            x = int(x)
            y = int(y)
            
            pred = results[c].reshape(w,w)
            pred = pred / np.max(pred)
            nbh = col[y:y+w,:]
            if y==0:
                nbh = nbh + vert_combiner_top * pred
            elif y==ytile_positions[-1]:
                nbh = nbh + vert_combiner_bottom * pred
            else:
                nbh = nbh + vert_combiner_center * pred
                
            col[y:y+w,:] = nbh
            c = c + 1
                
        col_list.append(col)
    
    # combine columns into full image
    
    rg = np.linspace(0,1,ol).reshape(1,ol)
    hor_combiner_right = np.hstack( (np.matlib.repmat(rg,Ncr,1), np.ones((Ncr, w-ol))) )
    hor_combiner_left = np.hstack( (np.ones((Ncr, w-ol)), np.matlib.repmat(np.fliplr(rg),Ncr,1)) )
    hor_combiner_center = np.hstack( (np.matlib.repmat(rg,Ncr,1), np.ones((Ncr, w-2*ol)), np.matlib.repmat(np.fliplr(rg),Ncr,1)) )
    
    prediction = np.zeros((Ncr,Mcr))
    c = int(0)
    for x in xtile_positions:
        x = int(x)
        nbh = prediction[:,x:x+w]
        if x==0:
            nbh = nbh + hor_combiner_left*col_list[c]
        elif x==xtile_positions[-1]:
            nbh = nbh + hor_combiner_right*col_list[c]
        else:
            nbh = nbh + hor_combiner_center*col_list[c]
        prediction[:,x:x+w] = nbh
        c = c+1
    
    return prediction

#%%
def make_nn_prediction_from_fused(well_folder, well, w, ol, nr_color_channels, nr_classes, unet_folder, model_name):

    print('Reading the fused image...')
    img_path = os.path.join(well_folder, well + '_fused_RGB.tif')
    img = io.imread(img_path)
    
    N,M,C = img.shape
    print('Read the image')
    
    xmin, xmax, xnumtiles = crop_img(M, w, ol)
    ymin, ymax, ynumtiles = crop_img(N, w, ol)
    cr_img = img[ymin:ymax, xmin:xmax, :]
    tiles = create_tiles(cr_img, xnumtiles, ynumtiles, w, ol)
    
    # predict the tiles
    results = predict_tiles(tiles, unet_folder, model_name, nr_classes)
    #saveResult(predict_path, results)
    
    # combine predicted tiles into cols
    print('\nCombining overlapping tiles into a fused image with linear blending...')
    prediction = combine_overlapping_tiles(results, w, ol, xnumtiles, ynumtiles, cr_img)
    
    # Make fused image into original size
    full_prediction = np.zeros((N,M))
    full_prediction[ymin:ymax, xmin:xmax] = prediction
    
    return full_prediction / np.max(full_prediction)

