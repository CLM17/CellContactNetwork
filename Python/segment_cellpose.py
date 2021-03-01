# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 14:02:51 2021

@author: lukasvandenheu
"""
import os
from skimage import io
from Helpers.segment_cellpose_helpers import segment_fused_image_with_cellpose

# ----------------------------SPECIFY PARAMETERS-------------------------------

well = 'B02'
root = r'M:\tnw\bn\dm\Shared\Lukas\NetworkAnalysis\CellContactNetwork\Cellpose'
well_folder = os.path.join(root,well)
[patch_height, patch_width] = [512,512]
edge_thickness = 60           # size of edge region in pixels
similarity_threshold = 0.8    # overlapping cells which are similar will be merged
cell_size_threshold = 100     # minimal cell size in pixels

# DEFINE CELLPOSE PARAMETERS
# model_type='cyto' or model_type='nuclei'

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
# channels = [0,0]
# IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
# channels = [0,0] # IF YOU HAVE GRAYSCALE
# channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
# channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus

# or if you have different types of channels in each image
# channels = [[2,3], [0,0], [0,0]]

# if diameter is set to None, the size of the cells is estimated on a per image basis
# you can set the average cell `diameter` in pixels yourself (recommended) 
# diameter can be a list or a single number for all images

channels = [1,2] # R=cytoplasm and G=nucleus
diameter = 97

# -------------------------------START CODE------------------------------------

segmented = segment_fused_image_with_cellpose(root, well, diameter, channels,
                                              edge_thickness=edge_thickness, similarity_threshold=similarity_threshold,
                                              cell_size_threshold=cell_size_threshold, patch_height=patch_height, patch_width=patch_width)

output_path = os.path.join(root, well+'_fused_mask.tif')
io.imsave(output_path, segmented)