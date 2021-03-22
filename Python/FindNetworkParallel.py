# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 10:04:40 2021

@author: lukasvandenheu
"""


# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 17:50:27 2021

@author: lukasvandenheuvel
"""

#%%
import numpy as np
from skimage import io
from skimage import measure
# from skimage.morphology import binary_closing, binary_opening, binary_dilation, binary_erosion
# from skimage.draw import line, circle, circle_perimeter
# from skimage.segmentation import watershed
# from skimage import img_as_ubyte, img_as_uint
# from skimage.filters import gaussian
# from skimage.feature import peak_local_max

from scipy import ndimage as ndi
from scipy import sparse
from scipy.io import savemat

from numba import jit
from joblib import Parallel, delayed
import multiprocessing
import os
# from timeit import default_timer as timer
# import matplotlib.pyplot as plt

import tkinter as tk
import tkinter.filedialog
#from tkinter import Label, W, IntVar, Button, Checkbutton, mainloop

from cellpose import models
# from progress.bar import Bar

use_GPU = models.use_gpu()
print('>>>> GPU activated? %d'%use_GPU)

#%%

def enlarge_fused_image(fused, patch_height=512, patch_width=512, overlap=256):
    
    '''
    This function makes the fused image larger s.t. an uneven integer number of
    patches fits into it.
    
    INPUTS
    ------
        fused (np array)
        Fused (whole well) image.
        
        overlap (int)
        Overlap of subimages in pixels. Default is 256.
        
        patch_width, patch_height
        Size of patches. Default is 512.
        
    OUTPUTS
    -------
        new_fused (np array)
        Enlarged fused image.
        
        patch_locations (list with 2 elements)
        patch_locations[0] is a numpy array with the y-coordinates of the subimages in pixels.
        patch_locations[1] is a numpy array with the x-coordinates of the subimages in pixels.
    '''
    
    [m,n] = [patch_height, patch_width]
    [M,N,C] = np.shape(fused)
    
    num_m = np.ceil((M - m) / (m - overlap) + 1).astype(int) # number of patches that fit on y-axis
    num_n = np.ceil((N - n) / (n - overlap) + 1).astype(int) # number of patches that fit on x-axis
    
    # make sure that num_m and num_n are uneven,
    # s.t. the last patch is NOT an overlapping patch
    if (num_m%2 == 0):
        num_m = num_m + 1
    if (num_n%2 == 0):
        num_n = num_n + 1

    new_M = (num_m - 1) * (m - overlap) + m # new fused image height
    new_N = (num_n - 1) * (n - overlap) + n # new fused image width
    
    new_fused = np.zeros([new_M,new_N,C], dtype='uint8')
    new_fused[0:M,0:N,:] = fused
    
    patch_locations = []
    patch_locations.append(np.arange(0,new_M-m+1,m-overlap)) # y-locations of patches (in pixels)
    patch_locations.append(np.arange(0,new_N-n+1,n-overlap)) # x-locations of patches (in pixels)
    
    return [new_fused, patch_locations]

#%%
def fused_to_patches(fused, patch_locations, patch_height=512, patch_width=512):
    '''
    This function takes as input an (enlarged) fused image.
    It outputs a list of patches which are ordered in a column-to-column grid.
    The locations of the patches in the fused image are specified by patch_locations.
    '''
    [n,m] = [patch_height, patch_width]
    patch_list = []
    for c_n in patch_locations[1]:     # c_n is n-coordinate in pixels
        for c_m in patch_locations[0]: # c_m is m-coordinate in pixels
            patch = fused[c_m:c_m+m,c_n:c_n+n,:]
            patch_list.append(patch)
    
    return patch_list

#%%
def calculate_img_similarity(img1, img2):
    '''
    This function calculates what percentage of
    img1 and img2 are the same.
    Img1 and img2 are boolean images.
    '''
    equal_pixels = np.logical_and(img1,img2)
    return np.sum(equal_pixels) / min([np.sum(img1), np.sum(img2)])

#%%
def find_cell_values_which_overlap(cell_on_overlap, combined_overlap, similarity_threshold):
    '''
    This function finds which cells on the combined_overlap image
    overlap with the boolean image cell_on_overlap.
    '''
    overlapping_cell_values = cell_on_overlap * combined_overlap
    unique = np.unique(overlapping_cell_values)
    values_with_sufficient_overlap = []
    
    # loop through cells which overlap with the cell on the egde
    # and calculate similarity
    for overlapping_value in unique[1:]:
        cell_on_mask = (combined_overlap==overlapping_value)
        if calculate_img_similarity(cell_on_mask, cell_on_overlap) > similarity_threshold:
            values_with_sufficient_overlap.append(overlapping_value)
            
    return values_with_sufficient_overlap

#%%
def store_overlapping_cells(values_on_edge, patch_ol, combined_ol, similarity_threshold):
    '''
    This function stores which cells on combined_ol overlap with the
    cells in the edge region of patch_ol (in a dictionary).
    Keys in the 'overlap_dict' dictionary are the values of cells which lie
    in the edge region of patch_ol.
    Corresponding values are the cells on combined_ol which sufficiently
    overlap with the cells.
    Example: overlap_dict = {23: [34,45]}
             means that there is a cell (nr 23) on path_ol in the edge region.
             It overlaps with cells 34 and 45. These are probably two half cells.
    '''
    # Initialise dictionary to store old and new values
    overlap_dict = {}
    # Collect values to be removed in the replace dictionary
    for value in values_on_edge:
        cell_shape = (patch_ol == value)
        values_overlapping = find_cell_values_which_overlap(cell_shape, combined_ol, similarity_threshold)
        overlap_dict[value] = values_overlapping   
        
    return overlap_dict

#%%
def replace_overlapping_cells(overlap_dict, patch_ol, combined_ol, max_value):
    '''
    This function replaces cells on the edge of combined_ol with 
    overlapping cells on the edge of patch_ol.
    The cells on patch_ol are correctly predicted, so they should
    replace the wrongly predicted cells on combined_ol.
    '''
    combined_ol_new = np.copy(combined_ol)
    for value, values_to_remove in overlap_dict.items():
        # remove value which overlaps with a cell on patch_ol
        for old_value in values_to_remove:
            combined_ol_new[np.where(combined_ol_new==old_value)] = 0
        # add the cells on mask_ol to the image
        if len(values_to_remove) > 0:     
            max_value = max_value + 1
            combined_ol_new[np.where(patch_ol==value)] = max_value
            
    return combined_ol_new

#%%
def align_patches(patch1, patch2, patch_ol, edge_thickness, similarity_threshold, max_value, axis):
    
    '''
    This function aligns patch1 and patch2 based on the overlap patch_ol.
    '''
    
    m,n = np.shape(patch1)
    edge_size = int(edge_thickness/2)
    
    # Get max value (values of new overlapping cells are always larger than max_value)
    # max_value = np.max(patch2)
    
    # Get a list of cells that lie on the edge region of the overlapping patch (values_on_edge)
    # Create combined overlap by pasting the 2 patches together (combined_ol)
    if (axis==0):
        center = int(m / 2)
        values_on_edge = np.unique( patch_ol[center-edge_size:center+edge_size+1,:] )
        combined_ol = np.concatenate([patch1[center:m,:],patch2[0:center,:]],axis=0)
    elif (axis==1):
        center = int(n / 2)
        values_on_edge = np.unique( patch_ol[:,center-edge_size:center+edge_size+1] )
        combined_ol = np.concatenate([patch1[:,center:n],patch2[:,0:center]],axis=1)
    else:
        raise ValueError('Invalid choice for axis. Please choose either 0 or 1.')
    
    # Remove 0 (=background) from the list
    values_on_edge = np.delete(values_on_edge, np.where(values_on_edge==0))    

    # Find overlapping cells 
    overlap = store_overlapping_cells(values_on_edge, patch_ol, combined_ol, similarity_threshold)
    combined_ol_new = replace_overlapping_cells(overlap, patch_ol, combined_ol, max_value)
    
    # Update patches
    if (axis==0):
        patch1_new = np.concatenate([patch1[0:center,:],combined_ol_new[0:center,:]],axis=0)
        patch2_new = np.concatenate([combined_ol_new[center:m,:],patch2[center:m,:]],axis=0)
    elif (axis==1):
        patch1_new = np.concatenate([patch1[:,0:center],combined_ol_new[:,0:center]],axis=1)
        patch2_new = np.concatenate([combined_ol_new[:,center:n],patch2[:,center:n]],axis=1)
    
    return patch1_new, patch2_new

#%%
def create_overlapping_columns_parallel(nn, mask_list, patch_locations, edge_thickness, similarity_threshold):
    
    patches_in_column = []
    num_patches_mm = np.size(patch_locations[0])

    # loop over rows in steps of 2 (avoid the last patch in column)
    for mm in range(0,num_patches_mm-1,2):

        # Get patch1, patch2 and patch_ol
        patch_nr = nn * num_patches_mm + mm
        if mm==0: # if we are at the top of the column
            patch1 = np.copy(mask_list[patch_nr])         # patch1 (upper patch) is new patch from patch_list
            patches_in_column.append(patch1)
        else:
            patch1 = np.copy(patches_in_column[-1])       # patch1 (upper patch) is previously processed patch
        patch2 = np.copy(mask_list[patch_nr + 2])         # patch2 is lower patch
        patch_ol = np.copy(mask_list[patch_nr + 1])       # patch_ol is overlapping patch

        # Increase cell values on patch2 with the max value of patch1
        patch2[np.where(patch2>0)] = patch2[np.where(patch2>0)] + np.max(patch1)
        
        # Align patch1 and patch2 using overlap
        max_value = np.max(patch2)
        patch1_new,patch2_new = align_patches(patch1, patch2, patch_ol, edge_thickness, similarity_threshold, max_value, axis=0)
        patches_in_column[-1] = patch1_new     # overwrite first patch
        patches_in_column.append(patch2_new)   # append new patch

    # Combine patches into a column
    aligned_patches = np.concatenate(patches_in_column,axis=0)
        
    return aligned_patches

#%%
def make_all_cell_values_unique(overlapping_columns):
    '''
    This function increases the cell values on all even columns with the  
    max cell value of the previous column, s.t. all cell valyues are unique.
    '''
    num_columns = len(overlapping_columns)
    max_value = 0
    for nn in range(2,num_columns,2):
        col = overlapping_columns[nn]
        col[np.where(col>0)] = col[np.where(col>0)] + np.max(overlapping_columns[nn-2])
        overlapping_columns[nn] = col
        if np.max(col) > max_value:
            max_value = np.max(col)
        
    return overlapping_columns,max_value

#%%
def align_overlapping_columns_parallel(nn, overlapping_columns, edge_thickness, similarity_threshold, max_value):

    patch1 = np.copy(overlapping_columns[nn])      # left column
    patch2 = np.copy(overlapping_columns[nn+2])    # right column
    patch_ol = np.copy(overlapping_columns[nn+1])  # overlapping column
    
    # Align patch1 and patch2 using overlap
    patch1_new,patch2_new = align_patches(patch1, patch2, patch_ol, edge_thickness, similarity_threshold, max_value, axis=1)
    # Paste the two new patches into one
    combined_patch = np.concatenate([patch1_new, patch2_new], axis=1)
            
    return combined_patch

#%%
def concatenate_overlapping_columns(aligned_overlapping_columns, new_M, new_N, patch_width):
    '''
    This function creates a fused image from a list of aligned overlapping columns.
    '''
   
    fused_mask = np.zeros((new_M,new_N))
    
    x = 0
    for i,col in enumerate(aligned_overlapping_columns):
        if (i==0): # First column
            new_x = x + int(3*patch_width/2)
            fused_mask[:,x:new_x] = col[:,0:int(3*patch_width/2)]
        elif (i==len(aligned_overlapping_columns)-1): # last column
            new_x = x + int(3*patch_width/2)
            fused_mask[:,x:new_x] = col[:,int(patch_width/2):]
        else: # columns in between
            new_x = x + patch_width
            fused_mask[:,x:new_x] = col[:,int(patch_width/2):int(3*patch_width/2)]
        x = new_x
    
    return fused_mask
        
#%%
def find_edges(mask):
    '''
    This function finds the edges of labeled objects in the mask.
    '''

    padded_mask = np.pad(mask,1,mode='edge')

    center = padded_mask[1:-1,1:-1]
    up = padded_mask[0:-2,1:-1]
    up_left = padded_mask[0:-2,0:-2]
    left = padded_mask[1:-1,0:-2]

    compare = np.array((center!=up,center!=up_left,center!=left))
    edges = np.logical_or.reduce(compare)
    
    return edges

#%%
def split_cells_on_mask(mask):
    '''
    This function separates objects on the mask
    based on edges.
    '''
    edges = find_edges(mask)
    compare = np.array((mask > 0, ~edges))
    segmented_mask = np.logical_and.reduce(compare)
    
    return segmented_mask

#%%
def remove_small_cells(fused_mask, cell_size_threshold):
    '''
    This function removes cells smaller than cell_size_threshold.
    '''
    # Separate cells on fused mask
    separated_cells = split_cells_on_mask(fused_mask)

    # Remove cells smaller than cell_size_threshold
    label_objects, nb_labels = ndi.label(separated_cells)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > cell_size_threshold
    mask_sizes[0] = 0
    filtered_cells = mask_sizes[label_objects]

    # Label again
    filtered_fused_mask, nb_labels = ndi.label(filtered_cells)
    
    return filtered_fused_mask

#%%
def segment_fused_image_with_cellpose(fused, diameter, channels,
                                      edge_thickness=60, similarity_threshold=0.8,
                                      cell_size_threshold=100, patch_height=512, patch_width=512,
                                      num_cpu_cores=None):

    overlap = int(patch_height/2)
    print('>>>> CREATING PATCHES.')
    # Make fused image larger s.t. an integer number of patches fit inside it.
    [M,N,C] = np.shape(fused)
    [new_fused, patch_locations] = enlarge_fused_image(fused, patch_height=patch_height, patch_width=patch_width, overlap=overlap)
    num_patches_nn = np.size(patch_locations[1]) # number of patches on horizontal axis
    
    # Make a list of patches
    patch_list = fused_to_patches(new_fused, patch_locations, patch_height=patch_height, patch_width=patch_width)
    print('>>>> Number of patches to predict: %d'%len(patch_list))
    
    # Predict patches with cellpose
    print('>>>> STARTING CELLPOSE.')
    model = models.Cellpose(gpu=use_GPU, model_type='cyto')
    mask_list, flows, styles, diams = model.eval(patch_list, diameter=diameter, flow_threshold=None, channels=channels)
    
    # Align vertical patches into overlapping columns
    if (num_cpu_cores==None):
        num_cpu_cores = multiprocessing.cpu_count()
    print('>>>> CREATING OVERLAPPING COLUMNS...')
    print('>>>> Using a parallel pool with ', num_cpu_cores, ' CPU workers.')
    overlapping_columns = Parallel(n_jobs=num_cpu_cores)(delayed(create_overlapping_columns_parallel)(nn, mask_list, patch_locations, edge_thickness, similarity_threshold) for nn in range(0,num_patches_nn))
    overlapping_columns,max_value = make_all_cell_values_unique(overlapping_columns)
    
    # Align overlapping columns
    print('>>>> ALIGNING OVERLAPPING COLUMNS...')
    print('>>>> Using a parallel pool with ', num_cpu_cores, ' CPU workers.')
    aligned_overlapping_columns = Parallel(n_jobs=num_cpu_cores)(delayed(align_overlapping_columns_parallel)(nn, overlapping_columns, edge_thickness, similarity_threshold, max_value) for nn in range(0,num_patches_nn-1,2))
    
    print('>>>> CREATING FUSED MASK...')
    # Combine columns into the fused mask
    [new_M,new_N,new_C] = np.shape(new_fused)
    fused_mask = concatenate_overlapping_columns(aligned_overlapping_columns, new_M, new_N, patch_width)
    fused_mask = fused_mask[0:M,0:N]
    
    # Remove cells smaller than cell_size_threshold
    filtered_fused_mask = remove_small_cells(fused_mask, cell_size_threshold)
    
    return filtered_fused_mask

#%%
@jit(nopython=True) # Function is compiled s.t. it runs in machine code
def find_network(segmented,R=4):
        
    # Initialize network matrix:
    M,N = segmented.shape # segmented image has M rows and N columns
    num_cells = np.max(segmented) # number of cells in the segmented image
    network = np.zeros((num_cells, num_cells))
    
    # Loop over all pixels of interest
    for c in np.arange(R,N-R):
        for r in np.arange(R,M-R):
            # Get the neighborhood (square with radius R)
            nbh = segmented[r-R:r+R+1,c-R:c+R+1]
            poi = nbh[R,R] # Pixel of interest (poi) is the center of nbh
            
            if poi > 0:
                # Loop through pixels in neighborhood
                # and connect the values to poi in the network.
                # Numba JIT likes loops!
                for row in nbh:
                    for value in row:
                        if (value!=0 and value!=poi): 
                            network[poi-1, value-1] = 1 # center is connected to the object labeled by value
                            network[value-1, poi-1] = 1 # object labeled by value is connected to center
    return network

#%%
def choose_images_with_gui(initial_dir):

    fused_img_list = []

    # find out whether the user wants one file processed,
    # or multiple files.    
    master = tk.Tk()
    tk.Label(master, text="How many wells do you want to process?").grid(row=0, sticky=tk.W)
    var1 = tk.IntVar()
    tk.Checkbutton(master, text="just one", variable=var1).grid(row=1, sticky=tk.W)
    var2 = tk.IntVar()
    tk.Checkbutton(master, text="multiple (batch processing)", variable=var2).grid(row=2, sticky=tk.W)
    tk.Button(master, text='Continue', command=master.destroy).grid(row=3, sticky=tk.W, pady=4)
    tk.mainloop()
    
    # There are 4 different choices:
    if (var1.get() and var2.get()):              # user made 2 choices, not 1
        raise ValueError('Please choose only one option, not both.')
    
    if (not(var1.get()) and not(var2.get())):    # user made no choice
        raise ValueError('Please choose one of both options.')
    
    if var1.get() and not(var2.get()):           # user chose to process only 1 file
        root = tk.Tk()
        root.filename =  tk.filedialog.askopenfilename(initialdir = initial_dir,title = "Select file",filetypes = [("tif files",".tif")])
        fused_img_list.append(root.filename)
        root.destroy() # close GUI
    
    if not(var1.get()) and var2.get():           # user chose to process multiple files
        root = tk.Tk()
        root.filename = tk.filedialog.askdirectory(initialdir = initial_dir, title = "Select root directory")
        root.destroy() # close GUI
        
        # Specify wells
        well_chooser = tk.Tk()
        well_chooser.geometry("500x200")
        wellstr = tk.StringVar(value='Enter wells do you want to process, separated by commas.')
        entry = tk.Entry(well_chooser, text='Which wells do you want to process?', textvariable=wellstr)
        entry.grid(row=0, sticky=tk.W)
        entry.place(x=10,y=10,width=400)
    
        allwells = tk.IntVar()
        allwellbutton = tk.Checkbutton(well_chooser, text='All wells in root', variable=allwells)
        allwellbutton.grid(row=1, sticky=tk.W)
        allwellbutton.place(x=10,y=50)
        
        button = tk.Button(well_chooser, text='Continue', command=well_chooser.destroy)
        button.grid(row=2, sticky=tk.W, pady=4)
        button.place(x=10,y=90)
        tk.mainloop()
        
        # Make a list of well names
        if allwells.get():  # the user wants to do all wells in the root directory
            well_list = next(os.walk(root.filename))[1]
        else:               # the user has specified the wells
            wells = wellstr.get()
            well_list = wells.split(',')
        
        # Make a list of filenames
        for well in well_list:
            filename = os.path.join(root.filename, well, well+'_fused_RGB.tif')
            
            # check if file exists
            if os.path.isfile(filename):
                fused_img_list.append(filename)
            else:
                raise ValueError('Sorry, the file '+filename+' does not exist.')
                
    return fused_img_list

#%%
# ----------------------------SPECIFY PARAMETERS-------------------------------

# General
[patch_height, patch_width] = [1024,1024]
edge_thickness = 60                     # size of edge region in pixels
similarity_threshold = 0.7              # overlapping cells which overlap more than similarity_threshold will be merged
cell_distance = 8                       # number of pixels that separate two cells (must be 2 or larger)
minimal_cell_size = 500               # minimal cell size in pixels. 100 for HeLa at 10x, 200 for HeLa at 20x

# Cellpose parameters (see below for more info)
channels = [1,3]                        # [1,3] for R=cytoplasm and B=nucleus
cell_diameter = 45                      # 125 for HeLa cells at 20x, 97 for HeLa cells at 10x

# Cell properties to measure
properties = ['label', 'area', 
              'centroid', 'orientation', 
              'minor_axis_length', 
              'major_axis_length', 
              'eccentricity', 'perimeter']

# Machine-specific
initial_dir = r"M:\tnw\bn\dm\Shared"    # Starting directory
num_cpu_cores = None                    # set to 'None' to use the maximum possible number of cores

# CELLPOSE PARAMETERS
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
# channels = [1,3] # IF YOU HAVE R=cytoplasm and B=nucleus

# or if you have different types of channels in each image
# channels = [[2,3], [0,0], [0,0]]

# if diameter is set to None, the size of the cells is estimated on a per image basis
# you can set the average cell `diameter` in pixels yourself (recommended) 
# diameter can be a list or a single number for all images

#%%
# -------------------------------START CODE------------------------------------

# tkinter GUI
fused_img_list = choose_images_with_gui(initial_dir)
    
# loop through files
for filename in fused_img_list:

    directory = os.path.split(filename)[0]
    print(directory)
    
    print('>>>> READING FUSED IMAGE IN ' + directory + '...')
    fused = io.imread(filename)
    segmented = segment_fused_image_with_cellpose(fused, cell_diameter, channels,
                                                  edge_thickness=edge_thickness, similarity_threshold=similarity_threshold,
                                                  cell_size_threshold=minimal_cell_size, patch_height=patch_height, patch_width=patch_width,
                                                  num_cpu_cores=num_cpu_cores)
    
    print('>>>> FINDING NETWORK...')
    M,N = segmented.shape # segmented image has M rows and N columns
    
    # Do cell measurements
    measurements = measure.regionprops_table(segmented, properties=properties)
    
    # Find network
    network = find_network(segmented,R=cell_distance)
    sparse_contact_matrix = sparse.lil_matrix(network)
    
    # Store network and cell measurements in dictionary
    mdic = {'contact_matrix': sparse_contact_matrix,'img_size': [M, N]}
    for key,value in measurements.items():
        key = key.replace('-','') # Fields containing '-' are invalid in Matlab
        mdic[key] = value
    
    print('>>>> SAVING OUTPUT...')
    segmentation_output_path = os.path.join(directory, 'fused_segmentation.tif')
    io.imsave(segmentation_output_path, segmented)
    
    network_output_path = os.path.join(directory, 'network.mat')
    savemat(network_output_path, mdic)  
    print('Graph is saved succesfully as .mat file.')
    print('\n\n')
