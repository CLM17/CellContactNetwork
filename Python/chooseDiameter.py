# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 13:13:30 2021

@author: lukasvandenheu
"""


#%%
import numpy as np
from skimage import io
from skimage import measure
from skimage import color
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
import random
# from timeit import default_timer as timer
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('TkAgg') 

import tkinter as tk
import tkinter.filedialog
#from tkinter import Label, W, IntVar, Button, Checkbutton, mainloop

from cellpose import models
from cellpose import utils
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
def choose_fused_image_with_gui(initial_dir):
    '''
    This function starts up a tkinter GUI which asks the user
    to open a tiff image.
    It outputs the path to the fused image (as a string).
    '''
    root = tk.Tk()
    root.filename =  tk.filedialog.askopenfilename(initialdir = initial_dir,title = "Select fused RGB imge",filetypes = [("tif files",".tif")])
    path_to_fused_image = root.filename
    root.destroy() # close GUI
    
    return path_to_fused_image

#%%
def get_display_size():
    '''
    This function outputs the height and width of the computer display in pixels.
    '''
    root = tkinter.Tk()
    root.update_idletasks()
    root.attributes('-fullscreen', True)
    root.state('iconic')
    height = root.winfo_screenheight()
    width = root.winfo_screenwidth()
    root.destroy()
    return height, width

#%%
def update_list(list_, new_entry):
    '''
    This function removes the first entry from the beginning of a list, 
    and adds new_entry at the end.
    Example: if list_ = [1,2,3] and new_entry = 4, then this function outputs
    list_ = [2,3,4].
    '''
    list_ = list_[1:]        # Remove first entry from list
    list_.append(new_entry)  # Add new entry to list
    return list_

#%%
def choose_random_patch_nr(patch_list):
    '''
    This function takes as input a list of patches. It outputs a random nr
    (between 0 and length(patch_list)) which is the index of some patch with a
    maximum pixel intensity larger than 10.
    '''
    num_patches = len(patch_list)
    nr = np.random.randint(0,num_patches)
    while np.max(patch_list[nr]) < 10:  # max intensity of patch must be larger than pixel value 10
        nr = np.random.randint(0,num_patches)
    return nr

#%%
def init_mask_list(patch_list,nr,model,diameter_list,channels):
    '''
    This function segments the same patch with different diameters.
    
    INPUTS
    ------
        patch_list (list)
        A list of all available patches.
        
        nr (int)
        Index of patch of interest.
        
        model (Cellpose model object)
        Defines the Cellpose model used for segmentation.
        
        diameter_list (list)
        List of different diameters to test.
        
        channels (list)
        Specifies channels for Cellpose model
    
    OUTPUT
    ------
        masks (list)
        List with the same length as diameter_list, with the segmentations of 
        the patch at different diameters.
    '''
    masks = []
    for diam in diameter_list:
        m, flows, styles, diams = model.eval(patch_list[nr], diameter=diam, flow_threshold=None, channels=channels)
        masks.append(m)
    
    return masks

#%%
def show_segmentations(patch_list,model,diameter,channels,fig_width,fig_height,outline_color):
    '''
    This function shows segmented masks.
    
    LOCAL INPUTS
    
    '''
    
    global diameter_list
    global masks
    
    model = models.Cellpose(gpu=use_GPU, model_type='cyto')
    patch = patch_list[nr]
    
    if (diameter != diameter_list[-1]):
        diameter_list = update_list(diameter_list, diameter)
        new_mask, flows, styles, diams = model.eval(patch, diameter=diameter, flow_threshold=None, channels=channels)
        masks = update_list(masks, new_mask)
        
    plt.close('all')
    plt.figure(figsize=(fig_width,fig_height))
    count = 0
    for mask,diam in zip(masks,diameter_list):
        count = count + 1
        plt.subplot(2,num_images,count)
        plot_outlines(mask, patch, outline_color)
        plt.title('diam = '+str(diam))
        
        plt.subplot(2,num_images,count+num_images)
        plt.imshow(mask)
    plt.show()

#%%
def change_patch_nr(patch_list,model,diameter,diameter_list,channels):
    '''
    
    '''
    global nr
    global masks
    nr = choose_random_patch_nr(patch_list)
    masks = init_mask_list(patch_list,nr,model,diameter_list,channels)
    show_segmentations(patch_list,model,diameter,channels,fig_width,fig_height,outline_color)

#%%
def plot_outlines(mask, img, rgb_color):
    outlines = utils.masks_to_outlines(mask)
    outX, outY = np.nonzero(outlines)
    imgout= img.copy()
    imgout[outX, outY] = np.array(rgb_color)
    plt.imshow(imgout)

#%%
    
[patch_height, patch_width] = [1024,1024]
num_images = 3
channels = [1,3]
outline_color = [0,255,0]

tk_width = 200
tk_height = 200
dpi = 100
screen_height,screen_width = get_display_size()

fig_width = int( ( screen_width - 3*tk_width/2 ) / dpi)
fig_height = int( screen_height / dpi )

initial_dir = r"M:\tnw\bn\dm\Shared"    # Starting directory
path_to_fused = choose_fused_image_with_gui(initial_dir)

print('>>>> READING FUSED IMAGE...')
fused = io.imread(path_to_fused)
[M,N,C] = np.shape(fused)

overlap = int(patch_height/2)

print('>>>> CREATING PATCHES.')
# Make fused image larger s.t. an integer number of patches fit inside it.
[new_fused, patch_locations] = enlarge_fused_image(fused, patch_height=patch_height, patch_width=patch_width, overlap=overlap)
patch_list = fused_to_patches(new_fused, patch_locations, patch_height=patch_height, patch_width=patch_width)

# Init
diameter_list = [50,50,75]
diameter = 100  # New 'dummy' diameter to add to the list
nr = choose_random_patch_nr(patch_list)
model = models.Cellpose(gpu=use_GPU, model_type='cyto')  
masks = init_mask_list(patch_list,nr,model,diameter_list,channels)  

#%%
show_segmentations(patch_list,model,diameter,channels,fig_width,fig_height,outline_color)

gui = tk.Tk()
x_left = int(screen_width - 3*tk_width/2)
y_top = int(screen_height - 3*tk_height/2)
gui.geometry("{}x{}+{}+{}".format(tk_width,tk_height,x_left,y_top))

diam_str = tk.IntVar(value='130')
entry = tk.Entry(gui, text='Diameter:', textvariable=diam_str)
entry.grid(row=0, sticky=tk.W)
entry.place(x=55,y=10,width=90)

B1 = tk.Button(gui, text = "Update", command = lambda: show_segmentations(patch_list,model,int(diam_str.get()),channels,fig_width,fig_height,outline_color))
B1.place(x=55,y=90,width=90)

B2 = tk.Button(gui, text = "Show different patch", command = lambda: change_patch_nr(patch_list,model,int(diam_str.get()),diameter_list,channels))
B2.place(x=30,y=150,width=140)

gui.mainloop()
