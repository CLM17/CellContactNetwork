# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 14:31:56 2021

@author: lukasvandenheu
"""

from skimage import io
import numpy as np
import os
from skimage import measure
from scipy.io import savemat

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
def get_mean_intensities(stack, segmented):
    '''
    This function outputs the mean intensity values of cells
    on the segmented image at each timepoint.
    '''
    
    # Get an array with all cell values on segmented image
    cell_values = np.unique(segmented)
    cell_values = np.delete(cell_values, 0)
    
    # Number of cells and number of timepoints
    num_cells = len(cell_values)
    num_timepoints = np.shape(stack)[0]
    
    # Init
    mean_intensities = np.zeros((num_cells, num_timepoints))
    std_intensities = np.zeros((num_cells, num_timepoints))

    # Loop over cells and get the intensity values of pixels inside
    # the cell region over time.
    # Then, calculate the mean and std of the intensity values at
    # each timepoint.
    for ii,value in enumerate(cell_values):
        intensity_values = stack[:,segmented==value]
        for t in range(0,num_timepoints):
            mean_intensities[ii,t] = np.mean(intensity_values[t])
            std_intensities[ii,t] = np.std(intensity_values[t])
                               
    return mean_intensities, std_intensities

#%%
def get_baseline_intensity(stack):
    '''
    This function calculates the average intensity of the whole
    image at all timepoints.
    '''
    num_timepoints = np.shape(stack)[0]
    baseline_intensity = np.zeros((num_timepoints))
    std_baseline = np.zeros((num_timepoints))

    for t in range(0,num_timepoints):
        baseline_intensity[t] = np.mean(stack[t,:,:])
        std_baseline[t] = np.std(stack[:,:,t])
    
    return baseline_intensity, std_baseline

#%%
def make_time_axis(num_timepoints, sample_period):
    
    total_time = num_timepoints * sample_period
    return np.linspace(0, total_time, num_timepoints)

#%%
def timelapse_to_intensity_levels(stack, segmented, sample_period, properties):
    '''
    This function takes as input an image timelapse (stack)
    and a segmented image with the regions of interest (segmented).
    '''
    
    num_timepoints = np.shape(stack)[0]
    time_axis = make_time_axis(num_timepoints, sample_period)
    mean_intensities, std_intensities = get_mean_intensities(stack, segmented)
    baseline_intensity, std_baseline = get_baseline_intensity(stack)
    measurements = measure.regionprops_table(segmented, properties=properties)

    results = {}
    results['time_axis'] = time_axis
    results['mean_intensities'] = mean_intensities
    results['std_intensities'] = std_intensities
    results['baseline_intensity'] = baseline_intensity
    results['std_baseline'] = std_baseline

    for key,value in measurements.items():
        key = key.replace('-','') # Fields containing '-' are invalid in Matlab
        results[key] = value

    roi = split_cells_on_mask(segmented)
    
    return results, roi

#%%
# ----------------------------SPECIFY PARAMETERS-------------------------------
sample_period = 2      # take 1 image every # seconds
properties = ['label', 'area', 'centroid', 'orientation', 'minor_axis_length', 'major_axis_length', 
              'eccentricity', 'perimeter']

root = r'M:\tnw\bn\dm\Shared\Lukas\Experiments\CLM001\HCA'
well  = 'B03'

#%%
# -------------------------------START CODE------------------------------------

well_folder = os.path.join(root, well)
path_to_stack = os.path.join(well_folder, 'B03_stack.tif')
path_to_segmented = os.path.join(well_folder, 'B03_stack_mask.png')

# Read input images
print('>>>> READING INPUT IMAGE STACK...')
stack = io.imread(path_to_stack)
segmented = io.imread(path_to_segmented)

print('>>>> CALCULATING MEAN INTENSITY LEVELS OF ALL CELLS OVER TIME...')
# Convert timelapse to intensity levels
results, roi = timelapse_to_intensity_levels(stack, segmented, sample_period, properties)

print('SAVING OUTPUT...')
# Save results
roi_path = os.path.join(well_folder, well+'_ROI.tif')
results_path = os.path.join(well_folder, well+'_intensity_levels.mat')
savemat(results_path, results)  
io.imsave(roi_path, roi)
print('Succesfully saved a ROI image (%s_ROI.tif) and intensity levels (%s_intensity_levels.mat) in %s.'%(well,well,well_folder))
