# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:46:11 2020

@author: LocalAdmin
"""

import numpy as np
from skimage import io
from skimage.morphology import binary_closing, binary_opening, binary_dilation, binary_erosion
from skimage.draw import line, circle, circle_perimeter
from skimage.segmentation import watershed

from scipy import ndimage
from scipy import sparse
from scipy.io import savemat
import os
from numba import jit
import pandas as pd
import networkx as nx

from skimage import img_as_ubyte, img_as_uint
from skimage.filters import gaussian
from skimage.feature import peak_local_max

#%%
def get_well_locations(folder, well):
    '''
    This function extracts the locations of the well from the excel file.
    '''
    # Import Well locations.xlsx as pandas dataframe and extract the locations
    well_locs_df = pd.read_excel( os.path.join(folder, 'Well locations.xlsx'), index_col=0 )
    xc_well = int(well_locs_df['xc'][well])
    yc_well = int(well_locs_df['yc'][well])
    diameter_well = int(well_locs_df['diameter'][well])
    
    # Raise error if the information is not complete:
    if np.any( np.isnan([xc_well, yc_well, diameter_well]) ):
        raise NameError('The well locations are incomplete. Please update the excel file "Well locations.xlsx"')
    
    return xc_well, yc_well, diameter_well

#%%
def apply_threshold(ch, th_method, well_folder, well):
    '''
    This function applies a threshold on an input image. The thresholding 
    method is specified by the method parameter.
    
    If "unet" is the specified method, the function will use the prediction 
    of the unet convolutional neural network to threshold on. Only works if 
    there is a prediction (wellname_prediction.png) in the well folder.

    Parameters
    ----------
    ch : 2D numpy array
        Normalized image channel.
    method : string
        Specification of thresholding method. Possibilities are:
        otsu, multiotsu, li, minimum, mean, isodata, unet.

    Returns
    -------
    2D boolean array
        Thresholded image.

    '''
    
    if th_method == 'otsu':
        from skimage.filters.thresholding import threshold_otsu
        threshold = threshold_otsu(ch)
    elif th_method == 'li':
        from skimage.filters.thresholding import threshold_li
        threshold = threshold_li(ch)
    elif th_method == 'minimum':
        from skimage.filters.thresholding import threshold_minimum
        threshold = threshold_minimum(ch)
    elif th_method == 'mean':
        from skimage.filters.thresholding import threshold_mean
        threshold = threshold_mean(ch)
    elif th_method == 'isodata':
        from skimage.filters.thresholding import threshold_isodata
        threshold = threshold_isodata(ch)
    elif th_method == 'unet':
        try:
            ch = io.imread(os.path.join(well_folder, well+'_prediction.png'))
            ch = ch / np.max(ch)
            threshold = 0.5
        except:
            raise ValueError('You need to make a U-Net prediction and save it as "wellName_predict" before you can threshold with a unet.')
    else:
        raise ValueError('Invalid thresholding method.')
    

    return ch > threshold
    
    
#%%
def remove_outside_well(img_list, xmin, ymin, xc_well, yc_well, diameter_well):
    '''
    Function excludes all pixels of the image outside the circle with center (xc_well, yc_well) 
    and diameter diameter_well.

    Parameters
    ----------
    img_list : list of 2D numpy arrays.
        list of 2D images of which you want to remove the outside of the well.
    xc_well : int
        x-coordinate of the well center.
    yc_well : int
        y-coordinate of the well center.
    diameter_well : int
        diameter of the well center.

    Returns
    -------
    new_img_list : list of 2D numpy arrays.
        list of images where the outside of the well is removed.

    '''
    
    new_img_list = []
    
    for img in img_list:
        # Create new black image with the same size:
        new_img = np.zeros(img.shape)
        # Fill the new image with the input image inside the well
        rr, cc = circle(yc_well - ymin, xc_well - xmin, int(diameter_well / 2))
        new_img[rr,cc] = img[rr,cc]
        new_img_list.append(new_img)
    
    return new_img_list

#%%
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
def find_nuclei_com(bin_DAPI, N):
    '''
    Function finds the center of mass (COM) of nuclei.
    Input: img = a grayscale image of DAPI channel, N = size of structuring element used for opening the image.
    Output: com = array with centers of mass, lbl = labels of objects.
    '''

    # Open the image with a N x N structuring element and label result
    opened = binary_opening(bin_DAPI, selem=np.ones([N,N]))
    lbl, n_obj = ndimage.label(opened)

    # find center of mass of labels
    com_list = ndimage.measurements.center_of_mass(opened, lbl, np.arange(1,n_obj+1))
    
    # store the coms in a numpy array (is compatible with numba)
    com = np.zeros((len(com_list), 3))
    for i in range(len(com)):
        com[i,:] = [ int(i+1), com_list[i][0], com_list[i][1] ]
    
    return com, lbl

#%%
def distance_transform_watershed(th, connectivity=np.ones((3,3)), sigma=1, min_distance=10):
    '''
    This function segments all cells in a mask with a distance transform watershed.
    Th is a thresholded (binary) image.
    Sigma is the sigma used for smoothing the distance transform.
    Min_distance is the minimal distance between two cell centers in pixels.
    '''
    dist = ndimage.distance_transform_edt(th)
    dist_smooth = gaussian(dist, sigma=1)

    local_maxi = peak_local_max(dist_smooth, threshold_abs=0, min_distance=min_distance, indices=False)

    markers = ndimage.label(local_maxi)[0]
    
    connectivity = np.ones((3,3))
    labels = watershed(dist_smooth, markers, connectivity=connectivity, mask=th, watershed_line=True)
    
    return labels, markers

#%%
def make_markers(com, N, M):
    '''
    Function loops over the locations of the centers of mass of nuclei.
    It outputs an NxM image markers with white pixels on the COM locations.
    '''
    markers = np.zeros((N, M))
    for i,x,y in com:
        markers[int(x)-1, int(y)-1] = i
    
    return markers

#%%
def seeded_watershed(bin_MEM, markers):
    '''
    Function labels individual cells in MEMBRITE channel with the seeded watershed algorithm.
    It first performs a histogram equalization and minimum thresholding.
    
    Input: 
        img: 2D image of membrite channel.
        makers: image of same dimensions as img with seeds.
    Output: labeled image.
    '''
    wts = watershed(bin_MEM, markers=markers, mask=bin_MEM, watershed_line=False)
    
    return wts

#%%
@jit(nopython=True) # compile to C++ for faster performance
def find_network(wts, connectivity=8):
    '''
    Function finds which cells are connected to each other.
    
    Input:
        com = a list of tuples indicating the locations of the centers of mass of nuclei (row, colum).
        wts = watershed-labeled MEMBRITE image.
    
    Output:
        network = a dictionary with the cell number as key and a list as value.
        The first item in value is a tuple indicating the center of mass of the cell's nucleus (row, column).
        The other items in value are the cell numbers of the cells connected to the key cell.  
    '''
    # Define structuring element:
    selem = np.ones((2,2))
    if connectivity == 4:
        selem[0,0] = 0
    
    # Initialize network dict with COM locations:
    N,M = wts.shape
    num_cells = np.max(wts)
    network = np.zeros((num_cells, num_cells))
        
    # Loop through the pixels in the image:
    for y in np.arange(1,N):
        for x in np.arange(1,M):
            
            # Get neighborhood (2x2 structuring element) and center pixel
            nbh = wts[y-1:y+1, x-1:x+1]
            prod = selem * nbh
            center = int(prod[1,1])

            if center > 0: # if the center pixel belongs to an object

                north = int(prod[0,1]) # pixel above center pixel
                west = int(prod[1,0])  # pixel left of center pixel
                northwest = int(prod[0,0]) # pixel left above center

                if north > 0 and center != north:
                    network[center-1, north-1] = 1 # center is connected to the object above it
                    network[north-1, center-1] = 1

                if west > 0 and center != west:
                    network[center-1, west-1] = 1 # center is connected to the object left of it
                    network[west-1, center-1] = 1
                    
                if northwest > 0 and center != northwest:
                    network[center-1, northwest-1] = 1 # center is connected to the object left above of it
                    network[northwest-1, center-1] = 1
    
    return network

#%%
def matrix_to_graph(matrix, com, N):
    '''
    Function converts the network matrix to a dictionary
    Input:
        matrix: a sparse matrix indicating the network
        com: numpy array with locations of centers of mass
        N: height of image
    Output: dictionary with network
    '''
    G = nx.Graph()
    
    for row in com:
        i = int(row[0])
        com_pos = tuple([row[2], N-row[1]]) # convert r,c to x,y
        G.add_node(i, pos=com_pos)
    
    i_array, j_array = matrix.nonzero()
    for i,j in zip(i_array, j_array):
        if i==j:
            matrix[i,j] = 0
            continue
        if (i+1,j+1) not in G.edges:
            G.add_edge(i+1, j+1)
                
    return G

#%%
def list_nodes_on_edge(wts, xc_well, yc_well, diameter_well):
    '''
    This function takes as input a network dictionary and the corresponding
    watershed-segmented image.
    It then makes a list of all cells that are on the edge.
    
    Input: 
        G = network graph
        wts = watershedded image.
        xc_well, yc_well, diameter_well: specifications of well position
        
    Output:
        on_edge: list of cells that touch the edge.
    '''
    
    rr_on_edge, cc_on_edge = circle_perimeter(yc_well, xc_well, int(diameter_well / 2))
    on_edge = []
    
    # Loop over the edges of the image, and keep track
    # of objects on the edge: 
    for r, c in zip(rr_on_edge, cc_on_edge):
        if wts[r,c] > 0 and wts[r,c] not in on_edge:
            on_edge.append(wts[r,c])
    
    return on_edge

#%%
def draw_network(G, img_DAPI, img_MEMBRITE):
    '''
    Function draws the network as nodes and edges on the image.
    '''
    # Create RGB output image:
    N,M = img_DAPI.shape
    output = np.zeros((N,M,3))
    output[:,:,0] = img_MEMBRITE / np.max(img_MEMBRITE) # membrite is RED channel
    output[:,:,2] = img_DAPI / np.max(img_DAPI)         # DAPI is BLUE channel
    
    radius = 5 # size of nodes (DAPI COM)
    pos = nx.get_node_attributes(G,'pos')

    # Loop over nodes:
    for node, loc in pos.items():
        x,y = loc
        c = int(x)
        r = N - int(y)
        # get indices of the circle
        rr, cc = circle(r, c, radius)
        # make sure the indices are not larger than the image edge
        rr = np.minimum(N-1, rr)
        cc = np.minimum(M-1, cc)
        
        # Draw node as white circles:
        output[rr, cc, :] = [1, 1, 1]

    # Loop over edges:
    for edge in G.edges:
        i,j = edge
        xi, yi = pos[i]  # x,y position node i
        xj, yj = pos[j]  # x,y position node j
        ci = int(xi)     # column position of node i
        cj = int(xj)     # column position of node j
        ri = N - int(yi) # row position of node i
        rj = N - int(yj) # row position of node j
        rr, cc = line(int(ri),int(ci), int(rj),int(cj))

        # Draw edge as while line:
        output[rr, cc, :] = [1, 1, 1]

    return output

#%%
def create_HELA_network(folder, well, th_method_MEM):
    '''
    This function finds the network of the well.
    '''
   
    well_folder = os.path.join(folder, 'well '+well)
    full_img = io.imread( os.path.join( well_folder, well+'_fused.tif') )

    print('I read the image. \n Busy...')

    ch_DAPI = full_img[:,:,2] / np.max(full_img[:,:,2])
    ch_MEM = full_img[:,:,0] / np.max(full_img[:,:,0])
    
    N_selem = 5 # size of structuring element for opening + closing
    # original image shape
    N0, M0 = ch_DAPI.shape
    
    ymin, xmin = [0,0]
    # Crop image in right shape (for unet)
    if th_method_MEM == 'unet':
        w = 512 # tile width
        ol = 100 # tile overlap
        ymin, ymax, ynumtiles = crop_img(N0, w, ol)
        xmin, xmax, xnumtiles = crop_img(M0, w, ol)
        ch_DAPI = ch_DAPI[ymin:ymax, xmin:xmax]
        ch_MEM = ch_MEM[ymin:ymax, xmin:xmax]
    print('cropped')
    
    # new image shape
    N, M = ch_DAPI.shape
    
    # Apply Otsu thresholding on DAPI channel
    th_DAPI = apply_threshold(ch_DAPI, 'otsu', well_folder, well)
    # apply multiotsu thresholding on MEM channel
    th_MEM = apply_threshold(ch_MEM, th_method_MEM, well_folder, well)
    print('threshold applied')
    
    # remove the outside of the well
    xc_well, yc_well, diameter_well = get_well_locations(folder, well)
    th_DAPI, th_MEM = remove_outside_well([th_DAPI, th_MEM], xmin, ymin, xc_well, yc_well, diameter_well)
    print('outside of well is removed')

    # Find the DAPI centers of mass and the labeled image
    com, lbl_DAPI = find_nuclei_com(th_DAPI, N_selem)
    print('DAPI COM found')
    # Make markers for seeded watershed
    markers = make_markers(com, N, M)
    print('Seeded watershed markers are made')
    # Perfrom seeded watershed
    wts = seeded_watershed(th_MEM, markers)
    print('Seeded watershed is done')
    # Find the network with connected component labeling
    contact_matrix = find_network(wts)
    print('Connected component labeling is done')
    # Convert matrix to graph and create the output img
    G = matrix_to_graph(sparse.lil_matrix(contact_matrix), com, N)
    on_edge = list_nodes_on_edge(wts, xc_well, yc_well, diameter_well)
    print('Matrix converted to graph')
    output_img = draw_network(G, ch_DAPI, ch_MEM)
    print('Network is drawn')
    io.imsave( os.path.join(well_folder,well+'_network.tif'), img_as_ubyte(output_img) )
    io.imsave( os.path.join(well_folder,well+'_watershed.tif'), img_as_uint(wts) )
    print('Watershed and network images are saved succesfully.')
    
    sparse_contact_matrix = sparse.lil_matrix(contact_matrix)
    mdic = {'contact_matrix': sparse_contact_matrix,
            'nuclei_com': com,
            'img_size': [N, M],
            'on_edge': on_edge}

    savemat(os.path.join(well_folder,'graph_'+well+'.mat'), mdic)
    print('Graph is saved succesfully to .mat file.')
    
    return com, contact_matrix, G, on_edge