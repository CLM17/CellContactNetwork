# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:39:43 2020

@author: LocalAdmin
"""
from prepare_stitch_helpers import prepare_dataset_for_stitching

target_folder = 'Data\\WKS024\\2020-10-05'
well_rows = ['B']
well_cols = [2,3,4,5,6,7]
w = 8 # number of stitched images

prepare_dataset_for_stitching(target_folder, well_rows, well_cols, w)