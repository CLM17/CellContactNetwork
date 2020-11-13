# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:43:27 2020

@author: LocalAdmin
"""

from timeit import default_timer as timer
import matplotlib.pyplot as plt
import networkx as nx
from Helpers.network_creator_helpers import create_HELA_network

start = timer()

well_list = ['C02','C03','C04','C05','C06','C07',
             'D02','D03','D04','D05','D06','D07']

root = 'Experiments/WKS024'

for well in well_list:
    
    com, contact_matrix, G, on_edge = create_HELA_network(root, well)
    
    print('Total runtime is ', timer()-start)
    
    pos = nx.get_node_attributes(G,'pos')
    ax = plt.figure(figsize=(10,10))
    nx.draw(G,pos,node_size=.5)



