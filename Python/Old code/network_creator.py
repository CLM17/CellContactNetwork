# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:43:27 2020

@author: LocalAdmin
"""

from timeit import default_timer as timer
import matplotlib.pyplot as plt
import networkx as nx
from network_creator_helpers import create_HELA_network

start = timer()

well = 'B03'
folder = 'Data/WKS023/2020-09-09'
th_method_MEM = 'unet'           # thresholding method for membrite channel.

com, contact_matrix, G, on_edge = create_HELA_network(folder, well, th_method_MEM)

print('Total runtime is ', timer()-start)

pos = nx.get_node_attributes(G,'pos')
ax = plt.figure(figsize=(30,30))
nx.draw(G,pos,node_size=.5)
