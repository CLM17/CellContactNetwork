# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 12:21:05 2020

@author: clmvandenheuve
"""

import matplotlib.pyplot as plt
import pickle

#dict = {'Python' : '.py', 'C++' : '.cpp', 'Java' : '.java'}
#f = open("models/file.pkl","r")
#pickle.dump(dict,f)
#f.close()
output = 'cos_soma_final'
summary_file = 'models/' + output + "_summary.pkl"
with open(summary_file, 'rb') as f:
    summary = pickle.load(f)

plt.figure()
plt.plot(summary['acc'])
plt.title(summary['name']+' training accuracy')
plt.xlabel('epoch')
plt.ylabel('accuracy')

plt.figure()
plt.plot(summary['loss'])
plt.title(summary['name']+' training loss')
plt.xlabel('epoch')
plt.ylabel('loss')