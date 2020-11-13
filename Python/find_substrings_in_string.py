# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 08:35:26 2020

@author: clmvandenheuve
"""


a = [0]
string = 'N:\tnw\BN\dm\Shared\Kasper\Experiments\WKS024\20x\MFGTMP_14_12_01'
start = 0
c = 0

while True:
    c = 0
    end_string = string[start:]
    i = end_string.find("\\")
    if i == -1:
        break
    start = start + i + 1
    a.append(a[c-1] + i)