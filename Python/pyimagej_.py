# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 09:14:33 2020

@author: clmvandenheuve
"""

import os

os.environ["JAVA_HOME"] = "C:/Users/clmvandenheuve/Anaconda3/envs/ijenv/Library"
os.environ["PYJNIUS_JAR"] = "C:/Users/clmvandenheuve/Anaconda3/envs/ijenv/share/pyjnius/pyjnius.jar"
#os.environ["PATH"] = os.environ["PATH"] +":/anaconda3/opt/maven/bin"

os.getenv("JAVA_HOME")
os.getenv("PYJNIUS_JAR")
os.getenv("PATH")

import imagej