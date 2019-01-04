#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 20:49:21 2018

@author: dominicburrows
"""

import os
import glob 
import numpy as np

#paths
datapath = "/Users/dominicburrows/Documents/PhD/Imaging/181203_grin_Gc6s_SA/F1_suite2p"
os.chdir(datapath)   #set current directory to datapath
planes = sorted(glob.glob("plane*"))   #list all planes in each folder


  
for i in range(len(planes)):            #loop through each plane
    os.chdir(datapath + "/" + planes[i]) #set new directory as each plane folder in turn 
    x = np.load("iscell.npy")            #new variable which loads iscell.npy (contains prob of being a cell) file into array
    F = np.load("F.npy") [x [:,1] > 0.5, :] #same as above for F.npy, but calls only the row numbers which in x are above 0.5
    stats = np.load("stat.npy") [x [:,1] > 0.5]  #same as above for stat.npy
    centers = np.zeros((len(stats),2))     #create matrix of zerox with shape size of F/stats arrau
    for  j in range (len(stats)):          #loop through length of stats
        centers [j,] = stats [j] ['med']   #fills centers matrix with info about cell number (j) and corresponding med values for each cell
    
    
        
    
    

            

            
            
    
              
        
