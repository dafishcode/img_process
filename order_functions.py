#imports

import numpy as np
from skimage import io
from skimage.io import imread
#from skimage.external.tifffile import imread, TiffFile
import os
import glob



#paths

def tiff_split(data_path, tif):   #define function tiff_split, and with two inputs
    os.mkdir(tif[tif.find('F'):tif.find('run')+6])         #create directory named with first two letters of file in tif
    img = imread(tif)             #function which reads the tif variable as a numpy array
    counter = 1                   #keeps running count of file number for saving/not saving
    names = 1                     #keeps running tab of file number for naming
   
    for i in range(img.shape [0]): #i loop through number of frames in each tif file
        if counter%11 == 0:         #if counter is divisible by six, do not save
            counter +=1            #counter iterate when it is at 6
        else:
            io.imsave(arr =img [i,:,:],  fname= data_path + tif[tif.find('F'):tif.find('run')+6] + "/" + "im_" +  str(names).zfill(6)  + ".tif")
            counter+= 1 #else if not divisible by 6,all other values save the ith image in array and name
            names+=1    #iterate counter and also name, but not name if i=6 so we get consecutive named files
       
       #i counts from 1 to max frames, and then is plugged in to io.imsave to
       #ensure that the ith frame is then saved as an array with correct name


datapath = "/Volumes/Dominic 10tb/Data2process/Project/PTZCON-WILDTYPE/200317-WILDTYPE/"
os.chdir(datapath)   #set current directory to datapath
tifs = sorted(glob.glob("*tif"))  
for i in range(len(tifs)):    #loop through number of files in tifs folder
    tiff_split(data_path=datapath, tif = tifs [i])  
   
     
datapath = "/Volumes/Dominic 10tb/Data2process/Project/PTZCON-WILDTYPE/210125-WILDTYPE/"
os.chdir(datapath)   #set current directory to datapath
tifs = sorted(glob.glob("*tif"))  
for i in range(len(tifs)):    #loop through number of files in tifs folder
    tiff_split(data_path=datapath, tif = tifs [i])  
   

   
     
