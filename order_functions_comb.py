import numpy as np
from skimage import io
from skimage.io import imread
#from skimage.external.tifffile import imread, TiffFile
import os
import glob


#paths
datapath = "/Volumes/Dominic 10tb/Data2process/Project/PTZ-WILDTYPE/201201-FOXG1+-/"
num = 'F02'
new_name = 'BLN-PTZ05-PTZ20-WILDTYPE'
os.chdir(datapath)   #set current directory to datapath
tif_list = sorted(glob.glob(num  + "*tif"))  


def tiff_split(data_path, tif_list, cut):   #define function tiff_split, and with two inputs
    curr_tif = tif_list[0]
    new_dir = num + '-' + new_name + '-' + curr_tif[curr_tif.find('2p'):curr_tif.find('run')+6]
    os.mkdir(new_dir)         #create directory named with first two letters of file in tif
    counter = 1                   #keeps running count of file number for saving/not saving
    names = 1                     #keeps running tab of file number for naming
    for e in range(len(tif_list)):
        tif = tif_list[e]
        img = imread(tif)             #function which reads the tif variable as a numpy array

        for i in range(img.shape [0]): #i loop through number of frames in each tif file
            if counter%cut == 0:         #if counter is divisible by six, do not save
                counter +=1            #counter iterate when it is at 6
            else:
                io.imsave(arr =img [i,:,:],  fname= data_path + new_dir + "/" + "im_" +  str(names).zfill(6)  + ".tif")
                counter+= 1 #else if not divisible by 6,all other values save the ith image in array and name
                names+=1    #iterate counter and also name, but not name if i=6 so we get consecutive named files


tiff_split(datapath, tif_list, cut= 11)  