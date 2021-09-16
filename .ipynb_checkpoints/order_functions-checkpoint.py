#=======================================================================
def tiff_split(Fraw, n): # tiff_split takes two inputs
#=======================================================================    
# This function looks i


    import numpy as np
    from skimage.external.tifffile import imread, TiffFile 
    import os
    import glob 
    from skimage import io
    
    os.chdir(Fraw) 
    tif = sorted(glob.glob("*tif")) 
    print('Found ' + str(len(tif)) + ' tiffs')
    
        
    for i in range(len(tif)):    #loop through number of files in tifs folder
    
        # Create directory named with first two letters of file in tif
        # imread function reads tif
        #--------------------------------------------------------------
        os.mkdir(tif [i][0:2])          
        img = imread(tif[i])            
        counter = 1                 
        names = 1
    
        #Loop through all frames in tif, if divisible by nof frames +1 will not save - else save
        #i counts from 1 to max frames, and then is plugged in to io.imsave to
        #ensure that the ith frame is then saved as an array with correct name
        #---------------------------------------------------------------------------------------------
        for j in range(img.shape [0]): 
            if counter%n+1 == 0:        
                counter +=1          
            else:
                io.imsave(arr =img [j,:,:],  fname= Fraw[40:50] + tif[i][0:2] + "/" + "im_" +  str(names).zfill(6)  + ".tif")
                counter+= 1 
                names+=1    
        
            if counter((1000*img.shape[0]/10000)) == 0: print("Doing row" + str(j) + "of" +str(img.shape[0]))
       
 

    