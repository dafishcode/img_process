#PROCESS
#--------------
#---------------
#===============================================================================
def rotate(Freg, opslist,degree): #rotate all images into correct orientation for registration
#===============================================================================
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import ndimage
    
    rotimglist = list(range(len(opslist)))

    f, axarr = plt.subplots(2,5, sharey = True, sharex = True, figsize = (40,30))
    axarr = axarr.flatten()
    for i in range(len(opslist)):
        ops = np.load(Freg + opslist[i], allow_pickle=True)
        ops = ops[()]
        raw = ops['meanImg']
        rotated_img = ndimage.rotate(raw, degree, reshape=False)
        rotimglist[i] = rotated_img
        axarr[i].matshow(rotated_img)
    
    f.tight_layout()
    plt.show()

    return rotimglist

#===============================================================================
def rotate_point(origin, point, angle): #rotate all images into correct orientation for registration
#===============================================================================
    import math
    import numpy as np
    angle = np.radians(angle)
    dot1 = origin[0] + math.cos(angle) * (point[0] - origin[0]) - math.sin(angle) * (point[1] - origin[1])
    dot2 = origin[1] + math.sin(angle) * (point[0] - origin[0]) + math.cos(angle) * (point[1] - origin[1])
    return(dot1,dot2)

#===============================================================================
def fishspec(Fdata, prefx = ''):
#===============================================================================
    import os
    import re
    
    if prefx: prefx = '^' + prefx + '.*' 
    names = os.listdir(Fdata)
    r     = re.compile(prefx + 'ZFRR.*')
    folds = list(filter(r.match, names))
 
    Zfish = []
    for f in folds:
        cfld = next(os.walk(Fdata + os.sep + f))[1]
        Cond = []
        for c in cfld:
            Tpaths = []
            tifs = os.listdir(Fdata + os.sep + f + os.sep + c)
            r    = re.compile('^.*[tif|tiff|TIF|TIFF]$')
            tifs = list(filter(r.match, tifs))
            Tpaths = []
            for t in tifs:
                Tpaths.append(Fdata + os.sep + f + os.sep + c + os.sep + t)
                
            Cond.append({'Name':c, 
                         'Path':Fdata + os.sep + f + os.sep + c, 
                         'Tifs':tifs,
                         'Tpaths':Tpaths})
            
        Zfish.append({'Cond':Cond, 'Name':f[len(prefx)-2:]})
    
    return Zfish

#===============================================================================
def meancalc(imgs, Fimg, noimages = 100, delfirst = True, crop = False, doplot = True):
#===============================================================================
    import numpy as np
    import ants
    import os

    print('I found ' + str(len(imgs)) + ' images')

    # Load subsection of tifs
    #---------------------------------------------------------------------------
    maxno   = np.min([len(imgs),noimages])
    loadi   = np.linspace(0,len(imgs)-1,maxno)
    loadi   = loadi.astype(int)
    print('Of these I''m loading ' + str(maxno))
    if delfirst: 
        loadi = np.delete(loadi, 0)
        print('I''m ignoring the first volume')
        
    # Load initial image for dimensions
    #---------------------------------------------------------------------------
    if type(imgs[0]) == str:     
        templ = ants.image_read(Fimg + os.sep + imgs[0])
   
    elif type(imgs[0]) == ants.core.ants_image.ANTsImage:                        
        templ = imgs[0]

    if crop:    
        templ = ants.crop_indices(templ, [0,0,1], templ.shape)
    
    mean_arr    = np.multiply(templ.numpy(), 0);
    imglist     = []
    
    for i in loadi:
        
        if type(imgs[0]) == str:     
            img = ants.image_read(Fimg + os.sep + imgs[i])
        elif type(imgs[0]) == ants.core.ants_image.ANTsImage:                        
            img = imgs[i]  
        if crop:    img = ants.crop_indices(img, [0,0,1], img.shape)

        mean_arr    = mean_arr + img.numpy() / maxno
        imglist.append(img)

    mimg = ants.from_numpy(mean_arr)
    if doplot: ants.plot(mimg, axis=2, slices = range(8), figsize=3)
    
    return mimg, imglist
        
#===============================================================================
def pointtrans(Fish, F):
#===============================================================================
    import os
    import ants
    import pandas as pd
    import numpy as np
    
    # Apply registration to the CMN identified cells
    #---------------------------------------------------------------------------
    for c in range(len(Fish["Cond"])):
        print('Transforming points to standard space for condition ' + str(c+1))
        Freg   = F["Freg"] + os.sep + Fish["Name"] + os.sep + Fish["Cond"][c]["Name"]
        Ff2cf  = Freg + os.sep + 'FUN2CF'
        os.listdir(Ff2cf)

        cs = pd.DataFrame(Fish["Cond"][c]["Pixels"])
        cs.columns = ['x', 'y', 'z'] 
        tcs = np.multiply(cs, (.3,.3,15))

        ncs = ants.apply_transforms_to_points(3, tcs, \
                                       [ Ff2cf +os.sep+ 'cf2fun_R.mat', 
                                         Ff2cf +os.sep+ 'cf2fun_S.mat', 
                                         Ff2cf +os.sep+ 'cf2fun_S.nii.gz'],  \
                                       whichtoinvert = [True, True, False])

        tcs = np.multiply(ncs, (1,1,1))
        nncs = ants.apply_transforms_to_points(3, tcs, \
                                        [ F["Ftrans"] +os.sep+ 'ref2cf_R.mat', 
                                          F["Ftrans"] +os.sep+ 'ref2cf_S.mat', 
                                          F["Ftrans"] +os.sep+ 'ref2cf_S.nii.gz'], \
                                        whichtoinvert = [True,True,False]) 
        
        Fish["Cond"][c]["ZBBCoord"] = nncs.values
    
    return Fish


#PLOT
#--------------
#---------------
#===============================================================================
def fishplot(img, overl = '', orient = 'axial', sliceno = 20, al = .5, col = 'magma'):
#===============================================================================
    # N.B This breaks frequently - I have no idea what is wrong with the implementation
    
    import ants
    import numpy as np
    r = img.shape
    
    if   orient == 'coronal':     axs = 0; ri = 0
    elif orient == 'sagittal':    axs = 1; ri = 1
    elif orient == 'axial':       axs = 2; ri = 2 
    
    # ALERT ALERT I HAVE NO IDEA WHAT THE SLICE INDEXING WANTS FROM ME
    # Does it live in the physical space?
    # Does it live in the voxel space?
    # I DON'T KNOW - so I'm fudging it
    
    sliceno = min([sliceno,r[ri]])
    
    if orient == 'axial': mx_dim = r[ri] - 1
    else: mx_dim = r[ri] * img.spacing[ri]-1
    sli     = list(map(int, np.ndarray.tolist(np.linspace(0, mx_dim, sliceno))))
        
    if not overl: ants.plot(img, axis = axs, slices = sli, figsize = 6)
    
    else: ants.plot(img, overl, overlay_cmap = col, overlay_alpha= al, axis = axs, slices = sli, figsize = 6)
        


#SAVE
#--------------
#---------------
#===============================================================================
def savemeanimg(Freg, opslist, Frotate, degree): #save mean image as hyperstack for registration
#===============================================================================
    import os
    from PIL import Image
    omlist = []

    for i in range(len(Frotate)):
        omlist.append(Image.fromarray(Frotate[i]))
    omlist[0].save(Freg + os.sep + opslist[0][:opslist[0].find('run')+6] + "_meanimgstack" + "_" + str(degree) + "deg.tif", save_all=True,
               append_images=omlist[1:])

    
    
    