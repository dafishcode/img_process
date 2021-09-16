import admin_functions as adfn
import ants


#PROCESS
#--------------
#---------------
#===============================================================================
def rotate(Freg, ops_list,degree): 
#===============================================================================
    """
    This function rotates mean images in a stack to the defined orientation for registration and plots the outcome. 
    
    Inputs:
    Freg (str): Path to mean image stack
    ops_list (list): list of ops files for each plane, each containing the mean image stack for a plane
    degree (int): angle of rotation
    
    Returns:
    rotimg_list (list): list containing rotated mean image for each plane
    
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import ndimage
    
    rotimg_list = list(range(len(ops_list)))

    f, axarr = plt.subplots(2,5, sharey = True, sharex = True, figsize = (40,30))
    axarr = axarr.flatten()
    for i in range(len(ops_list)): #Loop through each plane
        ops = np.load(Freg + ops_list[i], allow_pickle=True)
        ops = ops[()]
        raw = ops['meanImg'] #Load mean image
        rotated_img = ndimage.rotate(raw, degree, reshape=False) #Apply rotation
        rotimg_list[i] = rotated_img #Add each rotated image to list
        axarr[i].matshow(rotated_img)
    
    f.tight_layout()
    plt.show()

    return rotimg_list


#===============================================================================
def load_prereg(experiment, file, modality, condition, prefix, disk, Fdata):
#===============================================================================
    """
    This function loads prereg data in ants and also loads trace and coord data, and plots output.
    
    Inputs:
    experiment (str): experiment subfolder name
    file (str): fish name
    modality (str): imaging modality
    condition (str): imaging condition
    prefix (str): genotype, if None no genotype, if 'Yes' genotype is added
    disk (str): path to mean images
    Fdata (str): path to trace and coords data
    
    Returns:
    f (ants image): image stack to register
    Fish (dict): dictionary containing cell and trace data
    meanimg_fold (list): list containing the folder with the image tif
    """
    
    import numpy as np
    import os

    #Find fish image
    #--------------------------------------------------------
    img_type= adfn.return_files(disk + 'Project/', file + os.sep + modality + os.sep, '**')[0]
    meanimg_fold = adfn.return_files(disk + 'Project/', file + os.sep + modality + os.sep + img_type , '*' + condition + '*meanimgstack*')
    num = meanimg_fold[0][meanimg_fold[0].find(modality)-3:meanimg_fold[0].find(modality)-1]
    dpf = meanimg_fold[0][meanimg_fold[0].find('dpf')-1:meanimg_fold[0].find('dpf')+3]
    
    if prefix == 'Yes':
        pref = file[file.find(modality)-4:file[0].find(modality)-2]
        
    else:
        pref = ''

    #Find coords and trace files
    #--------------------------------------------------------
    name_list = ['*_realcoord.npy', '*_realtrace.npy', '*_binarised.npy', '*_deltaff.npy']
    dict_list = ['coord', 'trace', 'bind', 'dff']
    Fish = {}
    def grab(Fdata, experiment, condition, pref, number, age, name):
        out = adfn.return_files(Fdata, experiment, '*' + experiment + pref + '-' + number + '*' + age + '*' + condition + name )
        if len(out) > 1:
            print('More than one file found - check list length')
            return()

        if len(out) < 1:
            print('No files found found - check list length')
            return()
        return(out)
            
    for i, nom in enumerate(name_list):
        Fish.update({dict_list[i]: np.load(grab(Fdata, experiment, condition, pref, num, dpf, nom)[0])})
        
    # Visualise fish
    #---------------------------------------------------------------------------
    f = ants.image_read(disk + 'Project/' + file + os.sep + modality +os.sep + img_type + os.sep + meanimg_fold[0])
    f.set_spacing([.3, .3, 15.])
    f.set_direction(np.array([[-1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])) 
    ants.plot(f, axis = 2, nslices = 5, slices = np.arange(0, 5, 1), figsize = 10)
    
    savename = adfn.save_name(grab(Fdata, experiment, condition, pref, num, dpf, nom)[0])
    
    return(f, Fish, meanimg_fold, savename)



#===============================================================================
def rotate_coords(xyz, meanimg_list):
#===============================================================================
    """
    This applies the rotation already applied to the image, to the cell coordinates.
    
    Inputs:
    xyz (np array): X x Y x Z coordinates of cells
    meanimg_list (list): list of mean images
    
    Returns:
    rot_stack (np array): X x Y x Z coordinates of rotated cells
    to_del (np array): 1d vector of cell indeces that are outside image
    """
    
    import numpy as np
    
    def rotate_point(origin, point, angle): 
        import math
        angle = np.radians(angle)
        x = origin[0] + math.cos(angle) * (point[0] - origin[0]) - math.sin(angle) * (point[1] - origin[1])
        y = origin[1] + math.sin(angle) * (point[0] - origin[0]) + math.cos(angle) * (point[1] - origin[1])
        return(x,y)
    
    origin = 512/2,512/2
    angle = int(meanimg_list[0][meanimg_list[0].find('stack')+6:meanimg_list[0].find('deg')]) #Find angle of rotation in file name
    rot = np.zeros((xyz.shape[0], 2))
    for i in range(xyz.shape[0]):
        rot[i] = rotate_point(origin,xyz[i,:2], angle) #Apply rotation to points

    #Remove any neurons that fall outside of image due to squaring procedure
    rot_stack = np.column_stack((rot, xyz[:,2]))
    rot_stack = rot_stack.astype(int)
    mult = rot_stack[:,0]*rot_stack[:,1]
    zer = np.where(mult <= 0 )
    out_x = np.where(rot_stack[:,0] > 511)
    out_y = np.where(rot_stack[:,1] > 511)
    comb = []
    comb = np.append(zer, np.append(out_x, out_y)) #Neurons outside of image
    rot_stack[comb] = 0
    return(rot_stack, comb)


        #Build coordinate image from pixel coordinates in ants space

#===============================================================================
def lab_Img(ant_img, coords, to_del):
#===============================================================================
    """
    This function builds a coordinate image from the pixel coordinates in ants space. 
    
    Inputs:
    ant_img (ants image): pre-registered, rotated ants image
    coords (np array): X x Y x Z coordinates of rotated cells
    to_del (np array): 1d vector of cell indeces that are outside image
    
    Returns:
    lab_img (ants image): ants image with each pixel labelled by cell number
    """
    
    import numpy as np

    #Create labelled image from preregistered image
    lab_img = ant_img*0
    
    #This will create an image of the xy coordinates in the same space as the input ants image
    count=1 #start counting at 1 as all empty pixels = 0
    for i in range(coords.shape[0]):
        if sum(i == to_del) == 0: #Skip values that were removed in last step
            lab_img[coords[i][1], coords[i][0], coords[i][2]] = int(count) #Antsimage and suite2p coordinates are transpose of eachother - Transpose to match
        count+=1
        
    return(lab_img)




#===============================================================================
def match_pix2cells(pre_coord, warp_coord_img):
#===============================================================================
    """
    This function maps new pixel coordinates onto original cells. 
    
    Inputs:
    pre_coord (np array):  X x Y x Z coordinates of rotated cells pre registration
    warp_coord_img (ants image): ants image with each pixel labelled by cell number
    
    Returns:
    reg_coord (ants image):  X x Y x Z coordinates of cells post registration,
    """


    import numpy as np

    reg_coord = np.zeros((pre_coord.shape[0],3))
    warp_points_arr = np.array(warp_coord_img[:])
    xyz_clust = np.where(warp_points_arr > 0)

    count = 0
    for i in range(pre_coord.shape[0]):
        curr_ind = np.where(warp_points_arr[xyz_clust] == count+1)[0]
        if len(curr_ind) > 0:

            curr_x = xyz_clust[1][curr_ind] #untranspose - swap x and y back into coordinate space
            curr_y = xyz_clust[0][curr_ind]
            curr_z = xyz_clust[2][curr_ind]
            mean_xyz = np.mean(curr_x), np.mean(curr_y), np.mean(curr_z) #Calculate new coordinate location as pixel centroid
            reg_coord[count] = mean_xyz
        count+=1
    return(reg_coord)


#===============================================================================
def reg_label(fixed, moving, label, atlaslab, coord, trace, dff, bind, meanimglist, reg_type, mode):
#===============================================================================
    """
    This function builds a coordinate image from the pixel coordinates in ants space. 
    
    Inputs:
    ant_img (ants image): pre-registered, rotated ants image
    coords (np array): X x Y x Z coordinates of rotated cells
    to_del (np array): 1d vector of cell indeces that are outside image
    
    Returns:
    lab_img (ants image): ants image with each pixel labelled by cell number
    """


    import numpy as np
    import matplotlib.pyplot as plt
    import random
    from matplotlib import cm
    from matplotlib.axes._axes import _log as matplotlib_axes_logger
    matplotlib_axes_logger.setLevel('ERROR')

    if coord.shape[0] != trace.shape[0] or coord.shape[0] != dff.shape[0]  or coord.shape[0] != bind.shape[0]:
        print('Input files not the same shape - check you are loading the correct fish files')
        return()
    else:
    
        #Perform registration of fish image to atlas image
        warp_img = ants.registration(fixed, moving, type_of_transform = reg_type) 
        
        #Inspect registration
        fishplot(fixed,warp_img['warpedmovout'], orient = 'axial', al = 0.7)
        
        #Rotate coordinates
        rot_coords, to_del = rotate_coords(coord, meanimglist)
    
        #Build coordinate image from pixel coordinates in ants space
        coord_img = lab_Img(moving,rot_coords, to_del)
            
        #Apply transformation to coordinate image
        warp_coord_img = ants.apply_transforms(fixed, coord_img, warp_img['fwdtransforms'], interpolator='nearestNeighbor')
        
        #Map warped pixel coordinates to old pixel coordinates
        reg_coord = match_pix2cells(rot_coords, warp_coord_img)
        
        #Create final coordinate and trace arrays - remove cells outside of image
        loc = np.where(reg_coord[:,0] == 0)[0]
        fin_coord = np.delete(reg_coord, loc, 0)
        fin_trace = np.delete(trace, loc, 0)
        fin_dff = np.delete(dff, loc, 0)
        fin_bind = np.delete(bind, loc, 0)
        
        #label coordinates
        new_x = fin_coord[:,0]//2
        new_y = fin_coord[:,1]//2
        new_z = fin_coord[:,2]//2
        lab_xyz = (np.column_stack((new_x, new_y, new_z))).astype(int)
        coarse_reg = list(range(lab_xyz.shape[0]))
        gran_reg = list(range(lab_xyz.shape[0]))
        
        #Label each cell according to atlas labels
        for i in range(len(gran_reg)):
            curr_val = int(label[lab_xyz[i][1],lab_xyz[i][0],lab_xyz[i][2]]) #Transpose to match with label ants img
            gran_reg[i] = np.array(atlaslab[1])[curr_val]
            coarse_reg[i] = np.array(atlaslab[2])[curr_val]

        lab_coord = np.column_stack((fin_coord, gran_reg, coarse_reg))
        
        #Visualise ouputs
        if mode == 'check':
            print(meanimglist[0])
            print('Rotate suite2p coords onto pre-reg fish brain')
            plane_num= 6
            plt.figure(figsize = (10,10))
            plt.matshow(moving[:,:,plane_num])
            ploc = np.where(rot_coords[:,2] == plane_num)
            plt.scatter(rot_coords[:,0][ploc],rot_coords[:,1][ploc], s = 0.8, c = 'red', alpha = 1)
            plt.show()
            
            print(meanimglist[0])
            print('Build coordinate image and stack in z over pre-reg fish brain')
            pre_stackplot = np.zeros((coord_img[:,:,0].shape))
            for i in range(10):
                pre_stackplot +=coord_img[:,:,i]

            #Check that labelled image has been built correctly
            fig,axarr = plt.subplots(figsize = (5,5))
            axarr.matshow(moving[:,:,plane_num])
            axarr.scatter(np.where(pre_stackplot > 0)[1],np.where(pre_stackplot > 0)[0], s = 0.1, c = 'red') #Antsimage and suite2p coordinates are transpose of eachother - Transpose to match
            plt.show()

            
            print(meanimglist[0])
            print('Plot newly warped cell coordinates over atlas (left) and warped fish image (right)')
        
            #Check that all neurons are overlaid correctly over brain - postregistration
            xnum = 180
            curr_warped = warp_img['warpedmovout']
            fig,axarr = plt.subplots(1,2, figsize = (10,10))
            axarr[0].matshow(fixed[:,:,xnum])
            axarr[0].scatter(reg_coord[:,0],reg_coord[:,1], s = 3, alpha = 0.2, c = 'red')
            axarr[1].matshow(curr_warped[:,:,xnum])
            axarr[1].scatter(reg_coord[:,0],reg_coord[:,1], s = 2, alpha = 0.1, c = 'red')
            plt.show()
            
            print(meanimglist[0])
            print('Plot newly warped cell coordinates over atlas/warped fish brain, by individual planes')
            #Check neurons plane by plane
            curr_warped = warp_img['warpedmovout']
            curr_points = reg_coord.astype(int)
            znumlist = np.arange(100, 200, 20)
            xnumlist = np.arange(200, 300,20)
            for num in range(len(znumlist)):
                fig,axarr = plt.subplots(1,4, figsize = (15,15))
                axarr[0].matshow(fixed[:,:,znumlist[num]])
                axarr[0].scatter(curr_points[:,0][curr_points[:,2] == znumlist[num]],curr_points[:,1][curr_points[:,2] == znumlist[num]], s = 2, c = 'red')
                axarr[1].matshow(curr_warped[:,:,znumlist[num]])
                axarr[1].scatter(curr_points[:,0][curr_points[:,2] == znumlist[num]],curr_points[:,1][curr_points[:,2] == znumlist[num]] ,s = 1, c = 'red')
                axarr[2].matshow(fixed[:,xnumlist[num],:])
                axarr[2].scatter(curr_points[:,2][curr_points[:,0] == xnumlist[num]],curr_points[:,1][curr_points[:,0] == xnumlist[num]],s = 3, c = 'red')
                axarr[3].matshow(curr_warped[:,xnumlist[num],:])
                axarr[3].scatter(curr_points[:,2][curr_points[:,0] == xnumlist[num]],curr_points[:,1][curr_points[:,0] == xnumlist[num]], s = 2, c = 'red')
                plt.show()
                
                
            print(meanimglist[0])
            print('Plot same cells for pre-reg (left) and post-reg (right) coordinates - are cell positions retained?')  
            #Check that cell ids are correctly retained

            old_points = rot_coords#coord

            n_cells = 10
            xnum = 150

            fig,axarr = plt.subplots(1,2,figsize = (10,7))
            axarr[0].scatter(old_points[:,0],old_points[:,1], s = 2, c = 'k', alpha = 0.1)
            axarr[1].matshow(fixed[:,:,xnum])
            axarr[1].scatter(reg_coord[:,0],reg_coord[:,1], s = 2, alpha = 0.1, c = 'k')
            for i in range(n_cells):
                colors  = cm.Spectral(np.linspace(0,1,n_cells))
                choose = random.randint(1,len(old_points)+1)
                axarr[0].scatter(old_points[:,0][choose],old_points[:,1][choose], s = 40, c = colors[i], alpha = 1)
                axarr[1].scatter(reg_coord[:,0][choose],reg_coord[:,1][choose], s = 20, c = colors[i], alpha = 1)
            plt.show()

            

            print(meanimglist[0])
            print('Plot cells from the pre-reg matrix and traces (green, left) and post-reg matrix and traces (orange, right) - are traces mapped onto correct cells?')  
            old_points = coord
            n_cells = 5
            for i in range(n_cells):
                choose = random.randint(1,len(old_points)+1)
                if sum(choose == loc) == 0:
                    to_min = sum(choose >= loc)
                    old_val = choose
                    new_val = old_val - to_min
                    old_trace = trace[old_val]
                    new_trace = fin_trace[new_val]
                    print(old_val)
                    fig,axarr = plt.subplots(figsize = (8,1))
                    plt.plot(old_trace, c = 'green')
                    plt.show()
                    print(new_val)
                    fig,axarr = plt.subplots(figsize = (8,1))
                    plt.plot(new_trace, c = 'orangered')
                    plt.show()


                    fig,axarr = plt.subplots(1,2,figsize = (10,7))
                    axarr[0].scatter(old_points[:,0],old_points[:,1], s = 2, c = 'k', alpha = 0.1)
                    axarr[1].matshow(fixed[:,:,xnum])
                    axarr[1].scatter(fin_coord[:,0],fin_coord[:,1], s = 2, alpha = 0.1, c = 'k')
                    axarr[0].scatter(old_points[:,0][old_val],old_points[:,1][old_val], s = 40, c = 'green', alpha = 1)
                    axarr[1].scatter(fin_coord[:,0][new_val],fin_coord[:,1][new_val], s = 20, c = 'orangered',alpha = 1)
                    plt.show()

                
            #Check that all neurons are overlaid correctly over brain - postregistration
            xnum = 150
            curr_warped = warp_img['warpedmovout']
            fig,axarr = plt.subplots(1,2, figsize = (10,10))
            axarr[0].matshow(label[:,:,np.int(xnum/2)])
            axarr[0].scatter(lab_xyz[:,0],lab_xyz[:,1], s = 2, alpha = 0.2, c = 'red')
            axarr[1].matshow(label[:,:,np.int(xnum/2)])
            axarr[1].scatter(lab_xyz[:,0][np.array(coarse_reg) == 'Telencephalon'],lab_xyz[:,1][np.array(coarse_reg) == 'Telencephalon'], s = 2, alpha = 1, c = 'orange')
            axarr[1].scatter(lab_xyz[:,0][np.array(coarse_reg) == 'Diencephalon'],lab_xyz[:,1][np.array(coarse_reg) == 'Diencephalon'], s = 2, alpha = 1, c = 'green')
            axarr[1].scatter(lab_xyz[:,0][np.array(coarse_reg) == 'Midbrain'],lab_xyz[:,1][np.array(coarse_reg) == 'Midbrain'], s = 2, alpha = 1, c = 'cyan')
            axarr[1].scatter(lab_xyz[:,0][np.array(coarse_reg) == 'Hindbrain'],lab_xyz[:,1][np.array(coarse_reg) == 'Hindbrain'], s = 2, alpha = 1, c = 'violet')
            axarr[1].scatter(lab_xyz[:,0][np.array(coarse_reg) == 'nan'],lab_xyz[:,1][np.array(coarse_reg) == 'nan'], s = 2, alpha = 1, c = 'red')
            axarr[1].scatter(lab_xyz[:,0][np.array(coarse_reg) == 'Peripheral'],lab_xyz[:,1][np.array(coarse_reg) == 'Peripheral'], s = 2, alpha = 1, c = 'black')
            axarr[1].scatter(lab_xyz[:,0][np.array(coarse_reg) == 'Unspecified'],lab_xyz[:,1][np.array(coarse_reg) == 'Unspecified'], s = 2, alpha = 1, c = 'black')
            plt.show()


            #Check that all neurons are overlaid correctly over brain - postregistration
            xnum = 200
            curr_warped = warp_img['warpedmovout']
            fig,axarr = plt.subplots(1,2, figsize = (10,10))

            axarr[0].matshow(label[:,np.int(xnum/2),:])
            axarr[0].scatter(lab_xyz[:,2],lab_xyz[:,1], s = 2, alpha = 0.2, c = 'red')
            axarr[1].matshow(label[:,np.int(xnum/2),:])
            axarr[1].scatter(lab_xyz[:,2][np.array(coarse_reg) == 'Telencephalon'],lab_xyz[:,1][np.array(coarse_reg) == 'Telencephalon'], s = 2, alpha = 1, c = 'orange')
            axarr[1].scatter(lab_xyz[:,2][np.array(coarse_reg) == 'Diencephalon'],lab_xyz[:,1][np.array(coarse_reg) == 'Diencephalon'], s = 2, alpha = 1, c = 'green')
            axarr[1].scatter(lab_xyz[:,2][np.array(coarse_reg) == 'Midbrain'],lab_xyz[:,1][np.array(coarse_reg) == 'Midbrain'], s = 2, alpha = 1, c = 'cyan')
            axarr[1].scatter(lab_xyz[:,2][np.array(coarse_reg) == 'Hindbrain'],lab_xyz[:,1][np.array(coarse_reg) == 'Hindbrain'], s = 2, alpha = 1, c = 'violet')
            plt.show()
                
                
                
        return(fin_coord, lab_coord, fin_trace, fin_bind, fin_dff)

    
    


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

    """
    This function plots output of registration. 
    
    Inputs:
    im (ants image): Atlas to plot
    overl (ants image): warped image to overlay onto atlas
    orient (str): plotting orientation
    sliceno (int): number of slices to visualise
    al (int): alpha
    col (str): color map

    """
    
    import ants
    import numpy as np
    r = img.shape
    
    if   orient == 'coronal':     axs = 0; ri = 0
    elif orient == 'sagittal':    axs = 1; ri = 1
    elif orient == 'axial':       axs = 2; ri = 2 
    
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
def savemeanimg(Freg, ops_list, Frotate, degree): #save mean image as hyperstack for registration
#===============================================================================
    """
    Save rotated mean image as hyperstack for registration. 
    
    Inputs:
    Freg (str): Path to mean image stack
    ops_list (list): list of ops files for each plane, each containing the mean image stack for a plane
    Frotate (list): list containing rotated mean image for each plane
    degree (int): angle of rotation
    
    Saves:
    om_list (tif): a hyperstack of all images stacked together as a tif
    
    """

    import os
    from PIL import Image
    om_list = []

    for i in range(len(Frotate)):
        om_list.append(Image.fromarray(Frotate[i]))
    om_list[0].save(Freg + os.sep + ops_list[0][:ops_list[0].find('run')+6] + "_meanimgstack" + "_" + str(degree) + "deg.tif", save_all=True,
               append_images=om_list[1:])

    
    
    