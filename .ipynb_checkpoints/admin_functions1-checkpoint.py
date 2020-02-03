#SORT
#=============================
#=============================
import avalanches as crfn


#=============================
def name_zero(path, experiment, start, end): #add 0 to number if above 10
#=============================
    import os 
    count = 0
    listme = list(range(start, end+1))
    os.chdir(path + 'Project/' + experiment)
    for i in range(start, end+1):
        if i < 10:
            num = '0' + str(i)
        elif i >9:
            num = str(i)
        listme[count] = num
        count+=1
    return(listme)

#=============================
def name_list(path, experiment, inp, string): #return name list
#=============================
    import os 
    import glob
    os.chdir(path + 'Project/' + experiment)
    if inp < 10:
        out = '0' + str(inp)
    elif inp >9:
        out = str(inp)
    return(sorted(glob.glob('*E-' + str(out) + string)))

#=============================
def name_template(name, mode): #return name list
#=============================
    if mode == 'short':
        temp = name[:name.find('run')+6] 
    
    if mode == 'long':
        temp = name[:name.find('.npy')-3] 
        
    return(temp)

#==============================
def save_name(i, name_li): #find save name
#===============================
    return(name_li[i][:name_li[i].find('run')+6])

#PROCESS
#=============================
#=============================
#=====================================================================
def parallel(cores, datalist, func, paramlist): #make sure n_cores is divisible by total number
#=====================================================================
    from multiprocessing import Pool
    import numpy as np
    pool = Pool(cores)
    count = 0
    for i in range((np.int(len(datalist)/cores))):
        paramlist_levels = list(range(cores))
        for e in range(len(paramlist_levels)):
            newlist = datalist[count:count+1]
            newlist.extend(paramlist)
            paramlist_levels[e] = newlist
            count+=1
        output = pool.starmap(func, paramlist_levels)

#MATHS
#=============================
#=============================

