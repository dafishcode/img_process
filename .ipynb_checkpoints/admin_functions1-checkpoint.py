

#SORT
#=============================
#=============================
#=============================
def name_zero(pad, start, end, step): #add 0 to string of number - for saving 
#=============================
    import os 
    import numpy as np
    
    if pad == 'front': 
        count = 0
        listme = list(range(start, end+1, step))
        for i in range(start, end+1, step):
            if i < 10:
                num = '0' + str(i)
            elif i >9:
                num = str(i)
            listme[count] = num
            count+=1
        return(listme)

    if pad == 'back': 
        count, count1 = 0,0
        looplist = np.arange(start, end + step, step)
        listme = list(range(0, looplist.shape[0]))
        lenlist = list(range(looplist.shape[0]))
        for i in looplist:
            lenlist[count1] = len(str(round(i, len(str(step)))))
            count1 +=1
        for i in looplist:
            if len(str(round(i,len(str(step))))) < np.max(lenlist):
                num = str(round(i,len(str(step)))) + '0'
            else:
                num = str(round(i,len(str(step))))
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
def name_template(namelist, mode): #return name list
#=============================
    if mode == 'short':
        temp = namelist[0][:namelist[0].find('run')+6] 
    
    if mode == 'long':
        temp = namelist[0][:namelist[0].find('.npy')-3] 
        
    if mode == 'param':
        temp = namelist[0][:namelist[0].find('run')+12] + 'bin'  + namelist[1][namelist[1].find('run')+7:namelist[1].find('run')+14]
        
    return(temp)

#=============================
def repeat_list(name, length): #make list of same name repeated for given length
#=============================    
    itlist = list(range(length))
    for i in range(len(itlist)):
        itlist[i] = name
    return(itlist)


#==============================
def save_name(i, name_li): #find save name
#===============================
    return(name_li[i][:name_li[i].find('run')+6])

#==============================
def list_of_list(lengths): #expects a list of lists
#===============================
    listoflist = []
    for i in range(len(lengths)):
        listoflist = [listoflist for i in range(lengths[i] + 1)]
    return(listoflist)


#PROCESS
#=============================
#=============================
#=====================================================================
def parallel(cores, listlist, func, paramlist): #make sure n_cores is divisible by total number
#=====================================================================
    from multiprocessing import Pool
    import numpy as np
    pool = Pool(cores)
    count = 0

    if len(listlist) == 1:
        for i in range((np.int(len(listlist[0])/cores))):
            paramlist_levels = list(range(cores))
            for e in range(len(paramlist_levels)):
                newlist = listlist[0][count:count+1]
                newlist.extend(paramlist)
                paramlist_levels[e] = newlist
                count+=1
            output = pool.starmap(func, paramlist_levels)
            
    if len(listlist) > 1:
        for i in range(len(listlist)):
            if len(listlist[i]) != len(listlist[0]):
                print('Input lists must be the same length')
                exit()
            
        for i in range((np.int(len(listlist[0])/cores))):
            paramlist_levels = list(range(cores))
            for e in range(len(paramlist_levels)):
                for t in range(len(listlist)):
                    if t == len(listlist)-1:
                        newlist.extend(paramlist)
                        break
                    newlist = listlist[t][count:count+1]
                    newlist.extend(listlist[t+1][count:count+1])
                    paramlist_levels[e] = newlist
                    count+=1 
            output = pool.starmap(func, paramlist_levels)

#MATHS
#=============================
#=============================

