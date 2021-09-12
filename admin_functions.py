import avalanches as crfn


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
def name_list(path, experiment, num, string): #return name list
#=============================
    import os 
    import glob
    os.chdir(path + 'Project/' + experiment)
    if num < 10:
        out = '0' + str(num)
    elif num >9:
        out = str(num)
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
def list_of_list(rows, cols): #expects a list of lists
#===============================
    listoflist = [[[] for i in range(cols)] for j in range(rows)]
    return(listoflist)


#=======================================================================
def mean_distribution(distlist): #Generate mean distribution 
#=======================================================================
    import numpy as np
    comb_vec = []
    for i in range(len(distlist)):
        comb_vec = np.append(comb_vec, distlist[i])
    av = np.unique(comb_vec, return_counts=True)[0]
    freq = (np.unique(comb_vec, return_counts=True)[1]).astype(int)//len(distlist)
    mean_vec = []
    for e in range(freq.shape[0]):
        mean_vec = np.append(mean_vec, np.full(freq[e],av[e]))
    return(mean_vec)


#PROCESS
#=============================
#=====================================================================
def parallel_func(cores, savepath, iter_list, func, param_list, name, variables, mode): 
#=====================================================================
    """This function allows parallel pooling of processes using functions
    cores = number of cores 
    savepath = path for saving
    iter_list = list with parameter inputs that you will parallel process (inputs must be at start of function)
    func = function name
    param_list = list containing remaining function parameters 
    name = filename for saving, should be unique if mode = save_group
    variables = list containing name endings for each variable, if function returns multiple
    mode = output type:
        save_single - saves each variable of function output individually
        save_group - saves all batched function outputs in a list
        NA - returns all batched function outputs in a list, without saving
    """
    
    from multiprocessing import Pool
    import numpy as np
    pool = Pool(cores) #number of cores
    count = 0

    batch_list = list(range((np.int(len(iter_list)/cores)))) #define number of batches

    for i in range(len(batch_list)): #process each batch
        cores_inputs = list(range(cores)) #define input for each core
        for e in range(len(cores_inputs)):  
            sub_iter_list = iter_list[count:count+1] #Find current iter value - add to subset iter_list
            sub_iter_list.extend(param_list) #Append current iter value onto remaining parameter
            cores_inputs[e] = sub_iter_list 
            count+=1
        batch_list[i] = pool.starmap(func, cores_inputs) #pool process on each core

        if mode == 'save_single':
            for t in range(cores):  #loop through each core in current loop
                for f in range(len(batch_list[i][t])):
                    save_var = batch_list[i][t][f] #function output for current core in current batch
                    save_name = name + '-' + str(cores_inputs[t][0]) + '-' + variables #save name based on iterable parameter
                    np.save(savepath + save_name, save_var)
        
    if mode != 'save_single':
        #Append all calculated value together
        if isinstance(batch_list[0][0], int) or isinstance(batch_list[0][0], float) :
            return_me = np.hstack(np.array(batch_list))
        else:
            return_list = list(range(len(batch_list[0][0])))
            new_array = np.vstack(np.array(batch_list))
            return_me = [new_array[:,i] for i in range(new_array.shape[1])]

        if mode == 'save_group':
            save_name = name
            np.save(savepath + save_name, return_me)

        else:
            return(return_me)
      
    
#=====================================================================
def parallel_class(cores, savepath, iter_list, func, param_list, name, variables, mode): 
#=====================================================================
    """This function allows parallel pooling of processes using classes
    cores = number of cores 
    savepath = path for saving
    iter_list = list with parameter inputs that you will parallel process (inputs must be at start of function)
    func = function name
    param_list = list containing remaining function parameters 
    name = filename for saving, should be unique if mode = save_group
    variables = list containing name endings for each variable, if function returns multiple
    mode = output type:
        save_single - saves each variable of function output individually
        save_group - saves all batched function outputs in a list
        NA - returns all batched function outputs in a list, without saving
    """
    
    from multiprocessing import Pool
    import numpy as np
    pool = Pool(cores) #number of cores
    count = 0

    batch_list = list(range((np.int(len(iter_list)/cores)))) #define number of batches
    for i in range(len(batch_list)): #process each batch
        cores_inputs = list(range(cores)) #define input for each core
        for e in range(len(cores_inputs)):  
            sub_iter_list = iter_list[count:count+1] #Find current iter value - add to subset iter_list
            sub_iter_list.extend(param_list) #Append current iter value onto remaining parameter
            cores_inputs[e] = sub_iter_list 
            count+=1
        batch_list[i] = pool.starmap(func, cores_inputs) #pool process on each core
        
        

        
        if mode == 'save_single':
            import scipy.sparse
            for t in range(cores):  #loop through each core in current loop
                for s in range(len(variables)):
                    save_var = batch_list[i][t].__dict__[variables[s]] #function output for current core in current batch
                    save_name = name + '-' + str(cores_inputs[t][0]) + '-' + variables[s] #save name based on iterable parameter
                    sparse_matrix = scipy.sparse.csc_matrix(save_var)
                    scipy.sparse.save_npz(savepath + save_name, sparse_matrix)
                    #np.save(savepath + save_name, save_var)

    if mode != 'save_single':
    
        #Append all calculated values together
        if len(variables) == 1:
            if isinstance(batch_list[0][0].__dict__[variables[0]], int) or isinstance(batch_list[0][0].__dict__[variables[0]], float):
                count=0
                return_me = list(range(len(iter_list)))
                for first in range(len(batch_list)):
                    for second in range(len(batch_list[0])):
                        return_me[count] = batch_list[first][second].__dict__[variables[0]]
                        count+=1
            else:            
                count=0
                return_me = list_of_list(len(variables),len(iter_list))
                for first in range(len(batch_list)):
                    for second in range(len(batch_list[0])):
                        for third in range(len(variables)):
                            return_me[third][count] = batch_list[first][second].__dict__[variables[third]]
                        count+=1       

        if len(variables) > 1:
            count=0
            return_me = list_of_list(len(variables),len(iter_list))
            for first in range(len(batch_list)):
                for second in range(len(batch_list[0])):
                    for third in range(len(variables)):
                        return_me[third][count] = batch_list[first][second].__dict__[variables[third]]
                    count+=1
        
        if mode == 'save_group':
            save_name = name
            np.save(savepath + save_name, return_me)
        
        else:
            return(return_me)
        

        
#=======================================================================================        
def timeprint(per, r, numrows, name):
#=======================================================================================
    """ Print current time step every percentile
        per = how often you want to print (as percintiles)
        r = current iterator value
        numrows = total number of steps
        name = name to output
    """
    if r % round((per*numrows/100)) == 0: 
            print("Doing number " + str(r) + " of " + str(numrows) + " for " + name)
            
            
#MATHS
#=============================
#=============================
#=======================================================================================
def window(size, times): #make window of given size that is divisible of time series
#=======================================================================================
    """Returns the window size that is the closest divisor of a timeseries to given input
    Inputs:
    size - ideal window size
    times - overall trace shape
    
    Returns: 
    size - window size that is divisible by trace (rounds up)
    n_windows - number of windows that split up trace
    """
    for i in range(times):
        if times % size ==0:
            break
        else:
            size+=1
    n_windows = int(times/size)
    return(size, n_windows)

#=======================================================================================
def ttest(mydf, label, variable, comp_list, mode):
#=======================================================================================
    from scipy import stats 
    #Single comparison - label to compare to first element in list
    if mode == 'single':
        vals = list_of_list(len(comp_list)-1, 5)
        sig = 0.05/(len(comp_list)-1)
        base = comp_list[0]
        for i in range(len(comp_list)-1):
            vals[i][0], vals[i][1] = stats.ttest_rel(mydf[variable].where(mydf[label] == base).dropna(),mydf[variable].where(mydf[label] == comp_list[i+1]).dropna())[0],stats.ttest_rel(mydf[variable].where(mydf[label] == base).dropna(),mydf[variable].where(mydf[label] == comp_list[i+1]).dropna())[1]
            vals[i][2] = sig
            vals[i][4] = str(base) + ' - ' + str(comp_list[i+1])
            if vals[i][1] < sig:
                vals[i][3] = 'Significant'
            else:
                vals[i][3] = 'Not significant'
    
    
    if mode == 'multiple':
        vals = list(range(len(comp_list)))
        ncomp = 0
        for i in range(len(comp_list)):
            ncomp+= (len(comp_list)-1) - i
        sig = 0.05/ncomp
        
        for i in range(len(comp_list)):
            subval = list_of_list(len(comp_list), 5)
            for e in range(len(comp_list)):
                subval[e][0], subval[e][1] = stats.ttest_rel(mydf[variable].where(mydf[label] == comp_list[i]).dropna(),mydf[variable].where(mydf[label] == comp_list[e]).dropna())[0],stats.ttest_rel(mydf[variable].where(mydf[label] == comp_list[i]).dropna(),mydf[variable].where(mydf[label] == comp_list[e]).dropna())[1]
                subval[e][2] = sig
                subval[e][4] = str(comp_list[i]) + ' - ' + str(comp_list[e])
                if subval[e][1] < sig:
                    subval[e][3] = 'Significant'
                else:
                    subval[e][3] = 'Not significant'
            vals[i] = subval
    
    
    return(vals)

