import seaborn as sns
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
    

#=================================================           
def boxplot(Fdrop, experiment, plotlist, title):
#====================================================
    blnlist = list(range(np.int(len(plotlist)/3)))
    p5list = list(range(np.int(len(plotlist)/3)))
    p20list = list(range(np.int(len(plotlist)/3)))
    count = 0
    for i in range(len(plotlist)):
        if 'BLN' in plotlist[i]:
            blnlist[count] = np.load(Fdrop +'/Project/' + experiment + os.sep + plotlist[i])
        if 'PTZ05' in plotlist[i]:
            p5list[count] = np.load(Fdrop + '/Project/' + experiment + os.sep + plotlist[i])
        if 'PTZ20' in plotlist[i]:
            p20list[count] = np.load(Fdrop + '/Project/' + experiment + os.sep + plotlist[i])
            count+=1
            
    df = pd.DataFrame({'Baseline': blnlist, 'PTZ 5mM': p5list, 'PTZ 20mM': p20list})
    df
    sns.set_style("darkgrid")

    fig, ax = plt.subplots(figsize = (8,8))

    ax = sns.stripplot(data=df, jitter=True, color="0", size = 8)
    ax = sns.pointplot(data=df, linestyles = ['--'], color="#bb3f3f", size = 10, capsize = .15)
    plt.ylabel(title, size = 30, color = 'black')
    plt.xticks(np.arange(3), ('Baseline','PTZ 5mM', 'PTZ 20mM'), size = 15, color = 'black')
    plt.yticks(color = 'black')
    plt.show()
    sns.reset_orig()
    plt.style.use('dark_background')
    return(df)

#=======================================================================
def rasplot(namelist): # Select which fish data to visualise
#=======================================================================
    from matplotlib import pyplot as plt   
    nplot = len(namelist)
    if nplot == 1:
        fig, ax = plt.subplots(figsize= (15,15))
        plt.title(namelist[0][namelist[0].find('dpf')+4:namelist[0].find('run')-1], size = 30)
        ax.matshow(np.load(namelist[0]), cmap = 'tab20')
        plt.show()

    if nplot > 1: 
        f, axarr = plt.subplots(1,nplot,sharey=True, sharex=True, figsize = (15,15))
        f.subplots_adjust(hspace=0)

        for i in range(len(namelist)):
            axarr[i].set_title(namelist[i][namelist[i].find('dpf')+4:namelist[i].find('run')-1], size = 15)
            axarr[i].matshow(np.load(namelist[i]), cmap = 'tab20')
        plt.show()


#===============================================================================
def univars(mtype, d):
#===============================================================================
    # This function applies a certain calculation on data array(s) stored in the
    # list d. The function types can be:
    # 'p_firing' - expects the peak array as input
      
    from scipy import optimize
    import numpy as np
    
    if mtype == 'p_firing': return (np.sum(d[0], axis = 1) / d[0].shape[1]) 

#===============================================================================
def winslide(trace, binarise, win = 60 * 4, stp = 1, mtype = 'p_firing'):
#===============================================================================
    # This function estimates different univariate measures using a sliding 
    # window approach

    dat = trace
    pks = binarise
    
    starts = np.arange(0, dat.shape[1] - win - 1, win*stp)
    
    # N.B. The numpy stacking is the most confusing thing in the universe,
    # I have no idea why and when I am transposing anything here, so ERRORs
    # are likely 
    
    starts = np.arange(0, dat.shape[1] - win - 1, win*stp)
    wd = np.array([])
    for s in range(len(starts)-1):
        print('Working on time step ' + str(s+1) + ' of ' + str(len(starts)))
        d = dat[:,starts[s]:starts[s]+win]
        p = pks[:,starts[s]:starts[s]+win]
        if mtype == 'p_firing':  td = np.transpose(univars(mtype, [p]))
        wd = np.vstack((wd, td)) if wd.size else td
            
            
    wd = np.transpose(wd)
    
    return wd


#===============================================================================
def fishdot(fish, coord, cols, ax = None, cmap = 'Spectral', al = 0.8):
#===============================================================================
    # This function takes a fish (single fish, single condition) and a single 
    # vector of the same length as numbers of cells and returns a mapping of 
    # that vector onto the cells as a plot
    if ax == None: ax = plt.gca()
        
    outplot = ax.scatter(coord[:,0], coord[:,1], 100, c = cols, 
                         cmap = cmap, alpha = al)
  
    return outplot