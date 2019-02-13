#=======================================================================
def fish_load(Fs2p, Fdrop): # Load imaging datasets 
#=======================================================================
# This function looks in the Fdata folder for suite2p plane files and saves extracted cell traces/coordinates into Fdrop (a numpy file containing a ncells x ntimepoints array of data - from all active cells as defined in suite2p, ordered by plane)


    import os
    import glob 
    import re
    import numpy as np
    import pandas as pd

    dirlist = os.listdir(Fs2p)

    # Find planes of suite2p output
    #------------------------------
    
    r       = re.compile('^plane[0-9].*')
    planelist = list(filter(r.match, dirlist))
    planelist.sort()
    
    print('Found ' + str(len(planelist)) + ' planes') 
    
    # Compile coordinates and trace files for all planes into lists
    #---------------------------------------------------------------------
    
    coord = list((range(len(planelist))))
    signal = list((range(len(planelist))))
    cells = list((range(len(planelist))))
    
    
    for i in range(len(planelist)):
        os.chdir(Fs2p + os.sep + "plane" + str(i)) 
        allcells = np.load("iscell.npy")  
        fl = np.load("F.npy") [allcells [:,1] > 0.5, :] 
        stats = np.load("stat.npy") [allcells [:,1] > 0.5]  
        xy = np.zeros((len(stats),2))  
        
        for  j in range (len(stats)):          
            xy [j,] = stats [j] ['med']
    
        xyz = np.concatenate([xy, np.full((len(fl), 1), i)], axis = 1)
        coord[i] = xyz
        signal[i] = fl
    
    # Concatenate separate arrays for coordinate and xy file into two arrays 
    #--------------------------------------------------------------------

    com_coord = np.concatenate([coord[i] for i in range(len(planelist))])
    com_signal = np.concatenate([signal[i] for i in range(len(planelist))]) 
    
    print('Found ' + str(com_coord.shape[0]) + ' cells')
    
    # Save as three separate files (int file is for R) 
    #-------------------------------------------------

    os.chdir(Fdrop)
    p_id = Fs2p[40:49] + '_' + Fs2p[Fs2p.find('F'):Fs2p.find('F') + 5]
    np.save(p_id + '_'  'com_coord.npy', com_coord)
    np.save(p_id + '_' + 'com_signal.npy', com_signal)
    np.save(p_id + '_' + 'int_com_signal.npy', com_signal.astype("int"))
    
    print('Saved trace and coordinates in dom_data')
    
    
    
    

#========================================================================
def fish_filter(dat, highcut = 500 , lowcut = 0, fun='lin'):   # Filter out frequencies
#========================================================================
# This function removes high frequency components from the trace as noise
# It can also remove low frequency components - as a bandpass filter
# This is necessary for the max/min filtering step

    import numpy as np
    from scipy.fftpack import rfft, irfft, fftfreq
    
    # Filter specs
    #-----------------------------------------
    alld = np.zeros(dat.shape)
    highcut
    lowcut
    
    # loop through each time point and apply filter
    #------------------------------------------
    
    for i in range(dat.shape[0]):
        d = dat[i,:]
        
        f_signal = rfft(d)
        f_signal[0:(2*lowcut+1)] = 0
        f_signal[(2*highcut+1):len(f_signal)] = 0
       
        alld[i,:] = irfft(f_signal)

    return alld
    
    
#========================================================================
def fish_low(dat, fun='lin'):   # Filter low frequencies
#=======================================================================


    import numpy as np
    from scipy import stats, signal
    
    alld2 = np.zeros(dat.shape)
    
    if fun == 'lin':
        
        for i in range(dat.shape[0]):
            d2 = dat[i,:]
            slope, intercept, a,b,c = stats.linregress(np.linspace(0,len(d2)-1,len(d2)),d2)  
            d2 = d2 - (slope*d2 + intercept)   
            
            alld2[i,:]  = d2
            
    elif fun == 'filt':
        
        # Filter specs
        #----------------------------------------------------------------------------------------------
        fc  = .001           # cutoff frequency
        fs  = 2.7            # sampling frequency
        nfc = fc / (fs / 2)  # normalised cutoff 

        # Generate and apply filter
        #----------------------------------------------------------------------------------------------
        b, a   = signal.butter(5, nfc, 'high')
        alld2   = signal.filtfilt(b,a, dat, axis=1,method='gust')

    return alld2


    
#========================================================================
def fish_max_min(dat, window):   # Find max mins
#========================================================================
# This function calculates the minimum points across a sliding window for each cell - it then calculates the max of these mins - use this as a way to threshold out noisy cells (high minimum relative to overall maximum)
    
    import numpy as np
    
    maxmin_allcells =  np.zeros(dat.shape[0])
    windows = int(dat.shape[1]/window)
    
    # loop through all cells, sliding window over reshaped data to find the minimum value of each 9th of the data - then find the max of these mins
    #---------------------------------------------
            
    for i in range(dat.shape [0]):
        trace = dat [i,:]
        rshape = trace.reshape(windows, int(len(trace)/windows))
        minwin = np.apply_along_axis(min, 1, rshape)
        maxmin1c = max(minwin)
        maxmin_allcells [i] = maxmin1c
    
    return(maxmin_allcells)          

#========================================================================
def fish_thresh(trace, coord, mxmin, Ftrace, Fdrop, thresh):   # Remove noisy cells
#========================================================================

    import os
    import numpy as np

    kepttr = trace[mxmin > thresh]
    excltr = trace[mxmin < thresh]
    keptco = coord[mxmin > thresh]
    exclco = coord[mxmin < thresh]
    
    print( 'Kept ' + str(kepttr.shape[0]) + ' cells')
    print('Filtered ' + str(excltr.shape[0]) + ' cells')
    
    
    # Save non-noisy cell coords and traces
    #-------------------------------------------------
    os.chdir(Fdrop)
    p_id = Ftrace[56:71] + '_' + 'real'
    np.save(p_id + '_'  'com_coord.npy', keptco)
    np.save(p_id + '_' + 'com_signal.npy', kepttr)
    np.save(p_id + '_' + 'int_com_signal.npy', kepttr.astype("int"))
    

    return (kepttr, keptco, excltr, exclco)

