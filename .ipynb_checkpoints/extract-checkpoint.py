import os
import glob 
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#PROCESS
#--------------
#---------------
#=======================================================================
def fish_load(Fs2p, Fsave, fish, experiment, date): 
#=======================================================================

    """
     This function looks in the data folder for suite2p plane files and saves extracted cell traces/coordinates into savepath
    
    Inputs:
    Fs2p (string): path with suite2p folders inside
    Fsace (string): save path 
    fish (string): file name for fish from raw data
    experiment (string): parent path for fish raw data
    date (string): date of data collection
    
    Saves:
    ...all_trace: a numpy file containing a ncells x ntimepoints array of data - from all active cells as defined in suite2p, ordered by plane
    ...all_coords: a numpy file containing a ncells x X x Y x Z coordinates
    
    """

    dirlist = os.listdir(Fs2p)
    Fsave = Fsave

    # Find planes of suite2p output
    #------------------------------
    r       = re.compile('^plane[0-9].*')
    planelist = list(filter(r.match, dirlist))
    planelist.sort()
    print('Found ' + str(len(planelist)) + ' planes') 
    
    # Compile coordinates and trace files for all planes into lists
    #---------------------------------------------------------------------
    coord = list((range(len(planelist))))
    trace = list((range(len(planelist))))
    cells = list((range(len(planelist))))
    for i in range(len(planelist)):
        os.chdir(Fs2p + os.sep + "plane" + str(i)) 
        allcells = np.load("iscell.npy")  
        fl = np.load("F.npy") 
        
        stats = np.load("stat.npy", allow_pickle=True) 
        xy = np.zeros((len(stats),2))  
        for  j in range (len(stats)):          
            xy [j,] = stats [j] ['med']
        xyz = np.concatenate([xy, np.full((len(fl), 1), i)], axis = 1)
        coord[i] = xyz
        trace[i] = fl
    
    # Concatenate separate arrays for coordinate and xy file into two arrays 
    #--------------------------------------------------------------------
    com_coord = np.concatenate([coord[i] for i in range(len(planelist))])
    com_signal = np.concatenate([trace[i] for i in range(len(planelist))]) 
    print('Found ' + str(com_coord.shape[0]) + ' cells')
    
    # Save as three separate files (int file is for R) 
    #-------------------------------------------------


    p_id = experiment + '-' + fish[1:3] + '_' + fish[fish.find(date[7:]) + 1 + len(date[7:]):fish.find('dpf') - 2] + '_' + fish[fish.find('se'):fish.find('se')+8] + fish[fish.find('dpf')-1:fish.find('dpf')+3] + '_' + fish[fish.find(fish[1:3])+3:fish.find(date[7:])-1] + '_' + fish[fish.find('run'):fish.find('run')+6]

    np.save(Fsave + os.sep + p_id + '_'  'allcoord.npy', com_coord)
    np.save(Fsave + os.sep + p_id + '_' + 'alltrace.npy', com_signal)
    print('Saved ' + str(fish) + ' in '+ str(experiment))
    
    
#========================================================================
def fish_filter(dat, highcut, lowcut):   
#========================================================================

    """
     This function removes high frequency components from the trace as noise. It can also remove low frequency components - as a bandpass filter. This is necessary for the max/min filtering step. 

    Inputs:
    dat (np array): cells x timepoints of raw traces
    highcut (int): high frequency cutoff
    lowcut (int): low frequency cutoff
    nplt (int): number of traces to plot

    
    Returns:
    alld (np array): cells x timepoints of filtered traces
    norm (np array): cells x timepoints of normalised traces
    

    """

    from scipy.fftpack import rfft, irfft, fftfreq
    import random as rand
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
    
    # now normalise trace
    #------------------------------------------
    def fish_norm(alld):   
        return(alld/np.mean(alld))
    norm = np.apply_along_axis(fish_norm, 1, (alld + 3000))

    return (alld, norm)
    
    

#========================================================================
def fish_max_min(dat, window): 
#========================================================================

    """
     This function calculates the minimum points across a sliding window for each cell - it then calculates the max of these mins. We use this as a way to threshold out noisy cells (high minimum relative to overall maximum). 
    
    Inputs:
    dat (np array): cells x timepoints of normalised traces
    window (int): window size to calculate max mins

    Returns:
    maxmin_allcells (np array): vector of maxmins for each cell

    """

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
def fish_thresh(Fsave, name, trace, coord, mxmin, thresh, save):   
#========================================================================

    """
     This function removes noisy neurons with low max-mins according to defined threshold, and plots the results. Saving is an option. 
    
    Inputs:
    Fsave (str): save path
    name (str): save name 
    trace (np array): cells x timepoints of normalised traces
    coord (np array): cells x X x Y x Z coordinates
    mxmin (np array): vector of maxmins for each cell
    thresh (int): threshold for is cell vs not cell
    
    Saves:
    keptco (np array): cells x X x Y x Z coordinates
    kepttr (np array): cells x trace timepoints


    """
    
    
    kepttr = trace[mxmin > thresh]
    excltr = trace[mxmin < thresh]
    keptco = coord[mxmin > thresh]
    exclco = coord[mxmin < thresh]
    print( 'Kept ' + str(kepttr.shape[0]) + ' cells')
    print('Filtered ' + str(excltr.shape[0]) + ' cells')
    
    # Plot kept and excluded cells - adjust threshold
    #---------------------------------------------------------------------------
    print(name)
    plt.figure(figsize= (15,15))
    removed = plt.scatter(exclco[:,0], exclco[:,1], s=6)
    kept = plt.scatter(keptco[:,0], keptco[:,1], s=6)
    plt.legend((kept, removed), ('kept', 'removed'), loc = 'lower right', fontsize = 15, scatterpoints = 1, markerscale = 5)
    plt.show()
    
    if save == 'yes':
    
        # Save real cell coords and traces in new folder and save each plane as well
        #--------------------------------------------------------------------------------
        p_id = name[:name.find('run')+6]
        np.save(Fsave  + '/' + p_id + '_realcoord.npy', keptco)
        np.save(Fsave  + '/' + p_id + '_realtrace.npy', kepttr)
    

    
    
#=====================================================================================================
def bcl_function_parameters(wdt, savepath, experiment, name, data, lamb, varB, varC, Cmean, frequency, gausfilt, mode):
#=====================================================================================================
    """
     This function runs the BCL model on trace data and can either save or plot the results. BCL is Bayesian Cleaner - a Hidden Markov Model for estimating hidden states of the neuron from fluorescence traces. 
    
    Inputs:
    wdt (int): spiking probability
    lamb (int): calcium decay 
    varB (int): baseline variance
    varC (int): calcium variance
    Cmean (int): mean calcium
    frequency (int): frequency
    gausfilt (int): gaussian filter parameter for smoothing
    mode (str): define whether to run model and save output ('save') or just plot a few traces ('see')
    
    Returns:
    carray (np array): cells x time of estimated calcium 
    sksarray (np array): cells x time of binarised states of neurons
    Barray (np array): cells x time of estimated baseline for each neuron

    """

            
    def bcl_model(wdt, trace, lamb, varB, varC, Cmean, frequency, gausfilt):
        from scipy import fftpack
        import math
        from math import log, pi
        from scipy.ndimage import gaussian_filter1d
        import random as rand 

        #Preprocess Function
        t1 = ((trace + 500)/ (500 + np.mean(trace[np.where(trace < np.quantile(trace, 0.08, axis = 0))[0]]))) - 1 
        #normalise trace 
        y = gaussian_filter1d(t1, gausfilt, axis = 0)
        difft = diff(t1) 
        varX = get_variance_of_the_decreases(difft)

        #Declare Variables
        N = len(y)
        B = np.zeros(y.shape[0])
        c =  np.zeros(y.shape[0])
        sks = np.zeros(y.shape[0])
        B[0] = np.mean(y[0:500])  #make baseline starting point
        dff = np.zeros(y.shape[0])
        loglik = 0
        dt = float(1) / frequency

        #each time point, chance calcium event vs baseline shift
        #bcl outputs a timeseries for baseline and calcium

        for t in range(1,N):

            #new calcium value, if no calcium spike at current t
            #LET CALCIUM DECAY - IF CALCIUM AT T-1, THEN DECAY - IF NO CALCIUM AT T -1, CNEW REMAINS 0 
            #calcium at previous time point * exponential - if calcium = 0, then cnew = 0
            cnew = c[t - 1] * np.exp(-lamb* dt)  

            #new baseline value, if no calcium spike at current t - i.e. IF BASELINE IS DECREASING 
            #SIGNAL AT T - MODELLED CALCIUM (with variance) + BASELINE BEFORE (with variance)
            #(Baseline at t-1 * variance of decreases) + (Signal at t, - cnew * baseline variance + frequency) 
            #/ variances of overall data 
            Bnew = (varX * B[t - 1] + varB * dt * (y[t] - cnew)) / (varX + varB * dt)

            #p of timestep being explained by baseline, not spike
            logp0 = log(1 - wdt) - 0.5 * log(2 * pi) - 0.5 * log(varX + varB * dt) - (y[t] - cnew - B[t - 1]) ** 2 / (2 * varX + 2 * varB * dt)

            #new calcium value, if calcium spike
            # (signal at t, - baseline at t-1, - decayed calcium) + (mean calcium + calcium decay)
            # / (1 + variance of data)
            cspike = Cmean + cnew + (y[t] - cnew - B[t - 1]) / (1 + varB * dt / varC + varX / varC)
            cspike = np.clip(cspike, 0, 10000)

            #new baseline value, if calcium spike
            #Baseline at t-1, + variance of baseline*freq/variance of calcium * (calcium spike - decay - mean cal)
            Bspike = B[t - 1] + varB * dt / varC * (cspike - cnew - Cmean)

            #p of timestep being explained by spike
            logp1 = log(wdt) - 0.5 * log(2 * pi) - 0.5 * log(varX + varB * dt + varC) - (y[t] - cnew - B[t - 1] - Cmean)**2 / (2 * varX + 2 * varB * dt + 2 * varC)


            #compares logp1 vs logp0 
            if logp1 < logp0:
                c[t] = cnew
                B[t] = Bnew
                loglik = loglik + logp0

            else:
                c[t] = cspike
                B[t] = Bspike
                loglik = loglik + logp1

        dfftsks = diff(c)
        sks[np.where(np.asarray(dfftsks) > 0)] = 1
        return(c, sks, B, y)


    def diff(timeseries):
    #minus each timestep by the timestep before it 
    #so you can estimate variance from one step to next
        diff_timeseries = []
        for index in range(1,len(timeseries)):
            difference = timeseries[index] - timeseries[index - 1]
            diff_timeseries.append(difference)
        return diff_timeseries


    def lowpass_filter(trace,frequency_cutoff): 
    #fft gives symmetrical trace - half represented by real half by imaginary
    #ft represented with imaginary components, in 3d as 3d spiral - 2d section is sine wave
    #high frequencies represented in middle, low either side (reflected as real and imaginary)
    #so discard half of imaginary series by blocking out middle half

        from scipy import fftpack
        fast_fourier_transform = fftpack.fft(trace) # Take The Fourier Transform Of The Trace (power y axis, freq x axis)
        fast_fourier_transform[frequency_cutoff : len(fast_fourier_transform)-(frequency_cutoff-1)] = 0 
        #high frequencies represented in middle, low either side (reflected as real and imaginary) - so discard half of imaginary series by blocking out middle half
        filtered_signal = fftpack.ifft(fast_fourier_transform)  # Run The Inverse Fourier Transform To Get Back To A Signal
        real_filtered_signal = np.real(filtered_signal) #throw away imaginary
        return real_filtered_signal


    def get_variance_of_the_decreases(difft): 
    #variances of decreases, timepoint is whenever there is a drop 
    #variance of decrease says how much it decays ie. calcim signal
        import math
        from math import log, pi

        number_of_decreases = 0 
        squared_sum_of_decreases = 0
        for timepoint in difft:
            #if next point is a decrease ie. decay, model the variance of decay
            if timepoint < 0:
                number_of_decreases += 1
                squared_sum_of_decreases += (timepoint ** 2)
        variance = squared_sum_of_decreases / number_of_decreases
        variance = math.sqrt(variance)
        return variance

    
    import random as rand 
    c, N, B, sks, loglik, dt = [],[],[],[],[],[]
    
    if mode == 'save':
        Barray, carray, sksarray = np.zeros(data.shape),np.zeros(data.shape),np.zeros(data.shape)

        for i in range(data.shape[0]):
            trace = data[i]
            
            c, sks, B, y = bcl_model(wdt, trace, lamb, varB, varC, Cmean, frequency, gausfilt)
            Barray[i] = B
            carray[i] = c
            sksarray[i] = sks 
    
        #np.save(savepath + 'experiment + os.sep + name[:name.find('run')+6] + '_' + 'modelcal.npy', carray)  
        np.save(savepath + experiment + os.sep + name + '_' + 'regbinarised.npy', sksarray)  
        
        return carray,sksarray, Barray

    elif mode == 'see':
        from scipy.ndimage import gaussian_filter1d

        rdm = rand.sample(range(0, data.shape [0]), 10)
    
        for i in rdm:
            trace = data[i]
            C, sks, B, y = bcl_model(wdt, trace, lamb, varB, varC, Cmean, frequency, gausfilt) 
        
            ind = 2000
            plt.figure(figsize = (40,4))    
            plt.plot(y[:ind], c = 'deepskyblue', linewidth = 7)
            for x in np.where(sks[:ind] == 1)[0]: plt.vlines(x=x,  alpha = 1, ymin = -8, ymax = -2, color = 'red', linewidth = 2)
            plt.ylim(-5,np.max(y)+3)
            plt.show()

#=======================================================================
def deltaff(mat, percentile):
#=======================================================================
    """
     This function calculates delta f/f for the trace by estimating a baseline as the nth lowest percentile of the data.
    
    Inputs:
    old_trace (np array): cells x fluorescence time points, raw trace
    percentile (int): percentile of data to take as baseline
    
    Returns:
    new_trace (np array): cells x fluorescebce time points, normalised trace

    """
    def dff_calc(old_trace, percentile):
        mini = np.min(old_trace)
        trace = old_trace - np.min(old_trace)
        newtrace = trace
        baseline = np.quantile(trace, percentile, axis= 0)
        newtrace[np.where(trace <= baseline)] = 0 #set all values below baseline to 0
        pos_index = np.where(trace > baseline) #indeces of values above baseline
        newtrace[pos_index] = (trace[pos_index] - baseline)/baseline
        return(newtrace)
    
    dff_mat = np.zeros(mat.shape)
    for i in range(mat.shape[0]):
        dff_mat[i] = dff_calc(mat[i], percentile)

    return(dff_mat)



#SAVE
#--------------
#---------------
#=======================================================================
def fish_backup(Fs2p, backup, experiment, fish, date, add_day, dir_mode): 
#=======================================================================
    """
     This function backsup all datasets to hard drive. 
    
    Inputs:
    Fs2p (str): path to suite2p files
    backup (str): backup path
    experiment (str): name of experiment folder
    fish (str): name of current dataset
    date (str): date of current dataset
    add_day (str): add day to end of naming structure
    dir_mode (str): how to structure the directory when backing up, 'new' = creates new folder, 'add_condition' = adds new condition to existing folder, 'add_day' = adds day to existing folder

    """


    import shutil
    os.chdir(Fs2p)
    planes = sorted(glob.glob("plane*")) 

    #define subdirectories
    #--------------------------------------------------------

    day = fish[fish.find('dpf')-1:fish.find('dpf')+3]

    if add_day == 'yes':
        fold1 = experiment + '-' + day + '-' + fish[1:3] + '/'

    else:
        fold1 = experiment + '-' + fish[1:3] + '/'
    fold2 = fish[fish.find(date[7:]) + 1 + len(date[7:]):fish.find('dpf') - 2] 
    fold3 = fish[fish.find('se'):fish.find('se')+8] + fish[fish.find('dpf')-1:fish.find('dpf')+3]

    #create new subdirectories (only do once for each fish)
    #--------------------------------------------------------
    if dir_mode == 'new':
        os.chdir(backup + '/Project/')
        os.mkdir(fold1)
        os.chdir(backup + '/Project' + os.sep + fold1)
        os.mkdir(fold2)
        os.chdir(backup + '/Project' + os.sep + fold1 + os.sep + fold2)
        os.mkdir(fold3)

    elif dir_mode == 'add_condition':
        fold1 = fold1
        fold2 = fold2
        fold3 = fold3

    elif dir_mode == 'add_day':
        os.chdir(backup + '/Project' + os.sep + fold1 + os.sep + fold2)
        os.mkdir(fold3)

    #define subdirectories
    #--------------------------------------------------------            

    p_id = experiment + '-' + fish[1:3] + '_' + fish[fish.find(date[7:]) + 1 + len(date[7:]):fish.find('dpf') - 2] + '_' + fish[fish.find('se'):fish.find('se')+8] + day + '_' + fish[fish.find(fish[1:3])+3:fish.find(date[7:])-1] + '_' + fish[fish.find('run'):fish.find('run')+6]


    os.chdir(Fs2p + os.sep + 'plane0' + os.sep + 'reg_tif')
    tifs = sorted(glob.glob("*tif")) 
    for i in range(0,10):
        shutil.move(Fs2p + "/plane" + str(i) + "/ops.npy", backup + 'Project/' + fold1 + os.sep + fold2 + os.sep + fold3 + os.sep + p_id + '_plane' + str(i) + '_ops.npy' )
        shutil.move(Fs2p + "/plane" + str(i) + "/stat.npy", backup + 'Project/' + fold1 + os.sep + fold2 + os.sep + fold3 + os.sep + p_id + '_plane' + str(i) + '_stat.npy')
        for x in range(len(tifs)):
            shutil.move(Fs2p + "/plane" + str(i) + os.sep + 'reg_tif' + os.sep + tifs[x], backup + 'Project/' + fold1 + os.sep + fold2 + os.sep + fold3 + os.sep + p_id + '_plane' + str(i) + '_reg' + tifs[x][4:7] + '.tif')







