#=======================================================================
def diff(timeseries): # Load timeseries
#=======================================================================
#minus each timestep by the timestep before it 
#so you can estimate variance from one step to next
    diff_timeseries = []
    for index in range(1,len(timeseries)):
        difference = timeseries[index] - timeseries[index - 1]
        diff_timeseries.append(difference)
    return diff_timeseries


#=======================================================================
def lowpass_filter(trace,frequency_cutoff): # Load trace
#=======================================================================
#fft gives symmetrical trace - half represented by real half by imaginary
#ft represented with imaginary components, in 3d as 3d spiral - 2d section is sine wave
#high frequencies represented in middle, low either side (reflected as real and imaginary)
#so discard half of imaginary series by blocking out middle half

    fast_fourier_transform = fftpack.fft(trace) # Take The Fourier Transform Of The Trace (power y axis, freq x axis)
    fast_fourier_transform[frequency_cutoff : len(fast_fourier_transform)-(frequency_cutoff-1)] = 0 
    #high frequencies represented in middle, low either side (reflected as real and imaginary) - so discard half of imaginary series by blocking out middle half
    filtered_signal = fftpack.ifft(fast_fourier_transform)  # Run The Inverse Fourier Transform To Get Back To A Signal
    real_filtered_signal = np.real(filtered_signal) #throw away imaginary
    return real_filtered_signal


#=======================================================================
def get_variance_of_the_decreases(difft): # Load difference timeseries
#=======================================================================
#variances of decreases, timepoint is whenever there is a drop 
#variance of decrease says how much it decays ie. calcim signal

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