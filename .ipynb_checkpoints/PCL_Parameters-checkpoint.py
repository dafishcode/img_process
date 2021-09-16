import sys
import numpy as np
import math
import pandas
from scipy import fftpack
from math import log, pi
from matplotlib import pyplot as plt

sys.path.append(r"C:\Users\matth\OneDrive\Documents\Python Code\Matts Modules")
import matts_admin_functions
import matts_network_functions
import matts_matrix_functions


#R Functions Duplicated
def diff(timeseries):
    diff_timeseries = []
    for index in range(1,len(timeseries)):
        difference = timeseries[index] - timeseries[index - 1]
        diff_timeseries.append(difference)

    return diff_timeseries

def lowpass_filter(trace,frequency_cutoff):
    fast_fourier_transform = fftpack.fft(trace)             # Take The Fourier Transform Of The Trace
    fast_fourier_transform[frequency_cutoff : len(fast_fourier_transform)-(frequency_cutoff-1)] = 0
    filtered_signal = fftpack.ifft(fast_fourier_transform)  # Run The Inverse Fourier Transform To Get Back To A Signal
    real_filtered_signal = np.real(filtered_signal)
    return real_filtered_signal




def get_variance_of_the_decreases(trace):
    number_of_decreases = 0
    sqaured_sum_of_decreases = 0

    for timepoint in trace:
        if timepoint < 0:
            number_of_decreases += 1
            sqaured_sum_of_decreases += (timepoint ** 2)

    variance = sqaured_sum_of_decreases / number_of_decreases
    variance = math.sqrt(variance)
    return variance

def bcl_function_parameters(trace,wdt,lamb,varB,varC,initial_baseline,Cmean,frequency,fourier_cutoff):

    #Preprocess Function
    t1 = trace / np.mean(trace)
    y = lowpass_filter(t1, fourier_cutoff) #2000 too high, #200 too low

    difft = diff(t1)
    varX = get_variance_of_the_decreases(difft)

    #Declare Variables
    N = len(y)
    B = matts_admin_functions.create_empty_list(N)
    c = matts_admin_functions.create_list_of_zeros(N)
    sks = matts_admin_functions.create_list_of_zeros(N)
    B[0] = np.mean(y[0:100])
    loglik = 0
    dt = float(1) / frequency

    for t in range(1,N):
        cnew = c[t - 1] * np.exp(-lamb* dt)
        Bnew = (varX * B[t - 1] + varB * dt * (y[t] - cnew)) / (varX + varB * dt)

        logp0 = log(1 - wdt) - 0.5 * log(2 * pi) - 0.5 * log(varX + varB * dt) - (y[t] - cnew - B[t - 1]) ** 2 / (2 * varX + 2 * varB * dt)

        cspike = Cmean + cnew + (y[t] - cnew - B[t - 1]) / (1 + varB * dt / varC + varX / varC)
        Bspike = B[t - 1] + varB * dt / varC * (cspike - cnew - Cmean)

        logp1 = log(wdt) - 0.5 * log(2 * pi) - 0.5 * log(varX + varB * dt + varC) - (y[t] - cnew - B[t - 1] - Cmean)**2 / (2 * varX + 2 * varB * dt + 2 * varC)

        if logp1 < logp0:
            c[t] = cnew
            B[t] = Bnew
            loglik = loglik + logp0

        else:
            c[t] = cspike
            B[t] = Bspike
            sks[t] = 1
            loglik = loglik + logp1

    return c,B,sks,y





