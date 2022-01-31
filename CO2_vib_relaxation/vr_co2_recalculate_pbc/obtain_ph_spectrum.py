'''
Combine all output to a single output file
'''
import numpy as np
import sys
import glob
import os
from scipy import signal
import math
from scipy import fftpack

dtfs = 2.0

def auto_correlation_function_simple(x):
    n = x.size
    if n % 2 == 0:
        x_shifted = np.zeros(n*2)
    else:
        x_shifted = np.zeros(n*2-1)
    x_shifted[n//2 : n//2+n] = x
    # Convolute the shifted array with the flipped array, which is equivalent to performing a correlation
    autocorr_full = (signal.fftconvolve(x_shifted, x[::-1], mode='same')[-n:]/ np.arange(n, 0, -1))
    # Truncate the autocorrelation array
    autocorr = autocorr_full[0:n//2]
    return autocorr

def fft3(x, N=None):
    # Adding zeros to the end of x
    if N is not None:
        n = N
    else:
        n = np.size(x)
    lineshape = fftpack.dct(x, type=1, n=n)
    freq_au = np.linspace(0, 0.5/dtfs * 1e15, n)
    # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    # Calculate spectra
    #field_description =  freq_au**2
    field_description =  freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra

def calc_spectrum(data):
    phx, phy = data[:, -6], data[:, -2]
    dacf_x = auto_correlation_function_simple(phx)
    dacf_y = auto_correlation_function_simple(phy)
    dacf_x_freq, dacf_x_sp = fft3(dacf_x)
    dacf_y_freq, dacf_y_sp = fft3(dacf_y)
    return dacf_y_freq, dacf_x_sp + dacf_y_sp

def obtain_avg(path):
    print("working with", path)
    freq, sp = calc_spectrum(np.loadtxt(path+"/simu_1.out"))
    data_raw = np.zeros((np.size(sp), 2))
    files = glob.glob(path+"/simu_*.out")
    Ntraj = 0
    for fname in files:
        print(fname, os.path.getsize(fname))
        if os.path.getsize(fname) >= 3620730:
            freq, sp = calc_spectrum(np.loadtxt(fname))
            data_raw[:, 0] += freq
            data_raw[:, 1] += sp
            Ntraj += 1
    data_raw /= Ntraj
    print("%d trajs saved" %Ntraj)
    np.save(path+"/ph_sp_avg", data_raw)

if __name__ == '__main__':
    obtain_avg(path=sys.argv[-1])
