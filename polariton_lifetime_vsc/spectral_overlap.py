import numpy as np
import glob
import sys
import columnplots as clp
from scipy.optimize import curve_fit

'''
Starting from equilibrium trajectories, plotting the IR spectrum inside or outside the cavity;
and then capturing the spectral overlapping numerically from normalized Lorentzian integral

The overlapping integral J = int d\omega L_P(omega) * L_Dark(omega) / N
where N = int d \omega L_P(omega) int d\omega' L_Dark(omega')

'''

def calc_avged_IR(path="eq_Rabi/E0_0e-4", pattern="simu_*.dac.txt", freq_start=2000, freq_end=2750):
    filenames = glob.glob("%s/%s" %(path, pattern))
    N = len(filenames)
    print("reading %s with %d files" %(path, N))
    size = np.size(np.loadtxt(filenames[0])[:,0])
    freq, sp = np.zeros(size), np.zeros(size)
    for filename in filenames:
        data = np.loadtxt(filename)
        freq += data[:,5]
        sp += (data[:,6] + data[:,7]) / 2e28
    freq /= N
    sp /= N
    # Here is the total IR, I need to truncate to a subspace
    df = freq[1] - freq[0]
    nstart, nend = int(freq_start // df), int(freq_end // df)
    sub_freq = freq[nstart:nend]
    sub_sp = sp[nstart:nend]
    # normalized spectrum
    sub_sp = sub_sp / np.sum(sub_sp) / df
    return sub_freq, sub_sp

def obtain_avged_IR(path="eq_Rabi/E0_0e-4", pattern="simu_*.dac.txt", freq_start=2000, freq_end=2750):
    savefile = path+"/spectrum_%d_%d.txt" %(freq_start, freq_end)
    try:
        data = np.loadtxt(savefile)
        sub_freq, sub_sp = data[:,0], data[:,1]
        print("Loaded data from saved file")
    except:
        print("Loading data from raw files")
        sub_freq, sub_sp = calc_avged_IR(path=path, pattern=pattern, freq_start=freq_start, freq_end=freq_end)
        # save the calculated result to the file
        data = np.zeros((np.size(sub_freq), 2))
        data[:,0] = sub_freq
        data[:,1] = sub_sp
        np.savetxt(savefile, data)
    return sub_freq, sub_sp

def spectral_overlap(path_incav="eq_Rabi/E0_2e-4", freq_start=2000, freq_end=2750, freq_middle=2327, PlotIt=False, path_outcav="eq_Rabi/E0_0e-4"):
    '''
    Return: xs, ys, J_LP, J_UP,
    where xs = [freq]*3, ys = [sp_outcav, sp_incav, sp_overlap], and J denote integral [cm]
    '''
    # outside cavity
    freq_outcav, sp_outcav = obtain_avged_IR(path=path_outcav)
    # inside cavity
    freq_incav, sp_incav = obtain_avged_IR(path=path_incav)
    # obtain spectral overlap
    data = np.zeros((np.size(freq_incav), 2))
    data[:,0] = sp_outcav
    data[:,1] = sp_incav
    sp_overlap = np.min(data, axis=-1)
    # Whether to plot or not
    xs = [freq_outcav, freq_incav, freq_incav]
    ys = [sp_outcav, sp_incav, sp_overlap]
    if PlotIt:
        plot_multiple_IR(xs, ys)
    # Finally, we calculate the integrated spectral overlap at LP and UP region seperately
    df = freq_incav[1] - freq_incav[0]
    nmiddle = int((freq_middle - freq_start)//df)
    print("nmiddle = %d, freq_middle = %d" %(nmiddle, freq_middle))
    N_LP = 1.0 #df * np.sum(sp_incav[0:nmiddle]) * np.sum(sp_outcav[0:nmiddle]) * df
    J_LP = df * np.sum(sp_incav[0:nmiddle] * sp_outcav[0:nmiddle]) / N_LP
    N_UP = 1.0 #df * np.sum(sp_incav[nmiddle:]) * np.sum(sp_outcav[nmiddle:]) * df
    J_UP = df * np.sum(sp_incav[nmiddle:] * sp_outcav[nmiddle:]) / N_UP
    print("Path = %s, spectral integral for J_LP = %.3E [cm], J_UP = %.3E [cm]" %(path_incav, J_LP, J_UP))
    return xs, ys, J_LP, J_UP

def plot_multiple_IR(xs, ys):
    ax = clp.initialize(1, 1)
    clp.plotone(xs, ys, ax, showlegend=False)
    clp.adjust()

if __name__ == '__main__':
    spectral_overlap(path_incav="eq_Rabi/E0_1e-4", PlotIt=True)
