import numpy as np
import columnplots as clp
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exp_func(x, k, a):
    return a*np.exp(-x/k)

def fit_exp(x, y):
    popt, pocv = curve_fit(exp_func, x, y)
    return 1.0 / popt[0]

def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories for path %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def gather_IR_data(paths):
    freqs, sps = [], []
    for path in paths:
        data = obtain_avg_data(path=path, pattern="simu_*.dac.txt")
        freq, sp = data[:,5], (data[:,6] + data[:,7])/2e28
        freqs.append(freq)
        sps.append(sp)
    return freqs, sps

def gather_relaxation_trajectory(paths):
    ts, trajs = [], []
    for path in paths:
        data = obtain_avg_data(path=path, pattern="simu_*.VCO_summed_td.txt")
        t, traj = data[:,0], data[:, 1]
        ts.append(t*1e-3) # fs to ps
        trajs.append(traj / 10.0 - 0.000950037) # averaged over 10 molecules and minus 300K in a.u.
    return ts, trajs

'''
Here I plot the detuning dependence of vibrational energy relaxation rate
'''
paths_outcav = ["try_VER/Exc_10_E0_0e-4/"]
freqs = np.array([2000, 2100, 2150, 2200, 2250, 2300, 2320, 2350, 2400, 2450, 2500, 2600])
paths_incav  = ["try_VER/Freq_%d" %x for x in freqs]
paths_incav[6] = "try_VER/Exc_10_E0_2e-4/"

def plot_results():
    # obtain time-resolved trajectory
    ts_outcav, trajs_outcav = gather_relaxation_trajectory(paths=paths_outcav)
    ts_incav, trajs_incav = gather_relaxation_trajectory(paths=paths_incav)

    # calculate decay rates
    decay_rates_incav = np.array([fit_exp(x, y) for x, y in zip(ts_incav, trajs_incav)])
    decay_rate_outcav = np.array([fit_exp(x, y) for x, y in zip(ts_outcav, trajs_outcav)] * len(freqs))

    # obtain the bare molecule lineshape
    #freqs_ir, sps_ir = gather_IR_data(paths=["co2_eq/E0_0e-4/"])

    ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618, LaTeX=True, fontsize=12)

    clp.plotone([freqs-2327]*2, [decay_rates_incav, decay_rate_outcav], ax, colors=["r--o", "k--s"],
        labels=["inside cavity", "outside cavity"], xlabel="cavity mode detuning [cm$^{-1}$]",
        ylabel="fitted decay rate [ps$^{-1}$]", ylim=[0, 0.12])

    #ax2 = ax.twinx()
    #clp.plotone(freqs_ir, sps_ir, ax2, colors=["g"], showlegend=False, xlim=[1900, 2700], ylim=[-0.1, 1.])
    clp.adjust(tight_layout=True, savefile="detuning_dependence.pdf")


if __name__ == "__main__":
    plot_results()
