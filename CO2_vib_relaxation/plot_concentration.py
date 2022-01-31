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

def lin_func(x, k, a):
    return k*x + a

def fit_lin(x, y):
    popt, pocv = curve_fit(lin_func, x, y)
    return lin_func(x, *popt)


def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories in path %s" %(len(filenames), path))
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
Here I plot the concentration dependence of vibrational energy relaxation rate
'''
concentration = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
paths_incav = ["vr_co2_concentration/c_%.1f/E0_2e-4/" %x for x in concentration]
paths_incav[-1] = "try_VER/Exc_10_E0_2e-4/"
paths_outcav = ["vr_co2_concentration/c_%.1f/E0_0e-4/" %x for x in concentration]
paths_outcav[-1] = "try_VER/Exc_10_E0_0e-4"

paths_ir = ["vr_co2_concentration/eq_c_%.1f/E0_2e-4/" %x for x in concentration]
paths_ir[-1] = "co2_eq/E0_2e-4/"

def plot_results():
    # obtain time-resolved trajectory
    ts_incav, trajs_incav = gather_relaxation_trajectory(paths=paths_incav)
    ts_outcav, trajs_outcav = gather_relaxation_trajectory(paths=paths_outcav)

    # calculate decay rates
    decay_rates_incav = np.array([fit_exp(x, y) for x, y in zip(ts_incav, trajs_incav)])
    decay_rate_outcav = np.array([fit_exp(x, y) for x, y in zip(ts_outcav, trajs_outcav)])

    # obtain the IR spectrum inside cavity for different concentration of molecules
    freqs_ir, sps_ir = gather_IR_data(paths=paths_ir)
    sps_ir = [x + 2*i for i, x in enumerate(sps_ir)]

    axes = clp.initialize(1, 2, width=6.3, height=4.3*0.618, LaTeX=True, fontsize=12,  labelthem=True, labelthemPosition=[-0.02, 1.05])

    clp.plotone(freqs_ir, sps_ir, axes[0], xlim=[2000, 2600], showlegend=False,
        xlabel="frequency [cm$^{-1}$]", ylabel="$n(\omega)\\alpha(\omega)$ [arb. units]", ylim=[-0.1,19])
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(axes[0].lines))]
    for i,j in enumerate(axes[0].lines):
        j.set_color(colors[i])
    #axes[0].axvline(x = 2320, ymin=0.0, ymax=3.0, ls='--', c="b")
    for i in range(5):
        str = "$c$=%d" %(i*20+20)
        axes[0].text(2010, 2*i+0.2, str + "$\%$", fontsize=8)

    # plot the UP and LP lines
    dx = freqs_ir[0][1] - freqs_ir[0][0]
    low_end, middle_end, high_end = int(2188/dx), int(2320/dx), int(2600/dx)
    freq = freqs_ir[0]
    freq_UPrange = freq[middle_end:high_end]
    freq_LPrange = freq[low_end:middle_end]
    UP_freqs = np.array([freq_UPrange[np.argmax(sp[middle_end:high_end])] for sp in sps_ir])
    LP_freqs = np.array([freq_LPrange[np.argmax(sp[low_end:middle_end])] for sp in sps_ir])
    UP_signal = np.array([np.max(sp[middle_end:high_end]) for sp in sps_ir])
    LP_signal = np.array([np.max(sp[low_end:middle_end]) for sp in sps_ir])
    xs = [LP_freqs, UP_freqs]
    ys = [LP_signal, UP_signal]
    clp.plotone(xs, ys, axes[0], colors=["k--", "k--"], lw=0.8, showlegend=False)
    axes[0].text(2240, 9, "LP", c='r')
    axes[0].text(2450, 9, "UP", c='r')
    axes[0].text(2150, 0.5, "$^{14}$C=$^{18}$O")

    # Plot the Rabi splitting as a function of sqrt(c)
    Rabi = UP_freqs - LP_freqs
    #Rabi[0] -= 2
    ys = [Rabi, fit_lin(np.sqrt(concentration), Rabi)]
    xs = [np.sqrt(concentration)]*2
    # add inset to figure a.
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    inset_ax2 = inset_axes(axes[0], width="44%", height=0.42)
    clp.plotone(xs, ys, inset_ax2, colors=["ro", "k--"], lw=1., showlegend=False)
    inset_ax2.set_xlabel('$\sqrt{c/c_0}$', fontsize=7)
    inset_ax2.set_ylabel('$\Omega_N$ [cm$^{-1}$]', fontsize=7)
    inset_ax2.tick_params(axis="x", labelsize=7)
    inset_ax2.tick_params(axis="y", labelsize=7)

    clp.plotone([concentration]*2, [decay_rates_incav, decay_rate_outcav], axes[1], colors=["r--o", "k--s"],
        labels=["inside cavity", "outside cavity"], xlabel="$^{12}$C$^{16}$O$_2$ concentration $c/c_0$",
        ylabel="fitted decay rate [ps$^{-1}$]", ylim=[0, 0.14])

    clp.adjust(tight_layout=True, savefile="concentration_dependence.pdf")

if __name__ == "__main__":
    plot_results()
