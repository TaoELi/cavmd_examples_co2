import numpy as np
import glob
import sys
import columnplots as clp
from scipy.optimize import curve_fit
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories" %(len(filenames)))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def obtain_IR(path):
    data = obtain_avg_data(path=path, pattern="simu_*.dac.txt")
    x, y = data[:,5], (data[:,6] + data[:,7])/2e28
    return x, y

def func(x, k, a, b):
    return a * np.exp(-x/k) + b

def fit_exponential(x, y):
    popt, pocv =  curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def func_linear(x, k, b):
    return k*x + b

def fit_linear(x, y):
    popt, pocv =  curve_fit(func_linear, x, y)
    yfit = func_linear(x, *popt)
    print("popt i, s", popt)
    return yfit, popt[0]


def prepare_ph_data(path="pumping_Rabi/E0_2e-4_Exc_UP_Amp_6e-3", omega=2320, pattern="simu_*.out", dtfs=2):
    filenames = glob.glob("%s/%s" %(path, pattern))
    N = len(filenames)
    print("reading %s with %d files" %(path, N))
    size = np.size(np.loadtxt(filenames[0])[:,0])
    t, phx, phy = np.zeros(size), np.zeros(size), np.zeros(size)
    ph2au = 0.5 * 1.8897259885789**2 * omega**2 * 0.0000045563352812122295**2
    vph2au = 0.5 * 1.8897259885789**2 / 41.341374575751**2
    for filename in filenames:
        data = np.loadtxt(filename)
        t_loc, qx_loc, qy_loc = data[:,1], data[:,-6], data[:,-2]
        vx_loc = np.gradient(qx_loc, dtfs, axis=-1, edge_order=2)
        vy_loc = np.gradient(qy_loc, dtfs, axis=-1, edge_order=2)
        t += t_loc
        phx += qx_loc**2 * ph2au + vx_loc**2 * vph2au
        phy += qy_loc**2 * ph2au + vy_loc**2 * vph2au
    return t[:-1]/N, (phx[:-1] + phy[:-1])/N

def plot_ph_dynamics():

    # define cavloss
    # 30 60 100 300 600 1000 3000 5000 10000 1000000
    cav_lifetimes = [300, 600, 1000, 5000, 1000000]
    paths_LP_cavloss = ["pumping_cavloss/E0_2e-4_Exc_LP_CavLoss_%d_dilute" %t for t in cav_lifetimes]
    paths_UP_cavloss = ["pumping_cavloss/E0_2e-4_Exc_UP_CavLoss_%d_dilute" %t for t in cav_lifetimes]
    cav_loss_array = 1.0 / np.array(cav_lifetimes) * 1e3
    n_data = len(paths_LP_cavloss)

    au2omega0 = 0.01057069785241237 #  = 2320 cm-1
    # Calculate the LP data
    nstart, nend = 300, 20000
    xs_LP, ys_LP, xs_LP_fit, ys_LP_fit, taus_LP = [], [], [], [], []
    for path in paths_LP_cavloss:
        t, phtot = prepare_ph_data(path=path)
        yfit_ph, tau = fit_exponential(t[nstart:nend], phtot[nstart:nend])
        xs_LP += [t]
        ys_LP += [phtot/au2omega0]
        xs_LP_fit += [t[nstart:nend]]
        ys_LP_fit += [phtot[nstart:nend]/au2omega0]
        taus_LP.append(tau)
    taus_LP = np.array(taus_LP)

    # Calculate the UP data
    xs_UP, ys_UP, xs_UP_fit, ys_UP_fit, taus_UP = [], [], [], [], []
    for path in paths_UP_cavloss:
        t, phtot = prepare_ph_data(path=path)
        yfit_ph, tau = fit_exponential(t[nstart:nend], phtot[nstart:nend])
        xs_UP += [t]
        ys_UP += [phtot/au2omega0]
        xs_UP_fit += [t[nstart:nend]]
        ys_UP_fit += [phtot[nstart:nend]/au2omega0]
        taus_UP.append(tau)
    taus_UP = np.array(taus_UP)

    fig, axes = clp.initialize(1, 3, width=4.3*3, height=4.3*0.618, LaTeX=True, fontsize=12, labelthem=True, labelthemPosition=[-0.03, 1.03],  return_fig_args=True)

    # 1. Plot the photonic energy as a function of time after LP excitation
    colors0 = ['k']*n_data
    labels = ["$\\tau_{cav} = %.1f $ ps" %(t/1e3) for t in cav_lifetimes[0:-1]] + ["$\\tau_{cav} = \infty$"]

    clp.plotone(xs_UP, ys_UP, axes[0], colors=colors0, labels=labels,  lw=0.9, xlabel="time [ps]", ylabel="photonic energy [$\hbar\omega_c$]", xlim=[0, 5], ylim=[0, 0.023/au2omega0], showlegend=False)

    ax = axes[0]
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(ax.lines))]
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])

    axes[0].axvspan(0.1, 0.6, alpha=0.1, color='y')
    axes[0].text(0.2, 0.85, "excite UP", fontsize=12, transform=ax.transAxes, color='m')

    # 1.2. Add inset to plot the Rabi splitting
    axins = inset_axes(ax, width="30%", height="30%", loc=1)
    # plot the dilute CO2 IR
    freq, sp = obtain_IR(path="rerun_eq_density/box_60/E0_2e-4/")
    freq2, sp2 = obtain_IR(path="rerun_eq_density/box_60/E0_0e-4/")
    clp.plotone([freq2, freq], [sp2, sp+1.5], axins, colors=['k', 'r'], lw=1., showlegend=False, xlim=[2200, 2510], ylim=[-0.03, 4.6])
    # plot the liquid CO2 result for comparision
    freq, sp = obtain_IR(path="eq_Rabi/E0_2e-4/")
    freq2, sp2 = obtain_IR(path="eq_Rabi/E0_0e-4/")
    clp.plotone([freq2, freq], [sp2, sp+1.5], axins, colors=['k--', 'k--'], lw=1., showlegend=False, xlim=[2200, 2510], ylim=[-0.03, 4.6])
    axins.axvline(x=2320, color='b', linestyle="--")
    axins.set_xlabel("frequency [cm$^{-1}$]", fontsize=6)
    axins.set_ylabel("IR spectrum", fontsize=6)
    axins.set_xticks(ticks=[2200, 2400])
    axins.set_yticks(ticks=[])

    # 2. Plot the photonic energy as a function of time after UP excitation
    clp.plotone(xs_LP, ys_LP, axes[1], colors=colors0, labels=labels,  lw=0.9, xlabel="time [ps]", ylabel="photonic energy [$\hbar\omega_c$]", xlim=[0, 5], ylim=[0, 0.012/au2omega0], showlegend=False)
    ax = axes[1]
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(ax.lines))]
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])

    axes[1].axvspan(0.1, 0.6, alpha=0.1, color='y')
    axes[1].text(0.2, 0.85, "excite LP", fontsize=12, transform=ax.transAxes, color='c')

    # 3. Plot the polaritonic decay rate as a function of cavity loss rate
    labels = ["UP rate", "LP rate"]
    clp.plotone([cav_loss_array]*2, [1.0/taus_UP, 1.0/taus_LP], axes[2], colors=["mo", "c*"], labels=labels,
            xlabel="cavity loss rate ($1/\\tau_{cav}$) [ps$^{-1}$]", ylabel="polaritonic rate [ps$^{-1}$]",
         ylim=[0,2.4])
    # 3.1 plot the fitted linear line
    kUP_fit, slope_UP = fit_linear(cav_loss_array, 1.0/taus_UP)
    kLP_fit, slope_LP = fit_linear(cav_loss_array, 1.0/taus_LP)
    clp.plotone([cav_loss_array]*2, [kUP_fit, kLP_fit], axes[2], colors=["m--", "c--"], showlegend=False)
    # 3.2 add slope, which should equal to the molecular weight of polariton approx 0.5
    axes[2].text(0.5, 0.17, "slope %.2f" %slope_UP, fontsize=12, transform=axes[2].transAxes, color='m')
    axes[2].text(0.2, 0.5, "slope %.2f" %slope_LP, fontsize=12, transform=axes[2].transAxes, color='c')

    # add legend
    lgd = axes[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
          fancybox=False, shadow=False, ncol=5)

    clp.adjust(tight_layout=False, savefile="decay_dynamics_cavloss_dilute.pdf", includelegend=lgd)

if __name__ == '__main__':
    plot_ph_dynamics()
