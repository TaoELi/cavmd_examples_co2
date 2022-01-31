import numpy as np
import columnplots as clp
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json
import matplotlib.cm as cm
import matplotlib.patches as patches
from matplotlib import gridspec
import matplotlib.pyplot as plt


nlst = [1]
def obtain_avg_data(path="pumping_LP_N_216/E0_2e-4_Amp_6e-3/"):
    n = 500
    filenames = [path+"/simu_%d.xc.xyz.localIR_statistics_%d_%d.txt" %(k, n, n+2501) for k in nlst]
    print(filenames)
    data_lst = []
    for filename in filenames:
        data = np.loadtxt(filename)
        data_lst.append(data)
    return data_lst

def prepare_data(path="pumping_LP_N_216/E0_2e-4_Amp_6e-3/"):
    data_lst = obtain_avg_data(path=path)
    sp_lst = []
    for data in data_lst:
        freq, sp = data[:, 0], data[:, 1:]
        x = freq
        dx = x[2] - x[1]
        nstart, nend = int(2000 / dx), int(2500 / dx)
        x = x[nstart:nend]
        sp = np.abs(sp[nstart:nend,:] / 2e28)
        sp = sp[::-1]
        sp_lst.append(sp)
    return x, sp_lst

def obtain_IR_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))#[0:2]
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def obtain_IR(path, pattern="simu_*.dac.txt"):
    data = obtain_IR_avg_data(path=path, pattern=pattern)
    freq = data[:,5]
    sp = (data[:,6] + data[:,7]) / 2e28
    return freq, sp

def plot_IR():
    fig, gs = clp.initialize_gridSpec(4, 20, width=4.3*2.0, height=16.0/3.0, LaTeX=True, fontsize=10)
    ax1 = fig.add_subplot(gs[0, 0:5])
    ax2 = fig.add_subplot(gs[0, 7:])
    ax3 = fig.add_subplot(gs[1, 0:5])
    ax4 = fig.add_subplot(gs[1, 7:])
    ax5 = fig.add_subplot(gs[2, 0:5])
    ax6 = fig.add_subplot(gs[2, 7:])
    ax7 = fig.add_subplot(gs[3, 0:5])
    ax8 = fig.add_subplot(gs[3, 7:])
    axes = [ax2, ax4, ax6, ax8, ax1, ax3, ax5, ax7]

    # Plot IR spectrum at equilibrium
    # 2. Plot inset, the equilibrium IR spectroscopy
    E0_lst = ["2e-4", "2e-4", "4e-4", "4e-4"]
    for i in range(4):
        # obtain CO2 IR spectrocopy
        path_outcav = "equilibrium_simulation/eq_c_1.0_freq_2320.0/E0_0e-4/"
        path_incav = "equilibrium_simulation/eq_c_1.0_freq_2320.0/E0_%s/" %E0_lst[i]
        freq, sp1 = obtain_IR(path_outcav)
        freq, sp2 = obtain_IR(path_incav)
        sp2 += 1.7
        xs, ys = [freq]*2, [sp1, sp2]
        axins = axes[4+i]
        clp.plotone(xs, ys, axins, colors=["k", "r"], showlegend=False, lw=1.0,
                    xlim=[2050, 2600], ylim=[-0.04, 3.8])
        #axins.set_xlabel("IR freq [cm$^{-1}$]")
        #axins.set_ylabel("Intensity [arb. units]")

        if i == 0:
            axins.text(0.57, 0.1, "$^{12}$C=O", fontsize=9, transform=axins.transAxes)
            axins.text(0.18, 0.55, "LP", fontsize=9, transform=axins.transAxes, color='r')
            axins.text(0.79, 0.55, "UP", fontsize=9, transform=axins.transAxes, color='r')
        if i <= 2:    
            axins.set_xticklabels([])

        if i == 0:
            str = "$\widetilde{\\varepsilon}= 2\\times 10^{-4}$"
        elif i == 1:
            str = "$\widetilde{\\varepsilon}= 2\\times 10^{-4}$"
        elif i == 2:
            str = "$\widetilde{\\varepsilon}= 4\\times 10^{-4}$"
        else:
            str = "$\widetilde{\\varepsilon}= 4\\times 10^{-4}$"
        #axins.text(0.02, 0.85, str, fontsize=8, transform=axins.transAxes)
        axins.axvline(x=2320.0, color='blue', linestyle='--', alpha=0.5, lw=0.9)

    # Add 13CO2 frequency
    path_outcav_13CO2 = "equilibrium_simulation/eq_c_0.0_freq_2320.0/E0_0e-4/"
    freq, sp1 = obtain_IR(path_outcav_13CO2)
    clp.plotone([freq], [sp1], axes[5], colors=["0.50"], showlegend=False, lw=1.0)
    clp.plotone([freq], [sp1], axes[7], colors=["0.50"], showlegend=False, lw=1.0)
    axes[5].text(0.05, 0.1, "$^{13}$C=O\n(4.63\%)", fontsize=9, transform=axes[5].transAxes, color="0.50")

    ax7.text(0.5, -0.33, "IR freq [cm$^{-1}$]", ha='center', transform=ax7.transAxes)
    ax8.text(0.5, -0.33, "No. molecule", ha='center', transform=ax8.transAxes)

    fig.text(0.07, 0.5, "Intensity [arb. units]",  va='center', rotation='vertical')
    fig.text(0.35, 0.5, "local IR freq [cm$^{-1}$]",  va='center', rotation='vertical')

    fontsize=7
    x0, y0 = -0.02, 1.1
    ax1.text(x0, y0, "(a)", transform=ax1.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    ax2.text(x0, y0, "(b)", transform=ax2.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    ax3.text(x0, y0, "(c)", transform=ax3.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    ax4.text(x0, y0, "(d)", transform=ax4.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    ax5.text(x0, y0, "(e)", transform=ax5.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    ax6.text(x0, y0, "(f)", transform=ax6.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    ax7.text(x0, y0, "(g)", transform=ax7.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    ax8.text(x0, y0, "(h)", transform=ax8.transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')

    freq, sp_lst = prepare_data("pumping_LP_simulation/pump_c_1.0_freq_2320.0/Exc_2241.0_Amp_6e-3_E0_2e-4/")
    freq2, sp_lst2 = prepare_data("pumping_LP_simulation/pump_c_0.954_freq_2320.0/Exc_2241.0_Amp_6e-3_E0_2e-4/")
    freq3, sp_lst3 = prepare_data("pumping_LP_simulation/pump_c_1.0_freq_2320.0/Exc_2177.0_Amp_6e-3_E0_4e-4/")
    freq4, sp_lst4 = prepare_data("pumping_LP_simulation/pump_c_0.954_freq_2320.0/Exc_2177.0_Amp_6e-3_E0_4e-4/")
    extent = [0, 215, freq[0] , freq[-1]]
    from matplotlib.colors import LogNorm, Normalize
    vmax1 = np.max(np.max(sp_lst[0])) * 1.0
    vmin1 = vmax1 * 1e-3
    N=4

    # Plot local IR nonequilibrium
    im1 = axes[0].imshow(sp_lst[0], aspect='auto', extent=extent,  cmap=cm.hot, interpolation=None, norm=LogNorm(vmin=vmin1, vmax=vmax1))
    im2 = axes[1].imshow(sp_lst2[0], aspect='auto', extent=extent,  cmap=cm.hot, interpolation=None, norm=LogNorm(vmin=vmin1, vmax=vmax1))
    im3 = axes[2].imshow(sp_lst3[0], aspect='auto', extent=extent,  cmap=cm.hot, interpolation=None, norm=LogNorm(vmin=vmin1, vmax=vmax1))
    im4 = axes[3].imshow(sp_lst4[0], aspect='auto', extent=extent,  cmap=cm.hot, interpolation=None, norm=LogNorm(vmin=vmin1, vmax=vmax1))

    x = range(0, 215)
    xs = [x]
    ys = [np.ones(len(x)) * 2200]
    for i in range(N):
        clp.plotone(xs, ys, axes[i], showlegend=False, colors=["c--"], lw=1., alpha=0.)

    for i in range(N-1):
        axes[i].set_xticklabels([])

    arrow_properties = dict(facecolor="r", width=0.5, headwidth=4, shrink=0.1)
    excitations = [2241, 2241, 2177, 2177]
    E0_lst = ["pure $^{12}$CO$_2$, $\widetilde{\\varepsilon}= 2\\times 10^{-4}$ a.u.",
              "4.63\% $^{13}$CO$_2$@$^{12}$CO$_2$, $\widetilde{\\varepsilon}= 2\\times 10^{-4}$ a.u.",
              "pure $^{12}$CO$_2$, $\widetilde{\\varepsilon}= 4\\times 10^{-4}$ a.u.",
              "4.63\% $^{13}$CO$_2$@$^{12}$CO$_2$, $\widetilde{\\varepsilon}= 4\\times 10^{-4}$ a.u."]

    taus_LP =[0.2, 0.3, 2.9, 0.6]
    for i in range(N):
        axes[i].tick_params(color='c', labelsize='medium', width=2)
        axes[i].annotate("", xy=(0, excitations[i]), xytext=(-6, excitations[i]), arrowprops=arrow_properties)
        axes[i].text(5, 2388, E0_lst[i] + ", single traj 1-6 ps", color='w', fontsize=10)
        axes[i].text(5, 2035, "LP lifetime %.1f ps" %taus_LP[i], color='w', fontsize=10)
    axes[1].axvline(x=205.0, color='white', linestyle='--', lw=1.0)
    axes[3].axvline(x=205.0, color='white', linestyle='--', lw=1.0)

    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.96, 0.15, 0.02, 0.7])
    fig.colorbar(im1, cax=cbar_ax)

    clp.adjust(savefile="excitation_statistics_solvent.pdf")

if __name__ == "__main__":
    plot_IR()
