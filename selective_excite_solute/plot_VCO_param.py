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
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

'''
E0: incoming field amplitude in a.u.
change to
E0=6e-3 <-> I = 1.2636e12 W/cm^2 * 5E-13 s =  6.32E-1 J/cm2 = 632 mJ/cm2
'''

UseUnitAU=False
E0string_strong = "$E_0=6\\times 10^{-3}$" if  UseUnitAU else "F = 632 mJ/cm$^2$"

def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def obtain_VCO_trajectories(path, pattern="simu_*VCO_td.txt", nimpurity=1):
    data = obtain_avg_data(path=path, pattern=pattern)
    t = data[:,0] * 1e-3
    e1 = np.mean(data[:, 1:217-nimpurity], axis=-1) - 0.000950037256
    # last molecule is 13CO2
    e2 = np.mean(data[:, 217-nimpurity:], axis=-1) - 0.000950037256
    au2cminv = 219474.63
    return t, e1 * au2cminv, e2 * au2cminv

def obtain_VCO_NaCl_traj(path):
    e1_lst, e2_lst = [], []
    au2cminv = 219474.63
    for i in range(1, 41):
        print("processing traj %d" %i)
        VCO_filename = path + "/simu_%d.xc.xyz.VCO_td.txt" %i
        Vinter_filename = path + "simu_%d.xc.xyz.Vinter2NaCl_500_3001.txt" %i
        data = np.loadtxt(VCO_filename)
        t = data[:,0] * 1e-3
        VCO = data[:,1:]
        Vinter = np.loadtxt(Vinter_filename)

        Vinter = np.sum(Vinter, axis=-1)
        indx_sorted = np.argsort(Vinter) # most negative comes first
        ndivide = 2
        idx_vinter_max = indx_sorted[0:ndivide]
        idx_vinter_rest = indx_sorted[ndivide:]

        e2 = np.mean(VCO[:, idx_vinter_max], axis=-1) - 0.000950037256
        e1 = np.mean(VCO[:, idx_vinter_rest], axis=-1) - 0.000950037256
        e1_lst.append(e1 * au2cminv)
        e2_lst.append(e2 * au2cminv)

    e1_avg, e2_avg = e1_lst[0], e2_lst[0]
    for e1, e2 in zip(e1_lst[1:], e2_lst[1:]):
        e1_avg += e1
        e2_avg += e2

    return t, e1_avg / len(e1_lst), e2_avg / len(e2_lst)

def plot_total():
    fig, axes = clp.initialize(1, 3, width=4.3*2, height=4.3*0.618*0.9, LaTeX=True, fontsize=10,
                                labelthem=True, labelthemPosition=[-0.014, 1.03],  return_fig_args=True)

    # 1. plot VCO dynamics for impurity and also bath average
    t, e1, e2 = obtain_VCO_trajectories(path="pumping_LP_simulation/pump_c_0.954_freq_2320.0/Exc_2241.0_Amp_6e-3_E0_2e-4/", nimpurity=10)
    xs = [t]*2
    ys = [e1, e2]
    ymax = max([np.max(e1), np.max(e2)])
    colors = ["lightsteelblue", "tomato"]
    labels = ["solvent ($^{12}$CO$_2$)", "solute (4.63\% $^{13}$CO$_2$)"]

    arrowx = 2.5

    clp.plotone(xs, ys, axes[0], colors=colors, labels=labels, xlim=[0, 20],  ylim=[ymax*3e-3,  ymax*20],
            xlabel="time [ps]", ylabel="V(C=O) - k$_B$T [cm$^{-1}$/molecule]", alphaspacing=0.03, ylog=True)
    # Add an arrow to show amplification
    nsize = int(np.size(e1) * arrowx / 20.0)
    axes[0].annotate(s='', xy=(arrowx, e1[nsize]), xytext=(arrowx,e2[nsize]), arrowprops=dict(arrowstyle='<->'))
    axes[0].text(arrowx+0.5, ymax*0.25, '%dX' %(e2[nsize]/e1[nsize]))
    axes[0].text(0.06, 0.8, "$\widetilde{\\varepsilon}=2\\times 10^{-4}$ a.u.\nexcite LP", fontsize=10, transform=axes[0].transAxes)
    axes[0].axvspan(0.1, 0.6, color='yellow')

    # 2. plot VCO dynamics for impurity and also bath average
    t, e1, e2 = obtain_VCO_trajectories(path="pumping_LP_simulation/pump_c_0.954_freq_2320.0/Exc_2177.0_Amp_6e-3_E0_4e-4/", nimpurity=10)
    xs = [t]*2
    ys = [e1, e2]
    colors = ["lightsteelblue", "tomato"]
    labels = ["solvent ($^{12}$CO$_2$)", "solute (4.63\% $^{13}$CO$_2$)"]

    clp.plotone(xs, ys, axes[1], colors=colors, labels=labels, xlim=[0, 20],  ylim=[ymax*3e-3, ymax*20],
            xlabel="time [ps]", alphaspacing=0.03, ylog=True, showlegend=False)
    # Add an arrow to show amplification
    nsize = int(np.size(e1) * arrowx / 20.0)
    axes[1].annotate(s='', xy=(arrowx, e1[nsize]), xytext=(arrowx,e2[nsize]), arrowprops=dict(arrowstyle='<->'))
    axes[1].text(arrowx+0.5, ymax*0.3, '%dX' %(e2[nsize]/e1[nsize]))
    axes[1].text(0.06, 0.8, "$\widetilde{\\varepsilon}=4\\times 10^{-4}$ a.u.\nexcite LP", fontsize=10, transform=axes[1].transAxes)
    axes[1].axvspan(0.1, 0.6, color='yellow')

    # 3. plot VCO dynamics for impurity and also bath average
    t, e1, e2 = obtain_VCO_trajectories(path="pumping_LP_simulation/pump_c_0.954_freq_2320.0/Exc_2262.0_Amp_6e-3_E0_0e-4/", nimpurity=10)
    xs = [t]*2
    ys = [e1, e2]
    colors = ["lightsteelblue", "tomato"]
    labels = ["solvent ($^{12}$CO$_2$)", "solute (4.63\% $^{13}$CO$_2$)"]

    clp.plotone(xs, ys, axes[2], colors=colors, labels=labels, xlim=[0, 20],  ylim=[ymax*3e-3, ymax*20],
            xlabel="time [ps]", alphaspacing=0.03, ylog=True, showlegend=False)
    # Add an arrow to show amplification
    nsize = int(np.size(e1) * arrowx / 20.0)
    axes[2].annotate(s='', xy=(arrowx, e1[nsize]), xytext=(arrowx,e2[nsize]), arrowprops=dict(arrowstyle='<->'), size=10)
    axes[2].text(arrowx+0.5, ymax*0.15, '%dX' %(e2[nsize]/e1[nsize]))
    axes[2].text(0.06, 0.8, "cavity off\nexcite $^{13}$C=O asym.", fontsize=10, transform=axes[2].transAxes)
    axes[2].axvspan(0.1, 0.6, color='yellow')

    #for i in range(3):
    axes[0].legend(loc='lower right')

    clp.adjust(savefile="VCO_param.pdf", tight_layout=True)

if __name__ == "__main__":
    plot_total()
