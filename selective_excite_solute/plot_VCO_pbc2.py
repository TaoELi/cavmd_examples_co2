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
    filenames = glob.glob("%s/%s" %(path, pattern))#[0:10]
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def obtain_VCO_trajectories(path, pattern="simu_*VCO_td.txt", nimpurity=1, kT=0.000950037256):
    data = obtain_avg_data(path=path, pattern=pattern)
    t = data[:,0] * 1e-3
    e1 = np.mean(data[:, 1:217-nimpurity], axis=-1) - kT
    # last molecule is 13CO2
    e2 = np.mean(data[:, 217-nimpurity:], axis=-1) - kT
    au2cminv = 219474.63
    return t, e1 * au2cminv, e2 * au2cminv

def obtain_VCO_trajectories_usr_function(path, pattern="simu_*.out", kT=0.000950037256):
    data = obtain_avg_data(path=path, pattern=pattern)
    t = data[:,1]
    e1 = data[:, 3] - kT
    # last molecule is 13CO2
    e2 = data[:, 4] - kT
    au2cminv = 219474.63
    return t, e1 * au2cminv, e2 * au2cminv



def obtain_pbc_data():
    ts, e1s, e2s = [], [], []
    # original simulation cell: Nsub = 216
    t, e1, e2 = obtain_VCO_trajectories(path="pumping_LP_simulation/pump_c_0.954_freq_2320.0/Exc_2177.0_Amp_6e-3_E0_4e-4/", nimpurity=10)
    ts.append(t)
    e1s.append(e1)
    e2s.append(e2)
    # enlarged pbc cells
    Ns = [432, 864, 1728, 3456]
    for N in Ns:
        t, e1, e2 = obtain_VCO_trajectories_usr_function(path="pumping_LP_simulation_pbc/pump_N_%d_one13CO2/" %N)
        ts.append(t)
        e1s.append(e1)
        e2s.append(e2)
    return ts, e1s, e2s

def plot_total():
    fig, ax = clp.initialize(1, 1, width=4.3, LaTeX=True, fontsize=10, return_fig_args=True)

    # 1. plot VCO dynamics for impurity and also bath average
    ts, e1s, e2s = obtain_pbc_data()
    ymax = max([np.max(e1s[-1]), np.max(e2s[-1])])
    colors1 = ["lightsteelblue"]*len(ts)
    colors2 = ["tomato"]*len(ts)
    #labels = ["$N_{sub}$=%d" %x for x in [216, 432, 864, 1728, 3456]]
    labels = ["default cell"] + ["cell $\\times$%d" %x for x in [2, 4, 8, 16]]

    arrowx = 2.5

    clp.plotone(ts, e1s, ax, colors=colors1, labels=labels, xlim=[0, 20],  ylim=[ymax*3e-3,  ymax*20],
            xlabel="time [ps]", ylabel="V(C=O) - k$_B$T [cm$^{-1}$/molecule]", ylog=True, alphaspacing=0.15, lw=1)
    clp.plotone(ts, e2s, ax, colors=colors2, alphaspacing=0.15, ylog=True, showlegend=False, lw=1)

    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(ax.lines)//2)]*2
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])

    # Add an arrow to show amplification
    nsize = int(np.size(e1s[0]) * arrowx / 20.0)
    ax.annotate(s='', xy=(arrowx, e1s[1][nsize]), xytext=(arrowx,e2s[1][nsize]), arrowprops=dict(arrowstyle='<->', color=colors[1]))
    ax.text(arrowx+0.5, ymax*0.15, '%dX' %(e2s[1][nsize]/e1s[1][nsize]), color=colors[1])
    ax.text(0.06, 0.8, "$\widetilde{\\varepsilon}\propto 1/\sqrt{N_{sub}}$\nexcite LP", fontsize=10, transform=ax.transAxes)
    ax.axvspan(0.1, 0.6, color='yellow')

    #for i in range(3):
    # Put a legend to the right of the current axis
    lgd=ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #ax.legend(loc='lower right')

    ax.text(0.4, 0.06, "solvent ($^{12}$CO$_2$)", fontsize=10, transform=ax.transAxes, color='blue')
    ax.text(0.4, 0.7, "solute (10 $^{13}$CO$_2$)", fontsize=10, transform=ax.transAxes, color='blue')

    clp.adjust(savefile="VCO_pbc2.pdf", tight_layout=True, includelegend=lgd)


if __name__ == "__main__":
    plot_total()
