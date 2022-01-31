import numpy as np
import columnplots as clp
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json


def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories" %(len(filenames)))
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

def gather_relaxation_trajectory(paths, opposite=False):
    ts, trajs = [], []
    if opposite:
        pattern = "simu_*.VCO_summed_td_oppo.txt"
    else:
        pattern = "simu_*.VCO_summed_td.txt"
    for path in paths:
        data = obtain_avg_data(path=path, pattern=pattern)
        t, traj = data[1:,0], data[1:, 1]
        ts.append(t*1e-3) # fs to ps
        if opposite:
            trajs.append(traj / 206.0 - 0.000950037)
        else:
            trajs.append(traj / 10.0 - 0.000950037) # averaged over 10 molecules and minus 300K in a.u.
    return ts, trajs

def prepare_ph_data(path="pumping_Rabi/E0_2e-4_Exc_UP_Amp_6e-3", omega=2320, pattern="simu_*.ph*", dtfs=2):
    filenames = glob.glob("%s/%s" %(path, pattern))
    N = len(filenames)
    print("reading %s with %d files" %(path, N))
    size = np.size(np.loadtxt(filenames[0])[:,0])
    t, phx, phy = np.zeros(size), np.zeros(size), np.zeros(size)
    ph2au = 0.5 * 1.8897259885789**2 * omega**2 * 0.0000045563352812122295**2
    vph2au = 0.5 * 1.8897259885789**2 / 41.341374575751**2
    for filename in filenames:
        data = np.loadtxt(filename)
        t_loc, qx_loc, qy_loc = data[:,0], data[:,1], data[:,2]
        vx_loc = np.gradient(qx_loc, dtfs, axis=-1, edge_order=2)
        vy_loc = np.gradient(qy_loc, dtfs, axis=-1, edge_order=2)
        t += t_loc
        phx += qx_loc**2 * ph2au + vx_loc**2 * vph2au
        phy += qy_loc**2 * ph2au + vy_loc**2 * vph2au
    return t[1:-2]/N, (phx[1:-2] + phy[1:-2])/N

paths_IR = ["co2_eq/E0_0e-4/", "co2_eq/E0_2e-4/"]
paths_VR = ["try_VER_abe/Exc_10_E0_0e-4/", "try_VER_abe/Exc_10_E0_2e-4/"]

def plot_results():
    freqs, sps = gather_IR_data(paths=paths_IR)
    ts, trajs = gather_relaxation_trajectory(paths=paths_VR)
    ts2, trajs2 = gather_relaxation_trajectory(paths=paths_VR, opposite=True)
    # photonic energy
    t, ph = prepare_ph_data(path=paths_VR[1])
    ph -= 0.000950037 * 2.0 # reduce the thermal energy

    # convert all atomic energies to hbar*omega_c, where omega_c = 2320.0 cm-1
    au2omega_c = 219474.63 / 2320.0
    trajs = [x*au2omega_c for x in trajs]
    trajs2= [x*au2omega_c for x in trajs2]
    ph *= au2omega_c

    colors = ["k", "r"]
    labels = ["outside cavity", "inside cavity"]
    axes = clp.initialize(1, 4, width=12.0, height=4.*0.618, LaTeX=True, fontsize=12, labelthem=True, labelthemPosition=[0.15, 0.95])
    sps[1] += 1.5 # inverse the outside cavity lineshape
    clp.plotone(freqs, sps, axes[0], xlim=[2000, 2600], xlabel="frequency [cm$^{-1}$]",
        ylabel="$n(\omega)\\alpha(\omega)$ [arb. units]", colors=colors, labels=labels, showlegend=False, ylim=[-0.3, 4])
    clp.plotone(ts, trajs, axes[1], xlim=[-0.5, 40], ylim=[0, 0.03*au2omega_c], ylabel="V(C=O) - $k_B$T [$\hbar\omega_c$]",
        xlabel="time [ps]", colors=colors, labels=labels, lw=1.2, showlegend=False)
    clp.plotone(ts2, trajs2, axes[2], xlim=[-0.5, 40], ylim=[0, 0.001*au2omega_c], ylabel="V(C=O) - $k_B$T [$\hbar\omega_c$]",
        xlabel="time [ps]", colors=colors, labels=labels, showlegend=False, lw=1.2)
    clp.plotone([t[1::10]*1e-3], [ph[1::10]], axes[3], colors=["r"], showlegend=False, xlim=[-0.5,40],
        ylabel="$E_{ph}$ - 2$k_B$T [$\hbar\omega_c$]", xlabel="time [ps]", ylim=[0, 0.006*au2omega_c], lw=1.2)

    axes[0].text(2010, 0.1, "cavity off",  fontsize=12, color="k")
    axes[0].text(2010, 1.6, "cavity on",  fontsize=12, color="r")
    axes[1].text(0.5, 0.8, "hot molecules", transform=axes[1].transAxes, fontsize=12, color="b")
    axes[2].text(0.4, 0.1, "thermal molecules", transform=axes[2].transAxes, fontsize=12, color="b")
    axes[3].text(0.4, 0.8, "cavity photons", transform=axes[3].transAxes, fontsize=12, color="b")
    axes[0].axvline(x = 2320, ymin=0.0, ymax=3.0, ls='--', c="b")

    clp.adjust(tight_layout=True, savefile="cavity_effect.pdf")


if __name__ == "__main__":
    plot_results()
