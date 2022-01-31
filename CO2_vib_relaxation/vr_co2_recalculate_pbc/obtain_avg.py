'''
Combine all output to a single output file
'''
import numpy as np
import sys
import glob
import os

def calc_ph_energy(phx, phy, omega=2320, dtfs=2):
    ph2au = 0.5 * 1.8897259885789**2 * omega**2 * 0.0000045563352812122295**2
    vph2au = 0.5 * 1.8897259885789**2 / 41.341374575751**2
    qx_loc, qy_loc = phx, phy
    vx_loc = np.gradient(qx_loc, dtfs, axis=-1, edge_order=2)
    vy_loc = np.gradient(qy_loc, dtfs, axis=-1, edge_order=2)
    e_x = qx_loc**2 * ph2au + vx_loc**2 * vph2au
    e_y = qy_loc**2 * ph2au + vy_loc**2 * vph2au
    return e_x + e_y

def obtain_avg(path):
    print("working with", path)
    # data_raw save time [ps], temperature, VCO_hot, photonic energy
    data_raw = np.loadtxt(path+"/simu_1.out")
    data_raw = np.zeros((np.shape(data_raw)[0], 4))
    files = glob.glob(path+"/simu_*.out")
    Ntraj = 0
    for fname in files:
        print(fname, os.path.getsize(fname))
        if os.path.getsize(fname) >= 3620730:
            data = np.loadtxt(fname)
            # time
            data_raw[:, 0] += data[:, 1]
            # temperature
            data_raw[:, 1] += data[:, 2]
            # VCO_hot
            data_raw[:, 2] += data[:, 3]
            # photonic energy
            data_raw[:, 3] += calc_ph_energy(data[:, -2], data[:, -6])
            Ntraj += 1
    data_raw /= Ntraj
    print("%d trajs saved" %Ntraj)
    np.save(path+"/simu_avg.out", data_raw)

if __name__ == '__main__':
    try:
        obtain_avg(path=sys.argv[-1])
    except:
        print("Error occrured!!!")
