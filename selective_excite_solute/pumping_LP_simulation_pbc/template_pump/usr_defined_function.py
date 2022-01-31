import numpy as np

Nsub=216
Nimpurity=10

def VCO(x):
    x -= 2.196
    return 0.45106005 * x**2 - 0.67297933 * x**3 + 0.5857145 * x**4

def usr_defined_function(p, q, m):
    # q[0::3] --> x component of molecules
    # q[1::3] --> y component of molecules
    # q[2::3] --> z component of molecules
    qx, qy, qz = q[0::3], q[1::3], q[2::3]
    # get rid of the cavity photon modes
    qx, qy, qz = qx[:-2], qy[:-2], qz[:-2]
    # total number of atoms
    natoms = np.size(qx)
    # Let us assume that the molecules are arranged in the following order
    # C; O; O; C; O; O...
    # calculate the CO bond length in atomic units
    CO_1 = np.sqrt((qx[0::3] - qx[1::3])**2 + (qy[0::3] - qy[1::3])**2 + (qz[0::3] - qz[1::3])**2)
    CO_2 = np.sqrt((qx[0::3] - qx[2::3])**2 + (qy[0::3] - qy[2::3])**2 + (qz[0::3] - qz[2::3])**2)
    # calculate the C=O bond energy
    VCO_smoothed = VCO(CO_1) + VCO(CO_2)
    # sum over the selected molecules
    VCO_solvent = np.mean(VCO_smoothed[0:Nsub-Nimpurity])
    VCO_impurity = np.mean(VCO_smoothed[Nsub-Nimpurity:])
    
    return np.array([VCO_solvent, VCO_impurity])
