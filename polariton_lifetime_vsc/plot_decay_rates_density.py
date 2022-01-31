import numpy as np
import glob
import sys
import columnplots as clp
from scipy.optimize import curve_fit
from spectral_overlap import spectral_overlap
from extract_ph_decay_rates import extract_ph_decay_rates

'''
Under weak excitations
Fig. a plots IR spectrum under different effective coupling strengths
Fig. b plots (i) spectral overlap integral and (ii) polariton decay rates versus Rabi splitting
'''

def func(x, k, a):
    return k*x + a

def fit_linear(x, y):
    popt, pocv =  curve_fit(func, x, y)
    xfit = np.concatenate((x, np.array([0.0])))
    yfit = func(xfit, *popt)
    return xfit, yfit

def plot_decay_rates_density():
    BoxLength_lst = np.array([24.292, 30.0, 40.0, 50.0, 60.0])
    number_density = 216/BoxLength_lst**3

    def plot_sub(Amp="6e-3", Polariton="LP", ax=None, ylim=[0, 1.1]):
        paths = ["pumping_Rabi/E0_2e-4_Exc_%s_Amp_%s/" %(Polariton, Amp)] \
                + ["pumping_density/Box_%d_Exc_%s_Amp_%s/" %(l, Polariton, Amp) for l in BoxLength_lst[1:]]
        ks = extract_ph_decay_rates(paths=paths, omega_lst=[2320]*len(paths))
        xfit, ks_fit = fit_linear(number_density, ks)
        if Amp == "6e-3":
            label="$F = 632$ mJ/cm$^{2}$ \n excite %s" %Polariton
        else:
            label="$F = 6.32$ mJ/cm$^{2}$ \n excite %s" %Polariton
        color = "c" if Polariton == "LP" else "m"
        clp.plotone([number_density], [ks], ax, colors=[color+"o"], showlegend=False, xlim=[0, 0.016], ylim=ylim)
        clp.plotone([xfit], [ks_fit], ax, colors=["k--"], showlegend=False, lw=1)
        ax.text(0.3, 0.65, label, transform=ax.transAxes, color=color)

    # Define the plot framework
    axes = clp.initialize(2, 2, width=5.3,  fontsize=12, LaTeX=True,
            labelthem=True, labelthemPosition=[0.15, 0.95], sharex=True,
            commonX=[0.5, -0.05, "number density [$\AA^{-3}$]"],
            commonY=[-0.05, 0.5, "polariton decay rate [ps$^{-1}$]"])

    plot_sub(Amp="6e-4", Polariton="LP", ax=axes[0,0], ylim=[0, 1.1])
    plot_sub(Amp="6e-4", Polariton="UP", ax=axes[0,1], ylim=[0, 0.4])
    plot_sub(Amp="6e-3", Polariton="LP", ax=axes[1,0], ylim=[0, 10])
    plot_sub(Amp="6e-3", Polariton="UP", ax=axes[1,1], ylim = [0, 0.4])

    clp.adjust(tight_layout=True, savefile="decay_rates_density.pdf")


if __name__ == '__main__':
    plot_decay_rates_density()
