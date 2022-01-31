import numpy as np
import glob
import sys
import columnplots as clp
from scipy.optimize import curve_fit
from spectral_overlap import spectral_overlap
from extract_ph_decay_rates import extract_ph_decay_rates
from matplotlib.patches import ConnectionPatch

def func(x, k, a):
    return k*x + a

def fit_linear(x, y):
    popt, pocv =  curve_fit(func, x, y)
    xfit = np.concatenate((x, np.array([0.0])))
    yfit = func(xfit, *popt)
    return xfit, yfit

factor=1200

def plot_decay_rates(Amp="6e-4"):
    # Define the plot framework
    fig, axes = clp.initialize(1, 3, width=12, height=4*0.618, fontsize=12, LaTeX=True,
        labelthem=True, labelthemPosition=[0.95, 0.95], return_fig_args=True)

    def plot_Rabi(Amp="6e-4", axes=axes, path_outcav="eq_Rabi/E0_0e-4/", showtext=True, factor=factor):
        E0_lst = ["5e-5", "1e-4", "2e-4", "3e-4", "4e-4", "5e-4"]
        E0_lst_scientific = ["$5\\times 10^{-5}$", "$1\\times 10^{-4}$", "$2\\times 10^{-4}$",
                         "$3\\times 10^{-4}$", "$4\\times 10^{-4}$","$5\\times 10^{-4}$"]
        paths_incav_eq = ["eq_Rabi/E0_%s" %E0 for E0 in E0_lst]
        paths_incav_excUP = ["pumping_Rabi/E0_%s_Exc_UP_Amp_%s/" %(E0, Amp) for E0 in E0_lst]
        paths_incav_excLP = ["pumping_Rabi/E0_%s_Exc_LP_Amp_%s/" %(E0, Amp) for E0 in E0_lst]
        LP_freqs = np.array([2302, 2280, 2241, 2207, 2177, 2150])
        UP_freqs = np.array([2350, 2375, 2428, 2485, 2548, 2616])
        Middle_freqs = (LP_freqs + UP_freqs) / 2.0
        # gather IR spectrum
        # For the first figure, we plot IR spectrum
        Js_LP, Js_UP = [], []
        Polaritons_linshape = []
        for i, path_incav_eq in enumerate(paths_incav_eq):
            freq_middle = Middle_freqs[i]
            xs, ys, J_LP, J_UP = spectral_overlap(path_incav=path_incav_eq, freq_middle=freq_middle, path_outcav=path_outcav)
            Js_LP.append(J_LP)
            Js_UP.append(J_UP)
            ys = [y + 0.3 - i*0.05 for y in ys]
            Polaritons_linshape.append(ys[1])
            clp.plotone(xs[0:2], ys[0:2], axes[0], colors=["k", "r"], lw=1, xlim=[2000, 2750],
                    xlabel="frequency [cm$^{-1}$]",
                    ylabel="normalized IR spectrum [cm]",
                    showlegend=False)
            #axes[0].fill_between(xs[0], i*0.05, ys[2], color="tab:blue")
            if showtext:
                axes[0].text(2005, 0.3 - i*0.05 + 0.01, "$\widetilde{\\varepsilon} =$ %s" %(E0_lst_scientific[i]), fontsize=7.5)
        Js_LP = np.array(Js_LP)
        Js_UP = np.array(Js_UP)

        # add text for axes[0]
        xs_polariton = [LP_freqs, UP_freqs]
        dx = xs[0][2] - xs[0][1]
        nmiddle = int((2320.0 - 2000.0)/dx)
        ys_polariton = [np.array([np.max(y[:nmiddle]) for y in Polaritons_linshape])] + [np.array([np.max(y[nmiddle:]) for y in Polaritons_linshape])]
        #clp.plotone(xs_polariton, ys_polariton, axes[0], colors=["c--", "m--"],
        #            showlegend=False, lw=0.4)
        #axes[0].text(2159, 0.07, "LP", color='c')
        #axes[0].text(2552, 0.07, "UP", color='m')

        # 1. capture the poalriton decay rates
        ks_UP = extract_ph_decay_rates(paths=paths_incav_excUP, omega_lst=[2320]*len(paths_incav_excUP))
        ks_LP = extract_ph_decay_rates(paths=paths_incav_excLP, omega_lst=[2320]*len(paths_incav_excLP))
        #ax_twin = axes[1].twinx()
        ax_twin = axes[1]
        #ys_rates = [np.concatenate((ks_LP[::-1], ks_UP))]
        #xs_J = [np.concatenate((LP_freqs[::-1], UP_freqs))]
        ys_rates = [ks_LP[::-1], ks_UP]
        xs = [LP_freqs[::-1], UP_freqs]

        x_combined_overlap = np.concatenate((LP_freqs[::-1], UP_freqs))
        y_combined_overlap = np.concatenate((Js_LP[::-1], Js_UP))

        labels = ["LP", "UP", "$\\alpha|X_{\pm}^{(B)}|^2J_{\pm}$"]
        clp.plotone(xs, ys_rates, ax_twin, colors=["co", "ms"], labels=labels, ylog=False, lw=1.2,
                bothyticks=False, markersize=7)
        ax_twin.set_ylabel("polariton decay rate [ps$^{-1}$]")
        #ax_twin.tick_params(axis='y', labelcolor="r")

        # 2. plot spectral integral versus polariton freq
        print("factor for decay rates over overlap integral is %.3E" %factor)
        clp.plotone([x_combined_overlap], [y_combined_overlap * factor], axes[1], colors=["b--*"],
                labels=labels[2:],
                ylog=False, lw=1.2,
                xlabel="polariton frequency [cm$^{-1}$]", markersize=5)

        # change legend location
        axes[1].legend(loc="center right")
        #axes[1].set_ylabel("spectral overlap integral [cm]", color="b")
        #axes[1].tick_params(axis='y', labelcolor="b")

        # add lines between two figures
        connectionstyle="arc3,rad=-0.2"
        con = ConnectionPatch(xyA=(2150, 0.1482), xyB=(2150,0.062),
                      coordsA="data", coordsB="data",
                      axesA=axes[1], axesB=axes[0], arrowstyle="<-", connectionstyle=connectionstyle, lw=0.4, color="c")
        axes[1].add_artist(con)
        con = ConnectionPatch(xyA=(2617, -0.0428), xyB=(2617,0.062),
                      coordsA="data", coordsB="data",
                      axesA=axes[1], axesB=axes[0], arrowstyle="<-", connectionstyle=connectionstyle, lw=0.4, color='m')
        axes[1].add_artist(con)
        connectionstyle="arc3,rad=0.2"
        con = ConnectionPatch(xyA=(2302, 3.05), xyB=(2302,0.327),
                      coordsA="data", coordsB="data",
                      axesA=axes[1], axesB=axes[0], arrowstyle="<-", connectionstyle=connectionstyle, lw=0.4, color="c")
        axes[1].add_artist(con)
        con = ConnectionPatch(xyA=(2350, 2.10), xyB=(2350,0.333),
                      coordsA="data", coordsB="data",
                      axesA=axes[1], axesB=axes[0], arrowstyle="<-", connectionstyle=connectionstyle, lw=0.4, color='m')
        axes[1].add_artist(con)

    plot_Rabi(Amp=Amp, axes=axes[0:2], path_outcav="eq_Rabi/E0_0e-4/", showtext=True)

    def plot_detuning(Amp=Amp, Polariton="LP", ax=None, ylim=None, xlim=None, ax_twin=None, factor=factor):
        Freq_lst = np.array([2150.0, 2200, 2250, 2300, 2325, 2350, 2400, 2450, 2500])
        LP_lst = np.array([ 2129, 2167, 2202, 2232, 2244, 2255, 2272, 2284, 2294])
        UP_lst = np.array([2379, 2387, 2401, 2419, 2431, 2444, 2476, 2512, 2552])
        Middle_lst = (UP_lst + LP_lst) / 2.0
        paths = ["pumping_detuning/Freq_%d_Exc_%s_Amp_%s/" %(freq, Polariton, Amp) for freq in Freq_lst]
        paths_incav_eq = ["eq_detuning/Freq_%d/" %freq for freq in Freq_lst]
        path_outcav = "eq_Rabi/E0_0e-4/" if Amp == "6e-4" else "pumping_Rabi/E0_0e-4_Exc_Dark_Amp_4e-3/"
        # calculate decay rates
        ks = extract_ph_decay_rates(paths=paths, omega_lst=Freq_lst)
        # calculate spectral integral
        Js_LP, Js_UP = [], []
        for i, path_incav_eq in enumerate(paths_incav_eq):
            xs, ys, J_LP, J_UP = spectral_overlap(path_incav=path_incav_eq, freq_middle=Middle_lst[i], path_outcav=path_outcav)
            Js_LP.append(J_LP)
            Js_UP.append(J_UP)
        Js_LP = np.array(Js_LP)
        Js_UP = np.array(Js_UP)

        xs = [Freq_lst - 2327]
        ys_J = [Js_LP] if Polariton == "LP" else [Js_UP]
        ys_J = [factor * y for y in ys_J]
        clp.plotone(xs, ys_J, ax, colors=["b--*"] if Polariton == "LP" else ["b--*"], showlegend=False,
                yscientific=True, lw=1.2, ylog=False, bothyticks=False, markersize=5)

        clp.plotone(xs, [ks], ax_twin, colors=["co"] if Polariton == "LP" else ["ms"], lw=1.2, bothyticks=False,
            showlegend=False, ylog=False, markersize=7, xlabel="cavity mode detuning [cm$^{-1}$]",
            ylabel="polariton decay rate [ps$^{-1}$]")
        label="$\widetilde{\\varepsilon} = 2\\times 10^{-4}$ a.u."
        if Polariton == "UP":
            ax.text(0.04, 0.9, label, fontsize=12, color="k", transform=ax.transAxes)
        if Amp == "6e-4":
            ax_twin.set_ylim(-0.2, 6.2)

    plot_detuning(Amp=Amp, Polariton="UP", ax=axes[2], ax_twin=axes[2])
    plot_detuning(Amp=Amp, Polariton="LP", ax=axes[2], ax_twin=axes[2])

    clp.adjust(tight_layout=False, savefile="decay_rates_weak.pdf")

def plot_decay_rates_vs_density():
    BoxLength_lst = np.array([24.292, 30.0, 40.0, 50.0, 60.0]) # in units of Angstrom
    number_density = 216/(BoxLength_lst / 10.0)**3 # in units of nm-3
    # Here, I plot the decay rate divided by Jpm to gain the pure intermolecular interaction effect
    paths_outcav=["eq_Rabi/E0_0e-4/", "rerun_eq_density/box_30/E0_0e-4/",
        "rerun_eq_density/box_40/E0_0e-4/", "rerun_eq_density/box_50/E0_0e-4/", "rerun_eq_density/box_60/E0_0e-4/"]
    paths_incav_eq = ["eq_Rabi/E0_2e-4/", "rerun_eq_density/box_30/E0_2e-4/",
        "rerun_eq_density/box_40/E0_2e-4/", "rerun_eq_density/box_50/E0_2e-4/", "rerun_eq_density/box_60/E0_2e-4/"]
    Js_LP, Js_UP = [], []
    for path_outcav, path_incav_eq in zip(paths_outcav, paths_incav_eq):
        xs, ys, J_LP, J_UP = spectral_overlap(path_incav=path_incav_eq, freq_middle=(2241+2428)/2.0, path_outcav=path_outcav)
        Js_LP.append(J_LP)
        Js_UP.append(J_UP)
    Js_LP = np.array(Js_LP)
    Js_UP = np.array(Js_UP)
    axes = clp.initialize(1, 2, width=8.6, height=4.3*0.618,  fontsize=12, LaTeX=True,
        labelthem=True, labelthemPosition=[0.95, 0.95])

    Amp="6e-4"
    ylim=[0, 1.4]
    for Polariton in ["LP", "UP"]:
        paths = ["pumping_Rabi/E0_2e-4_Exc_%s_Amp_%s/" %(Polariton, Amp)] \
            + ["pumping_density/Box_%d_Exc_%s_Amp_%s/" %(l, Polariton, Amp) for l in BoxLength_lst[1:]]
        ks = extract_ph_decay_rates(paths=paths, omega_lst=[2320]*len(paths))
        xfit, ks_fit = fit_linear(number_density, ks)
        color = "co" if Polariton == "LP" else "ms"
        Js = Js_LP if Polariton == "LP" else Js_UP
        clp.plotone([number_density], [ks], axes[0], colors=[color], labels=[Polariton], xlim=[0, 16], ylim=[0, 1.4], markersize=7,
                xlabel="molecular number density [nm$^{-3}$]", ylabel="$1/\\tau_{\pm}$ [ps$^{-1}$]")
        clp.plotone([number_density], [ks / Js], axes[1], colors=[color], labels=[Polariton], xlim=[0, 16], ylim=[0, 15500], markersize=7,
                xlabel="molecular number density [nm$^{-3}$]", ylabel="$1/\\tau_{\pm} J_{\pm}$ [ps$^{-1}\cdot$cm$^{-1}$]")
        #clp.plotone([xfit], [ks_fit], ax, colors=["k--"], showlegend=False, lw=1)

    #clp.plotone([number_density]*2, [Js_LP, Js_UP], ax, colors=["co", "m^"])
    clp.adjust(tight_layout=True, savefile="decay_rates_weak_density.pdf")

if __name__ == '__main__':
    plot_decay_rates()
    plot_decay_rates_vs_density()
