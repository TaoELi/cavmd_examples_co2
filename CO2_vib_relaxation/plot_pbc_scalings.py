'''
Containing 3 subplots given a fixed N_hot/N_sub ratio:
(a) C=O bond energy versus N_sub
(b) Polaritonic spectrum versus N_sub
(c) VCO bond energy density distribution
'''
import numpy as np
import columnplots as clp
import glob
from scipy import signal
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# parameters for the first plot
Ns = np.array([216, 302, 324, 432, 648, 864, 1080, 1296, 1728, 2160, 2592, 3456, 4320, 5184, 6048, 6912, 8000, 12960])

paths_incav_percentage = ["vr_%d/Exc_%d_E0_*.*/" %(N, round(10.0/216.0*N)) for N in Ns]
paths_outcav_percentage = ["vr_%d/Exc_%d_E0_0e-4/" %(N, round(10.0/216.0*N)) for N in Ns]

paths_incav_10 = ["vr_%d/Exc_%d_E0_*.*/" %(N, 10) for N in Ns]
paths_outcav_10 = ["vr_%d/Exc_%d_E0_0e-4/" %(N, 10) for N in Ns]

paths_incav_increasinghot = ["vr_2160/Exc_%d_E0_*.*/" %N for N in [10,20,30,40,50,60,70,80,90,100]]
paths_outcav_increasinghot = ["vr_2160/Exc_%d_E0_0e-4/" %N for N in [10,20,30,40,50,60,70,80,90,100]]

def exp_func(x, k, a):
    return a*np.exp(-x*k)

def fit_exp(x, y):
    popt, pocv = curve_fit(exp_func, x, y)
    perr = np.sqrt(np.diag(pocv))
    print("error is", perr)
    return popt[0], exp_func(x, *popt), perr[0]

def lin_func(x, k, a):
    return k*x + a

def fit_linear(x, y):
    popt, pocv = curve_fit(lin_func, x, y)
    x_new = np.array([0.0] + x.tolist())
    #x_new = x
    return x_new, lin_func(x_new, *popt), popt

def get_data_statistics(path="E0_2e-4/", pattern="simu_avg.out.npy"):
    fname = glob.glob(path + pattern)[0]
    data = np.load(fname)
    t = data[1::10,0]
    VCO = data[1::10,2]
    return t, VCO

def get_results(paths_incav, N_excs, dir="vr_co2_recalculate_pbc/"):
    xs, ys = [], []
    ks_fit, ys_fit = [], []
    ks_fit_upper, ks_fit_lower = [], []
    for i, path in enumerate(paths_incav):
        N_MODIFIED = N_excs[i]
        print("Path = %s, modifying %d molecules" %(path, N_MODIFIED))
        t, data_avg = get_data_statistics(path=dir+path)
        xs.append(t)
        KT = 0.000950037
        ys.append(data_avg - KT) # minus k_BT, the thermal energy
    for t, y in zip(xs, ys):
        print("doing exp fit")
        k, y_fit, k_error = fit_exp(t, y)
        ks_fit.append(k)
        ys_fit.append(y_fit)
        ks_fit_upper.append(k + k_error)
        ks_fit_lower.append(k - k_error)
    ks_fit = np.array(ks_fit)
    ks_fit_upper = np.array(ks_fit_upper)
    ks_fit_lower = np.array(ks_fit_lower)
    return xs, ys, ks_fit, ys_fit, ks_fit_upper, ks_fit_lower

# parameters for the second plot
def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

Ns_subset = np.array([216, 864, 4320])
Nexcs_subset = np.array([10, 30, 90])

paths_incav_percentage_subset = ["vr_%d/Exc_%d_E0_*.*/" %(N, round(10.0/216.0*N)) for N in Ns_subset]
paths_incav_10_subset = ["vr_%d/Exc_%d_E0_*.*/" %(N, 10) for N in Ns_subset]
paths_incav_increasinghot_subset = ["vr_2160/Exc_%d_E0_*.*/" %N for N in Nexcs_subset]

def get_energy(path="E0_2e-4/", pattern="simu_avg.out.npy"):
    fname = glob.glob(path + pattern)[0]
    data = np.load(fname)
    t, ph = data[1:-2:10,0], data[1:-2:10,-1]
    # delete thermal energy
    KT = 0.000950037
    return t, ph - 2.0 * KT

def get_spectrum(path="E0_2e-4/", pattern="ph_sp_avg.npy"):
    fname = glob.glob(path + pattern)[0]
    data = np.load(fname)
    freq, sp = data[:,0], data[:,1]
    sp /= 1e32
    # smooth the spectrum
    sp = smooth(sp)
    df = freq[1] - freq[0]
    nstart, nend = int(2150 // df), int(2500 // df)
    freq_sub = freq[nstart:nend]
    sp_sub = sp[nstart:nend]
    integration = np.sum(sp_sub)*df
    print("integral is", integration)
    return freq_sub, sp_sub, integration

def get_results_sp(paths_incav, dir="vr_co2_recalculate_pbc/"):
    xs_sp, ys_sp = [], []
    for i, path in enumerate(paths_incav):
        freq, sp, __ = get_spectrum(path=dir+path)
        xs_sp.append(freq)
        ys_sp.append(sp)
    return xs_sp, ys_sp

def get_sp_area(paths_incav, dir="vr_co2_recalculate_pbc/"):
    areas = []
    for path in paths_incav:
        __, __, integration = get_spectrum(path=dir+path)
        areas.append(integration)
    areas = np.array(areas)
    return areas


def plot_main():

    axes = clp.initialize(2, 2, width=6.3, height=4.3*0.618*2, LaTeX=True, fontsize=12, labelthem=True, labelthemPosition=[0.15, 0.95])
    axes = [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]

    # The first plot
    def plot_rate_versus_Nsub(indx, N_excs, x, paths_incav, paths_outcav):
        # plot for 10 molecules are excited
        xs, ys_in, ks_fit, ys_fit, ks_fit_upper, ks_fit_lower = get_results(paths_incav, N_excs=Nexcs)
        x_new, yfit, popt = fit_linear(x, ks_fit)
        yfit_10 = yfit
        print("When N --> inf, the rate approaches %.5E inside the cavity and slope is %.5E" %(yfit[0], popt[0]))
        xs, ys = [x, x_new], [ks_fit, yfit]
        clp.plotone(xs, ys, axes[indx], colors=["ro", "r--"],
            showlegend=False, xlim=[0, np.max(x_new)*1.1],
            ylim=[0, np.max(ks_fit)*1.1],
            alpha=0.5, xlog=False, ylog=False)
        #clp.plotone([x, x], [ks_fit_upper, ks_fit_lower], axes[indx], colors=["r", "b"], showlegend=False)

        xs, ys_out, ks_fit2, ys_fit2, ks_fit_upper2, ks_fit_lower2 = get_results(paths_outcav, N_excs=Nexcs)
        ys_diff = [y2 - y1 for y1, y2 in zip(ys_in, ys_out)]
        x_new, yfit2, popt = fit_linear(x, ks_fit2)
        print("When N --> inf, the rate approaches %.5E outside the cavity and slope is %.5E" %(yfit[0], popt[0]))
        clp.plotone([x, x_new], [ks_fit2, yfit2], axes[indx], colors=["ks", "k--"], showlegend=False, alpha=0.5, xlog=False, ylog=False)
        #clp.plotone([x, x], [ks_fit_upper2, ks_fit_lower2], axes[indx], colors=["r", "b"], showlegend=False, alpha=0.5)

        ks_fit3 = ks_fit - ks_fit2
        x_new, yfit3, popt = fit_linear(x, ks_fit3)
        print(popt)
        #clp.plotone([x], [ks_fit3], axes[indx], colors=["g*", "g--"], showlegend=False)

        clp.plotone([x, x_new], [ks_fit3, yfit3], axes[indx], colors=["g*", "g--"], showlegend=False, alpha=0.5, xlog=False, ylog=False)
        return ks_fit3

    Nexcs = np.array([10]*len(paths_incav_10))
    ks_fit_part1 = plot_rate_versus_Nsub(indx=0, N_excs=Nexcs, x=1.0/Ns, paths_incav=paths_incav_10, paths_outcav=paths_outcav_10)
    axes[0].text(0.3, 0.85, "$N_{hot}$ = 10, change $N_{sub}$", transform=axes[0].transAxes, fontsize=10,)
    axes[0].set_xlabel("1/$N_{sub}$")
    axes[0].set_ylabel("fitted decay rate [ps$^{-1}$]")
    axes[0].set_ylim(0, 0.13)

    Nexcs = np.array([round(10.0/216.0*N) for N in Ns])
    ks_fit_part3 = plot_rate_versus_Nsub(indx=2, N_excs=Nexcs, x=1.0/Ns**0.7, paths_incav=paths_incav_percentage, paths_outcav=paths_outcav_percentage)
    axes[2].text(0.22, 0.85, "$N_{hot} / N_{sub}$ = 4.63\% (10/216)", transform=axes[2].transAxes, fontsize=10,)
    axes[2].set_xlabel("$1/N_{sub}^{0.7}$")
    axes[2].set_ylabel("fitted decay rate [ps$^{-1}$]")
    axes[2].set_ylim(0, 0.13)

    axes[0].text(0.17, 0.54, "inside cavity", fontsize=10, color='r', transform=axes[0].transAxes, alpha=0.5)
    axes[0].text(0.3, 0.05, "outside cavity", fontsize=10, color='k', transform=axes[0].transAxes, alpha=0.5)
    axes[0].text(0.65, 0.24, "difference\n(cavity effect)", fontsize=10, color='g', transform=axes[0].transAxes)

    # The second plot
    def plot_sp(indx, x, paths_incav, labels):
        xs_sp, ys_sp = get_results_sp(paths_incav)
        clp.plotone(xs_sp, ys_sp, axes[indx],  labels=labels, xlabel="frequency [cm$^{-1}$]", xlim=[2150, 2500], ylim=[0, 0.34], lw=1.2, alphaspacing=0.02)
        colormap = plt.cm.hot
        colors = [colormap(i) for i in np.linspace(0.0, 0.6, len(axes[indx].lines))]
        leg = axes[indx].get_legend()
        for i,j in enumerate(axes[indx].lines):
            j.set_color(colors[i])
            leg.legendHandles[i].set_color(colors[i])

    def plot_sp_inset(indx, x, paths_incav, xlabel):
        areas = get_sp_area(paths_incav)
        axins = axes[indx].inset_axes([0.27, 0.55, 0.3, 0.4])
        clp.plotone([x], [areas], axins, colors=["ko"], showlegend=False)
        axins.set_xlabel(xlabel, fontsize=8)
        ylabel="area [arb. units]"
        axins.set_ylabel(ylabel, fontsize=8)
        for tick in axins.xaxis.get_major_ticks():
            tick.label.set_fontsize(8)
            # specify integer or one of preset strings, e.g.
            tick.label.set_fontsize('x-small')
            #tick.label.set_rotation('vertical')
        for tick in axins.yaxis.get_major_ticks():
            tick.label.set_fontsize(6)
            # specify integer or one of preset strings, e.g.
            #tick.label.set_fontsize('x-small')
            #tick.label.set_rotation('vertical')

    labels = ["$N_{sub}$ = %d" %N for N in Ns_subset]
    plot_sp(indx=1, x=1.0/Ns_subset, paths_incav=paths_incav_10_subset, labels=labels)
    axes[1].legend(fontsize = 'x-small')
    axes[1].set_ylabel("photonic spectrum [arb. units]")
    plot_sp_inset(indx=1, x=Ns**(0.5), paths_incav=paths_incav_10, xlabel="$\sqrt{N_{sub}}$")

    labels = ["$N_{sub}$ = %d" %N for N in Ns_subset]
    plot_sp(indx=3,  x=1.0/Ns_subset, paths_incav=paths_incav_percentage_subset, labels=labels)
    #axes[1].set_title("$N_{hot} / N_{sub}$ = 4.63\% (10/216)")
    axes[3].legend(fontsize = 'x-small')
    axes[3].set_ylabel("photonic spectrum [arb. units]")
    plot_sp_inset(indx=3, x=Ns**(0.5), paths_incav=paths_incav_percentage, xlabel="$\sqrt{N_{sub}}$")


    clp.adjust(tight_layout=True, savefile="pbc_dependence_scaling.pdf")

if __name__ == "__main__":
    plot_main()
