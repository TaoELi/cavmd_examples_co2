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
from sklearn.linear_model import LinearRegression
from matplotlib.ticker import PercentFormatter
from matplotlib import animation
from matplotlib.patches import Rectangle

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

# parameters for the third plot
def obtain_tot_data_hist(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    tot_data = []
    for filename in filenames:
        data = np.loadtxt(filename)
        tot_data.append(data)
    tot_data = np.array(tot_data)
    print(np.shape(tot_data))
    return tot_data

def plot_main():

    axes = clp.initialize(1, 2, width=4.3*2, height=4.3*0.618, LaTeX=True, fontsize=12, labelthem=True, labelthemPosition=[-0.1, 1.05])

    # The third spectrum
    def plot_hist(paths, ax, adding_on_outcav = 0.0, ylabel=True):
        data_tot = []
        for i, path in enumerate(paths):
            data = obtain_tot_data_hist(path=path, pattern="simu_*.xc.xyz.VCO_td_oppo.txt") / 0.010603
            data_tot.append(data)

        #labels = ["inside cavity", "outside cavity"]
        nhist_lst, nbins_lst = [], []
        for h in [30]:
            for i, path in enumerate(paths):
                #data = np.reshape(data_tot[i][:,h,1:], -1)
                data = np.reshape(data_tot[i][:,:,1:], -1)
                print(np.shape(data))
                if i == 0:
                    ax.hist(data, range=(0, 2), bins=30, density=True, log=True, color = "orchid", ec="black")
                else:
                    ax.hist(data + adding_on_outcav, range=(0, 2), bins=30, density=True, log=True, color = "c", ec="black", alpha=0.5)
                ax.set_ylim(1e-3, 10)
                if ylabel:
                    ax.set_ylabel("Density distribution")
                #ax.text(0.3, 0.8, labels[i], fontsize=12, transform=ax.transAxes)
                nhist, bins = np.histogram(data, bins=30, range=(0,2), density=True)
                nhist_lst.append(np.array(nhist))
                nbins_lst.append(np.array(bins[:-1] + (bins[1] - bins[0]) ))
        return nhist_lst, nbins_lst

    paths2 = ["vr_co2_pbc/vr_N_2000/Exc_accordingly_E0_6.5727e-5", "vr_co2_pbc/vr_N_2000/Exc_accordingly_E0_0e-4"]

    handles = [Rectangle((0,0),1,1,color="orchid",ec="k"), Rectangle((0,0),1,1,color="c",ec="k", alpha=0.5)]
    labels= ["inside cavity","outside cavity"]
    axes[0].legend(handles, labels)

    nhist_lst, nbins_lst = plot_hist(paths2, axes[0], adding_on_outcav=0.0)

    axes[0].set_xlabel("C=O bond potential energy in thermal CO$_2$ [$\hbar\omega_0$]")
    axes[0].set_xlim(0, 2.0)
    label = "$N_{sub}=2000$\n$N_{hot}/N_{sub}$ = 4.63\%"
    axes[0].text(0.5, 0.35, label, fontsize=12, transform=axes[0].transAxes)

    xs = [nbins_lst[0]]
    ys = [(nhist_lst[0] / nhist_lst[1])]
    # Add additional plot to show only the difference
    clp.plotone(xs, ys, axes[1], colors=["k-o"], labels=["cavity on / cavity off"], ylabel="Density ratio", xlim=[0, 1.5])
    axes[1].set_xlabel("C=O bond potential energy in thermal CO$_2$ [$\hbar\omega_0$]")

    clp.adjust(tight_layout=True, savefile="nonlinear.pdf")

if __name__ == "__main__":
    plot_main()
