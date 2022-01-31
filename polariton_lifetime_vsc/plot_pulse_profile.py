import numpy as np
import glob
import sys
import columnplots as clp

def get_data(pattern="eq_Rabi/E0_2e-4/*dac.txt"):
    files = glob.glob(pattern)
    data = np.loadtxt(files[0])
    for i in range(1, len(files)):
        data += np.loadtxt(files[i])
    data /= len(files)
    return data

def get_Rabi_splitting():
    data = get_data()
    freq, sp = data[:,5], (data[:,6]+data[:,7])/1.2e29
    # we do truncation
    freq_start = 2000
    freq_end = 3000
    df = freq[1] - freq[0]
    nstart = int(freq_start//df)
    nend = int(freq_end//df)
    return freq[nstart:nend], sp[nstart:nend]

#Plot everything in atomic units
E0_au = 6e-4
cminv2au = 1.0 / 219474.63
omega_LP = 2241 * cminv2au
omega_UP = 2428 * cminv2au
fs2au = 41.341374575751

def get_Eprofile(t, E0, omega, t0, t1):
    E = E0 * np.cos(omega*t)
    E[t<t0] = 0.0
    E[t>t1] = 0.0
    return E

def get_spectrum(t, signal):
    sp = np.fft.fft(signal)
    dt = t[1] - t[0]
    freq = np.fft.fftfreq(np.size(t), d=dt)
    freq *= 2.0 * np.pi
    return freq, np.abs(sp)

t = np.linspace(0, 2e5, 10000)
t_start = 100 * fs2au
t_end = 600 * fs2au
E_LP = get_Eprofile(t, E0_au, omega_LP, t_start, t_end)

freq, sp_LP = get_spectrum(t, E_LP)
freq_cminv = freq / cminv2au

freq_Rabi, sp_Rabi = get_Rabi_splitting()

axes = clp.initialize(2, 1, width=4.3,   height=4.3*0.618*2, fontsize=12, LaTeX=True,
            labelthem=True, labelthemPosition=[0.12, 0.95])

clp.plotone([t/fs2au], [E_LP], axes[0], colors=["k"], showlegend=False, xlim=[0, 700],
    xlabel="time [fs]", ylabel="$E_x(t)$ [a.u.]", lw=1)

clp.plotone([freq_cminv], [sp_LP], axes[1], colors=["k"], showlegend=False,
    xlim=[2100, 2600],
    xlabel="frequency [cm$^{-1}$]", ylabel="$|E_x(\omega)|$ [a.u.]", lw=1)

ax2 = axes[1].twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'
ax2.set_ylabel('IR spectrum [arb. units]', color=color)  # we already handled the x-label with ax1
ax2.plot(freq_Rabi, sp_Rabi, color=color, lw=1)
ax2.tick_params(axis='y', labelcolor=color)

clp.adjust(tight_layout=True, savefile="pulse_profile.pdf")
