# CavMD simulations of vibratioal polariton relaxation experiments for liquid carbon dioxide

This is a step-by-step tutorial on how to reproduce all involved results published in:

- Li, T. E.; Nitzan, A.; Subotnik, J. E. Cavity Molecular Dynamics Simulations of Vibrational Polariton-Enhanced Molecular Nonlinear Absorption. [J. Chem. Phys. 2021, 154 (9), 094124](https://doi.org/10.1063/5.0037623).

All necessary input files and post-processing scripts are saved in the floders of the current path. NERSC users should be able to reproduce the results without modification, other linux cluster users may need to slightly modify the job submission scripts.

## 1. Equilibrium trajectories of CO<sub>2</sub>
We start from the CavMD simulations of liquid CO<sub>2</sub> under VSC. Before performing laser pumping simulations, we need to know the vibrational modes of CO<sub>2</sub> inside the cavity, which can be obtained from equilibrium CavMD simulations. Let us enter the path for this job:
<pre><code>cd CO2Experiments/CO2only_changeE0/
</code></pre>  

All sections with label **(can skip)** can be skipped without hurting the final results.

### 1.0. (can skip) Generate xyz input file for liquid CO<sub>2</sub>
<pre><code>cd generate_initial_config/
</code></pre>  
This folder contains two files: **co2.xyz** stores the geometry of a single CO<sub>2</sub> molecule (in angstrom), which can be generated by hand or using molecular structure drawing software like [*Avogadro*](https://avogadro.cc/); **create_config_packmol.inp** is the [*PACLMOL*](http://m3g.iqm.unicamp.br/packmol/home.shtml) input script to generate a box of 216 CO<sub>2</sub> molecules from a single-molecule geometry.
<pre><code>packmol < create_config_packmol.inp
</code></pre>
will generate the raw input xyz file **init_raw.xyz**. This raw input xyz file sometimes cannot be directly used as the initial geometry for a MD simulation because its energy might be too large. Therefore, we need to perform a geometry energy minimization.

### 1.1. (can skip) Minimize the xyz configuration
<pre><code>mv init_raw.xyz ../minimize_initial_config/
cd ../minimize_initial_config/
i-pi input_optimize.xml &
lmp < in.lmp
</code></pre>
We use the interface of i-pi to do this energy minimization. When the simulation is done, we use the optimized geometry as the initial geometry (**init.xyz**) for CavMD simulations.
<pre><code>tail -n 650 minimize.xc.xyz > init.xyz
sed -i "s/^648/650/" init.xyz
echo "       L  0.00000e+00  0.00000e+00  0.00000e+00" >> init.xyz
echo "       L  0.00000e+00  0.00000e+00  0.00000e+00" >> init.xyz
rm minimize* RESTART log.lammps init_raw.xyz
</code></pre>
Note that we need to append two cavity photo modes (x- and y- polarized) at the end of **init.xyz** to do CavMD simulations.

### 1.2. (can skip) Define the cavity photon controlling parameter
We now return to the path CO2Experiments/CO2only_changeE0/ for CavMD simulations.
<pre><code>mv init.xyz ../
cd ../
</code></pre>
In this folder, the **photon_params.json** file controls the parameters for cavity photons. This file can be generated by
<pre><code>python create_photon_params.py && mv test.json photon_params.json
</code></pre>
Very importantly, we need to explicitly define the partial charges of each nuclus in **photon_params.json**, which is controlled by the *charge_array* option. Please go to the [tutorial of Rabi splitting](../../tutorials/Rabi_splitting/) for more details about how to define this file.
### 1.3. Submit CavMD jobs for equilibrium trajectories

We can then submit equilibrium CavMD jobs by
<pre><code>./submit_jobs.sh
</code></pre>
If the job is killed due to exceeding the requested job time (48 h), please rerun the above command.

This command will run CavMD simulations when the light-matter coupling (E0 in **photon_params.json**, or the effective coupling strength) is zero (outside the cavity) or finite (inside the cavity). For each value of E0, it will run a 150 ps NVT simulation PLUS 40 consecutive trajectories of 20ps NVT simulations. IR spectrum should be calculated by averaging 40 NVT trajectories.

### 1.4. Post-processing
In order to obtain the Rabi splitting in the IR spectrum, we calculate the dipole auto-correlation function (DACF) of molecules in the x- and y- directions (noting tha the z-direction is decoupled from the cavity). Users can use whatever software or package they like to calculate the DACF, here we use the script similar to that in the [tutorial of Rabi splitting](../../tutorials/Rabi_splitting/) to do this job:
<pre><code>cd ../
python collect_all_data_equilibrum.py CO2only_changeE0/E0_*
</code></pre>

### 1.5. Detuning effect
We are also interested in the equilibrium Rabi splitting when the cavity mode frequency is detuned. The procedure is largely the same as above:
<pre><code>cd CO2only_changeFreq/
./submit_jobs.sh
</code></pre>
After the jobs are done, we obtain the DACF as well
<pre><code>cd ../
python collect_all_data_equilibrum.py CO2only_changeFreq/Freq*
</code></pre>

### 1.6. Plotting
Finally, we plot the Rabi splitting and detuning effect:
<pre><code>python plot_IR_CO2.py
</code></pre>
which is the first figure in the manuscript.

Note that the script plot_single_IR.py can also be used to plot the IR spectrum in a single path, e.g.,
<pre><code>python plot_single_IR.py CO2only_changeE0/E0_2e-4/
</code></pre>
will plot the IR spectrum in folder CO2only_changeE0/E0_2e-4/.

## 2. Nonequilibrium response of CO<sub>2</sub>
For nonequilibrium response of CO<sub>2</sub> after a pulse excitation, *after running the above equilibrium simulations*, we go to each folder starting with exciteCO2*, and submit jobs
<pre><code>./submit_all_excitation_jobs.sh
</code></pre>
After the simulation (or when a few trajectories are finished), in general we need to obtain the time-resolved dynamics of photonic energy, IR, and local IR:
<pre><code>python collect_all_data_tdIR.py -phIR exciteCO2*/Amp*/
</code></pre>
Then the python script will obtain these dynamics for each path in exciteCO2\*/Amp\*/ . Note that the user can also parallel this data processing procedure by running many python collect_all_data_tdIR.py -phIR path/ at the same time for different path/.

We also need to obtain the information of the C=O bond energies for each molecules:
<pre><code>python collect_all_data_tdIR.py -vco exciteCO2LP/Amp_6e-3/ exciteCO2LP_2200/Amp_6e-3/ exciteCO2LP_nocavity/Amp_6e-3/
</code></pre>

After obtaining these information, please plot the rest figures with python plotting scripts (plot_*.py).

The system temperature for different simulations can also be obtained by running (e.g.,)
<pre><code>python obtain_temperature.py exciteCO2LP/Amp_6e-3/
</code></pre>

# 3. Effects of Periodic Boundary Conditions
Finally, let us check if the above nonequilibrium response depends on periodic boundary conditions. For this part of simulations, due to the recent shutdown of NERSC Cori, the submission is run on a local cluster with PBS (not Slurm). Run the following:
<pre><code>cd PBC_check/
tar -xvf eq_number.tar.bz2
cd pumping_number/
qsub submit_diabatical_exc_polariton.sh
</code></pre>
Here, because we are interested in only the dynamics of cavity photons, we do not need to output the whole xyz configurations. Instead, in the third line of the input file input_traj.xml.bak:
<pre><code>[... atom_x(648){angstrom}, atom_x(649){angstrom}]
</code></pre>
This allows to output the positions of cavity photon mode at every time step in simu_*.out files. Note that the PBS submission file needs to be modified accordingly for different computation environments.

After the simulation, run
<pre><code>python plot_decay_rates_number.py
</code></pre>
to plot the decay rate as a function of molecular number in a single simulation cell.