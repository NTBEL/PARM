# NumPy, Pandas, and Matplotlib
import numpy as np
import pandas as pd

from simplepso.pso import PSO


from pysb.simulator import ScipyOdeSimulator # http://pysb.org/download
from swarm_it import SwarmIt, SwarmParam #  https://github.com/LoLab-VU/swarm_it
from parm import classic as model # path_to/PARM needs to be added to your PYTHONPATH
from parm import microM_to_num_per_pL, nM_to_num_per_pL, nM_2AT_to_num

# Experimental data is in PARM/exp_data
exp_data = pd.read_csv('../../exp_data/FRET_data_Kang2019_FigS3C.csv')

times = exp_data['Time'].values

def expand_times(times, expand_by=100):
    tsp = np.linspace(0,times.max(), expand_by)
    tsp = np.concatenate((tsp, times[1:-1]))
    tsp.sort()
    findex = list()
    for time in times:
        idx = np.where(tsp == time)[0][0]
        findex.append(idx)
    return tsp, findex

tspan, fidx = expand_times(times)

# PySB solver options.
integrator_options = {"rtol": 1e-6, "atol": 1e-6}
solver_args = {'integrator': 'lsoda', 'integrator_options':integrator_options}
solver = ScipyOdeSimulator(model, tspan=tspan, **solver_args)

data_sets = ['10', '31.6', '100', '316', '1000', '3160']

def fancy_index(swarm_param, parameters):
    idxs = list()
    for name in swarm_param.names():
        i = 0
        for parm in parameters:
            if parm.name == name:
                idxs.append(i)
                break
            i += 1
    return idxs

# Use SwarmParam instance to log parameters for sampling.
swarm_param = SwarmParam()
# Cell volume
v_cell = model.parameters['Vcell'].value
# cell membrane area factor
sa_cm = model.parameters['SAcell'].value
# PAR2, which is in the cell membrane
 # 1/micron^2 from Falkenburger et al. 2010 (muscarinic receptor)
 # 500/micron^2 from Falkenburger et al. 2010 (muscarinic receptor)
 # 2e3/cell low value from Brinkerhoff et al. 2008
 # 2e4/cell high value from Brinkerhoff et al. 2008
par2_numbers = [1*sa_cm, 500*sa_cm, 2e3, 2e4]
par2_low = np.log10(min(par2_numbers))
par2_high = np.log10(max(par2_numbers))
print("PAR2: ", 10**par2_low, 10**par2_high)
swarm_param(model.parameters['PAR2_0'], loc=par2_low, width=(par2_high-par2_low))
# G-protein
# 1e4/cell of Brinkerhoff et al. 2008
# 40/micron^2 endogenous G-protein density from Falkenburger et al. 2010
gp_numbers = [1e4, 40*sa_cm]
gp_low = np.log10(min(gp_numbers))
gp_high = np.log10(max(gp_numbers))
print("G-protein: ", 10**gp_low, 10**gp_high)
swarm_param(model.parameters['Gaq_0'], loc=gp_low, width=gp_high-gp_low)
# PLC and PIP2, which are also in the cell membrane
plc_low = np.log10(3*sa_cm) # 3/micron^2 endogenous PLCB1 expression from Falkenburger et al. 2010
plc_high = np.log10(10*sa_cm) # 10/micron^2 endogenous total-PLC from Falkenburger et al. 2010
swarm_param(model.parameters['PLC_0'], loc=plc_low, width=plc_high-plc_low)
pip_low = np.log10(49997) # basal level of PIP2 as per Lemon et al. 2003
pip_high = np.log10(5000*sa_cm) # free PIP2 of 5000 per micrometer^2 used by Falkenburger et al. 2013
swarm_param(model.parameters['PIP2_0'], loc=pip_low, width=pip_high-pip_low)
# IP3R, which is in the ER membrane
v_erm = model.parameters['Verm'].value
ip3r_low = np.log10(1*nM_to_num_per_pL*(v_erm/v_cell)) # 1 nM low-end for signaling molecule range from Albeck et al.
ip3r_high = np.log10(1*microM_to_num_per_pL*(v_erm/v_cell)) # 1 microM high-end for signaling molecule range from Albeck et al.
swarm_param(model.parameters['IP3R_0'], loc=ip3r_low, width=ip3r_high-ip3r_low)
# Ca2+ in the ER lumen
# 400-600 microM range for ER lumen of HEK-293 cells from Foyouzi-Youssefi et al.
# Broader, less specific, range of 100-1000 microM
vol_er = model.parameters['Ver'].value
er_ca_low = np.log10(100*microM_to_num_per_pL*vol_er)
er_ca_high = np.log10(1000*microM_to_num_per_pL*vol_er)
swarm_param(model.parameters['Ca_0'], loc=er_ca_low, width=er_ca_high-er_ca_low)
# Add the baseline cytosolic concentration for range 10-150 nM

cyt_ca_low = np.log10(10*nM_to_num_per_pL*v_cell)
cyt_ca_high = np.log10(150*nM_to_num_per_pL*v_cell)
swarm_param(model.parameters['Ca_C_0'], loc=cyt_ca_low, width=cyt_ca_high-cyt_ca_low)
# Set the Kd for agonist binding to PAR2
# PAR2 agonists in HEK 293T cells - LeSarge et al. https://doi.org/10.1021/acsmedchemlett.9b00094
kd_low = np.log10(38*nM_to_num_per_pL*v_cell) # Isox-Cha-Chg-ARK(Sulfo-Cy5)-NH2 has Kd = 38 nM with EC50 = 16 nM
kd_high = np.log10(430*nM_to_num_per_pL*v_cell) # 2f-LIGRLO(Sulfo-Cy5)-NH2 has Kd = 430 nM with EC50 = 296 nM
# Since 2AT has EC50 = 101.8 nM in Hek 293 cells we'll
# assume that the Kd for 2AT is somewhere between those two compounds.
swarm_param(model.parameters['Kd_PAR2_bind_TAT'], loc=kd_low, width=kd_high-kd_low)
# Add the kinetic parameters -- the try/except guards against Expressions which can't be calibrated.
# Set the search width to 6 orders of magnitude centered on nominal values.
kf_difflim = 1.66e-3 # 1/s*number diffusion limit of 1x10^-9 1/M*s for Vcell
kf_difflim_log10 = np.log10(kf_difflim)
diff_limited = ['kf_PAR2_bind_Gaq', 'kf_rgs_bind_gaq', 'kf_PLC_bind_Gaq', 'kf_PLC_bind_PIP2']
for rule in model.rules:
    if rule.rate_forward:
        try:
            if rule.rate_forward.name in diff_limited:
                nominal = rule.rate_forward.value
                k_low = np.log10(nominal) - 5
                k_upper = np.log10(nominal) + 5
                if k_upper > kf_difflim_log10: kf_upper = kf_difflim_log10
                swarm_param(rule.rate_forward, loc=k_low, width=(k_upper - k_low))
            else:
                swarm_param(rule.rate_forward, width=8)
        except:
            pass
    if rule.rate_reverse:
        try:
            swarm_param(rule.rate_reverse, width=8)
        except:
            pass
# Reset a few kinetic parameters with more specific intervals
v_extra = model.parameters['Vextra'].value
kf_upper = np.log10(kf_difflim*(v_cell/v_extra)) # upper diffusion limit adjusted for extracellular volume
kf_lower = np.log10(1e-9) # lower value of 1e-9 based on previous PSO tests.
swarm_param(model.parameters['kf_PAR2_bind_TAT'], loc=kf_lower, width=(kf_upper-kf_lower))

swarm_param['kcat_tranport_erCa'] = (2, 8) # 1e2-1e8 1/s

# Remove these parameters from the calibration
swarm_param -= model.parameters['k_gtp_to_gdp_auto']
swarm_param -= model.parameters['k_gtp_to_gdp_rgs']
swarm_param -= model.parameters['kdeg_ip3']

calibrate_mask = fancy_index(swarm_param, model.parameters)
print(calibrate_mask)

param_values = np.array([parm.value for parm in model.parameters])
# Mask for initial concentration of 2AT
twoat_mask = [parm.name=='TAT_0' for parm in model.parameters]

def cost_chi2_forwarddiff_peakbias(position):
    """Variance weighted chi-squared statistic.
    This function returns a goodness of fit metric based on the
    chi-squared statistic for use with chi-squared minimization.
    The cost returned by this function is given by
        chi**2/Nd = (1/Nd) sum_{i=1 to Nd} (y_exp,i - y_sim,i)**2 / (sigma_exp,i)**2 ,
    where Nd is number of experimental data points being fitted.
    A function value of zero corresponds to a perfect fit.
    """
    Y = 10**position
    param_values[calibrate_mask] = Y
    chi2_vals = list()
    Nd = 0
    for data_set in data_sets:
        sig_name = data_set + "_sig"
        err_name = data_set + "_err"
        tat_conc = float(data_set)*nM_2AT_to_num
        param_values[twoat_mask] = tat_conc
        sim = solver.run(param_values=param_values).all
        # Get the simulated FRET response
        y_sim = sim['FRET'][fidx][1:]
        # Get experimental measurement and variance
        y_exp = exp_data[sig_name].values[1:]
        sigma_exp = exp_data[err_name].values[1:]
        Nd += len(y_exp)
        # Compute error between simulation and experiment
        chi2 = np.sum((y_exp - y_sim)**2 / sigma_exp**2)
        # Now increase the relative weight of fitting the
        # peak value and the two points on either side of it.
        # This is to help improve recovery of the peak activation.
        peak_idx = np.argmax(y_exp)
        slice_lower = peak_idx - 1
        if slice_lower < 0:
            slice_lower = 0
        slice_upper = peak_idx + 2

        if slice_lower == peak_idx:
            # The peak is the first timepoint so there
            # are only 2 peak points to reweight.
            # With bias_factor of 9 the 2 peak points,
            # max and max+1, contribute roughly 36% of total chi2 weight.
            bias_factor = 6
            Nd += 2*bias_factor - 2
        else:
            # The peak is not the first timepoint so there
            # are three peak points to reweight.
            # With bias_factor of 6 the 3 peak points,
            # max, max-1, and max+1, contribute roughly 35% of total chi2 weight.
            bias_factor = 3
            Nd += 3*bias_factor - 3
        chi2 += np.sum((bias_factor-1)*(y_exp[slice_lower:slice_upper] - y_sim[slice_lower:slice_upper])**2 / sigma_exp[slice_lower:slice_upper]**2)
        y_exp = exp_data[sig_name].values
        sigma_exp = exp_data[err_name].values
        fd_exp = y_exp[1:] - y_exp[:-1]
        fd_exp_sigma = np.sqrt(sigma_exp[1:]**2 + sigma_exp[:-1]**2)/2
        y_sim = sim['FRET'][fidx]
        fd_sim = y_sim[1:] - y_sim[:-1]
        chi2 += np.sum((fd_exp - fd_sim)**2 / fd_exp_sigma**2)
        Nd += len(fd_exp)
        if np.isnan(chi2):
            return np.inf
        chi2_vals.append(chi2)
    total_chi2 = np.sum(chi2_vals)
    avg_total_chi2 = total_chi2 / Nd
    #if avg_total_chi2 > 4.:
     #   return np.inf
    return avg_total_chi2

# Setup the particle swarm optimization run

# Set the number of particles in the swarm.
num_particles = 160
# Set the number of iterations for PSO run.
num_iterations = 30
# Number of PyDREAM chains that we want to generate starting points for.
nchains = 8
# Number of processors for PSO multiprocessing.
nproc = 8
print("Doing PSO {} times with {} particles and {} iterations".format(nchains,num_particles, num_iterations))
print("Using {} processors".format(nproc))
costs = list()
for i in range(nchains):
    print("Running PSO for chain {}".format(i))
    # Construct the optimizer
    pso = PSO(save_sampled=False,
              verbose=True)
    starting_position = swarm_param.centers()
    pso.set_start_position(starting_position)
    pso.set_bounds(lower=swarm_param.lower(), upper=swarm_param.upper())
    # sets maximum speed that a particle can travel
    pso.set_speed(-0.5, 0.5)
    # run it
    pso.run(num_particles,
            num_iterations, cost_function=cost_chi2_forwarddiff_peakbias,
            stop_threshold=1e-3,
            num_processors=nprocs)
    best = pso.best.pos
    costs.append(pso.best.fitness)
    print("best at chain {} with cost {} :".format(i, pso.best.fitness))
    print(best)
    np.save("pso_best_chain_{}".format(i),best)
    np.savetxt("pso_best_costs.csv",np.array(costs),delimiter=',')
