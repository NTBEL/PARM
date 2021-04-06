# NumPy, Pandas
import numpy as np
import pandas as pd
# PyDREAM imports
from pydream.core import run_dream
from pydream.parameters import SampledParam
from pydream.convergence import Gelman_Rubin
# Scipy Distributions
from scipy.stats import norm, uniform, halfnorm
# PySB ODE simulator
from pysb.simulator import ScipyOdeSimulator # http://pysb.org/download
# pydream_it
from pydream_it import DreamParam #  https://github.com/LoLab-VU/pydream_it
# Import the mode
from parm import classic as model # path_to/PARM needs to be added to your PYTHONPATH
# Import the conversion factors from the model
from parm import microM_to_num_per_pL, nM_to_num_per_pL, nM_2AT_to_num

# Load the Experimental Data, which is in PARM/exp_data
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
# We'll just try calibrating against all 2AT doses
data_sets = ['10', '31.6','100','316','1000','3160']
# Solver
solver = ScipyOdeSimulator(model, tspan=tspan, **solver_args)


# Use DreamParam instance to log parameters for sampling.
dream_param = DreamParam()

v_cell = model.parameters['Vcell'].value
# PAR2, which is in the cell membrane
# cell membrane area factor
sa_cm = model.parameters['SAcell'].value
 # 1/micron^2 from Falkenburger et al. 2010 (muscarinic receptor)
 # 500/micron^2 from Falkenburger et al. 2010 (muscarinic receptor)
 # 2e3/cell low value from Brinkerhoff et al. 2008
 # 2e4/cell high value from Brinkerhoff et al. 2008
par2_numbers = [1*sa_cm, 2e3, 2e4]
par2_low = np.log10(min(par2_numbers))
par2_high = np.log10(max(par2_numbers))
print("PAR2: ", 10**par2_low, 10**par2_high)
dream_param(model.parameters['PAR2_0'], loc=par2_low, width=(par2_high-par2_low))
# G-protein
# 1e4/cell of Brinkerhoff et al. 2008
# 40/micron^2 endogenous G-protein density from Falkenburger et al. 2010
gp_numbers = [1e4, 40*sa_cm]
gp_low = np.log10(min(gp_numbers))
gp_high = np.log10(max(gp_numbers))
print("G-protein: ", 10**gp_low, 10**gp_high)
dream_param(model.parameters['Gaq_0'], loc=gp_low, width=gp_high-gp_low)
# PLC and PIP2, which are also in the cell membrane
plc_low = np.log10(3*sa_cm) # 3/micron^2 endogenous PLCB1 expression from Falkenburger et al. 2010
plc_high = np.log10(10*sa_cm) # 10/micron^2 endogenous total-PLC from Falkenburge et al. 2010
dream_param(model.parameters['PLC_0'], loc=plc_low, width=plc_high-plc_low)
pip_low = np.log10(49997) # basal level of PIP2 as per Lemon et al. 2003
pip_high = np.log10(5000*sa_cm) # free PIP2 of 5000 per micrometer^2 used by Falkenburger et al. 2013
dream_param(model.parameters['PIP2_0'], loc=pip_low, width=pip_high-pip_low) # 1000/cell to 1000000/cell
# IP3R, which is in the ER membrane
v_erm = model.parameters['Verm'].value
ip3r_low = np.log10(1*nM_to_num_per_pL*(v_erm/v_cell)) # 1 nM low-end for signaling molecule range from Albeck et al.
ip3r_high = np.log10(1*microM_to_num_per_pL*(v_erm/v_cell)) # 1 microM high-end for signaling molecule range from Albeck et al.
dream_param(model.parameters['IP3R_0'], loc=ip3r_low, width=ip3r_high-ip3r_low)
# Ca2+, which is in the ER lumen
# 400-600 microM range for ER lumen of HEK-293 cells from Foyouzi-Youssefi et al.
# More generally, 100-1000 microM.
vol_er = model.parameters['Ver'].value
er_ca_low = np.log10(100*microM_to_num_per_pL*vol_er)
er_ca_high = np.log10(1000*microM_to_num_per_pL*vol_er)
dream_param(model.parameters['Ca_0'], loc=er_ca_low, width=er_ca_high-er_ca_low)
# Add the baseline cytosolic concentration for range 10-150 nM
cyt_ca_low = np.log10(10*nM_to_num_per_pL*v_cell)
cyt_ca_high = np.log10(150*nM_to_num_per_pL*v_cell)
dream_param(model.parameters['Ca_C_0'], loc=cyt_ca_low, width=cyt_ca_high-cyt_ca_low)
# Set the Kd for agonist binding to PAR2
# PAR2 agonists in HEK 293T cells - LeSarge et al. https://doi.org/10.1021/acsmedchemlett.9b00094
kd_low = np.log10(38*nM_to_num_per_pL*v_cell) # Isox-Cha-Chg-ARK(Sulfo-Cy5)-NH2 has Kd = 38 nM with EC50 = 16 nM
kd_high = np.log10(430*nM_to_num_per_pL*v_cell) # 2f-LIGRLO(Sulfo-Cy5)-NH2 has Kd = 430 nM with EC50 = 296 nM
# Since 2AT has EC50 = 101.8 nM in Hek 293 cells we'll
# assume that the Kd for 2AT is somewhere between those two compounds.
dream_param(model.parameters['Kd_PAR2_bind_TAT'], loc=kd_low, width=kd_high-kd_low)
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
                k_low = np.log10(nominal) - 2
                k_upper = np.log10(nominal) + 2
                if 10**k_upper > kf_difflim_log10: kf_upper = kf_difflim_log10
                dream_param(rule.rate_forward, loc=k_low, width=(k_upper - k_low))
            else:
                dream_param(rule.rate_forward, width=6)
        except:
            pass
    if rule.rate_reverse:
        try:
            dream_param(rule.rate_reverse, width=6)
        except:
            pass
# Reset a few kinetic parameters with more specific intervals
v_extra = model.parameters['Vextra'].value
kf_upper = np.log10(kf_difflim*(v_cell/v_extra)) # upper diffusion limit adjusted for extracellular volume
kf_lower = np.log10(1e-9) # lower value of 1e-9 based on previous PSO tests.
dream_param(model.parameters['kf_PAR2_bind_TAT'], loc=kf_lower, width=(kf_upper-kf_lower))

dream_param['kcat_tranport_erCa'] = (uniform,2, 6) # 1e2-1e8 1/s
# Remove these parameters from the calibration
dream_param -= model.parameters['k_gtp_to_gdp_auto']
dream_param -= model.parameters['k_gtp_to_gdp_rgs']
dream_param -= model.parameters['kdeg_ip3']
# Get the fancy index for the parameters being calibrated
calibrate_mask = dream_param.fancy_index(model.parameters)
print('parms to calibrate: ',calibrate_mask)
param_values = np.array([parm.value for parm in model.parameters])

# Mask for initial concentration of 2AT
twoat_mask = [parm.name=='TAT_0' for parm in model.parameters]


like_fret = dict()
for data_set in data_sets:
    sig_name = data_set + "_sig"
    err_name = data_set + "_err"
    # Get experimental measurement and std. dev.
    y_exp = exp_data[sig_name].values[1:]
    sigma_exp = exp_data[err_name].values[1:]
    like_fret[data_set] = norm(loc=y_exp, scale=sigma_exp)

like_fret_slope = dict()
for data_set in data_sets:
    sig_name = data_set + "_sig"
    err_name = data_set + "_err"
    # Get experimental measurement and std. dev.
    y_exp = exp_data[sig_name].values
    sigma_exp = exp_data[err_name].values
    fd_exp = y_exp[1:] - y_exp[:-1]
    fd_exp_sigma = np.sqrt(sigma_exp[1:]**2 + sigma_exp[:-1]**2)
    like_fret_slope[data_set] = norm(loc=fd_exp, scale=fd_exp_sigma)

def loglikelihood(position):
    """log(Likelihood) function.
    This function includes components for fitting the training data at each 2AT
    concentration, as well as added weight for fitting the peak in the experimental
    curve, additional contribution from fitting the slope of the curve, and
    a contribution accouting for upward shifting in time of the peak at 330 nM 2AT
    after 50% of the PAR2 is inactivated.
    """
    Y = 10**position
    # Start with first soft constraint (data fit) data set so that if
    # the parameter set is really bad (ODE integrator problems)
    # we can exit having only tried to run one simulation.
    param_values[calibrate_mask] = Y
    logp_vals = list()
    for data_set in data_sets:
        sig_name = data_set + "_sig"
        err_name = data_set + "_err"
        tat_conc = float(data_set)*nM_2AT_to_num
        param_values[twoat_mask] = tat_conc
        sim = solver.run(param_values=param_values).all
        # Get the simulated FRET response
        y_sim = sim['FRET'][fidx][1:]
        logp_s = like_fret[data_set].logpdf(y_sim)
        logp = np.sum(logp_s)
        # Compute error between simulation and experiment
        if np.isnan(logp):
            return np.inf
        logp_vals.append(logp)
        # Add extra check for fitting the slope (forward difference)
        y_sim = sim['FRET'][fidx]
        fd_sim = y_sim[1:] - y_sim[:-1]
        logp = np.sum(like_fret_slope[data_set].logpdf(fd_sim))
        logp_vals.append(logp)
        # Now increase the relative weight of fitting the
        # peak value in the exp. data.
        y_exp = exp_data[sig_name].values[1:]
        peak_idx = np.argmax(y_exp)
        slice_lower = peak_idx - 1
        if slice_lower < 0:
            slice_lower = 0
        slice_upper = peak_idx + 2
        if slice_lower == peak_idx:
            # The peak is the first timepoint >0 so there
            # are only 2 peak points to reweight.
            bias_factor = (len(y_exp) - 2) / 2 # with this bias factor the two peak points
            # should contribute as much weight as fitting all other points.
        else:
            # The peak is not the first timepoint so there
            # are three peak points to reweight.
            bias_factor = (len(y_exp) - 3) / 3# with this bias factor the three peak points
            # should contribute as much weight as fitting all other points.

        logp = bias_factor * np.sum(logp_s[slice_lower:slice_upper])
        #print(logp)
        logp_vals.append(logp)
    # Constraint for timepoint of peak FRET ratio.
    # Most data suggest the peak doesn't shift up in time much after MH inactivation of PAR2.
    # In most cases where a shift is discernable it is only 1 measurement unit (9 s).
    # We'll incorporate this via a halfnorm.
    param_values[calibrate_mask] = Y
    tat_conc = 330
    tat_conc *= nM_2AT_to_num
    param_values[twoat_mask] = tat_conc
    #fracs = [0.5, 0.015, 0.02, 0.05]
    frac = 0.5 # inactivation of 50% of PAR2
    sim = solver.run(param_values=param_values).all
    # Get the simulated FRET response
    fret_sim = sim['FRET'][fidx]
    ref_peak_idx = np.argmax(fret_sim)
    ref_peak_val = fret_sim[ref_peak_idx]
    ref_peak_time = tspan[fidx][ref_peak_idx]
    param_values[calibrate_mask[0]] = Y[0]*frac
    sim = solver.run(param_values=param_values).all
    # Get the simulated FRET response
    fret_sim = sim['FRET'][fidx]
    peak_idx = np.argmax(fret_sim)
    peak_fret = fret_sim[peak_idx]
    peak_time = tspan[fidx][peak_idx]
    dt = peak_time - ref_peak_time
    # halfnorm centered at zero (no shift upwards) with scale 9 (i.e., the measurement interval of 9 s)
    logp = halfnorm.logpdf(dt, scale=9)
    logp_vals.append(logp)
    # cost = -log(Likelihood)
    return np.sum(logp_vals)

nchains = 4
niterations = 50000
params = dream_param.sampled_params()

print(len(params))

# Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
converged = False
total_iterations = niterations
print("About to run_dream...")
model_name = 'parm_dreamzs_{}chains'.format(nchains)
sampled_params, log_ps = run_dream(parameters=params, likelihood=loglikelihood,
                                   niterations=niterations, nchains=nchains, multitry=False,
                                   gamma_levels=4, adapt_gamma=True, history_thin=1,
                                   model_name=model_name, verbose=False)

# Save sampling output (sampled parameter values and their corresponding logps).
for chain in range(len(sampled_params)):
    np.save('pydream_results/'+model_name+'_sampled_params_chain_' + str(chain)+'_'+str(total_iterations), sampled_params[chain])
    np.save('pydream_results/'+model_name+'_logps_chain_' + str(chain)+'_'+str(total_iterations), log_ps[chain])

#Check convergence and continue sampling if not converged
# Using Gelman-Rubin statistic to monitor convergence
GR = Gelman_Rubin(sampled_params)
np.savetxt('pydream_results/'+model_name+'_GelmanRubin_iteration_'+str(total_iterations)+'.txt', GR)
print('At iteration: ',total_iterations,' GR = ',GR)
print('GR>1.2 = : ',GR>1.2)
