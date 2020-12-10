#!/usr/bin/env python
# coding: utf-8

# # PSO-based calibration of PARM against 316 nM 2AT data
#
# In this notebook we will try calibrating PARM against the FRET signal data for 316 nM 2AT case depicted Fig. S3 C of Kang et al. 2019 https://pubs.acs.org/doi/abs/10.1021/acsnano.9b01993.
#
# Note: After some initial tests, we determined that a degradation rule for cytosolic Ca2+ was needed to help better capture the decay of the FRET signal. We also determined that adding some additional weighting into the cost funtion for fitting the peak value helps recover an overall better fit of the peak.

# NumPy, Pandas
import numpy as np
import pandas as pd
# PyDREAM imports
from pydream.core import run_dream
from pydream.parameters import SampledParam
from pydream.convergence import Gelman_Rubin
# Scipy Distributions
from scipy.stats import norm, uniform
# PySB ODE simulator
from pysb.simulator import ScipyOdeSimulator # http://pysb.org/download
# pydream_it
from pydream_it import DreamParam #  https://github.com/LoLab-VU/pydream_it
# Import the mode
from parm import model # path_to/PARM needs to be added to your PYTHONPATH
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

# PAR2, which is in the cell membrane
# cell membrane area factor
sa_cm = model.parameters['SAcell'].value
par2_low = np.log10(1*sa_cm) # low value for area density 1/micron^2
par2_high = np.log10(100*sa_cm) # high value at overexpression density of 3000/micron^2
print("PAR2: ", 10**par2_low, 10**par2_high)
dream_param(model.parameters['PAR2_0'], loc=par2_low, width=(par2_high-par2_low))
# G-protein
gp_low = np.log10(10000) # 1e4/cell of Brinkerhoff et al.
gp_high = np.log10(100*sa_cm) # 3000/micron^2 overexpressed density from Falkenburger et al. 2010
print("G-protein: ", 10**gp_low, 10**gp_high)
dream_param(model.parameters['Gaq_0'], loc=gp_low, width=gp_high-gp_low)
# PLC and PIP2, which are also in the cell membrane
plc_low = np.log10(3*sa_cm) # 3/micron^2 endogenous PLCB1 expression from Falkenburger et al. 2010
plc_high = np.log10(10*sa_cm) # 10/micron^2 endogenous total-PLC from Falkenburge et al. 2010
dream_param(model.parameters['PLC_0'], loc=plc_low, width=plc_high-plc_low)
dream_param(model.parameters['PIP2_0'], loc=np.log10(1000), width=2) # 1000/cell to 100000/cell
# IP3R, which is in the ER membrane
#sa_er = model.parameters['SA_ER'].value
vol_er = model.parameters['Ver'].value
dream_param(model.parameters['IP3R_0'], loc=np.log10(1000), width=2) # 1000/cell to 100000/cell
# Ca2+, which is in the ER lumen
# 400-600 microM range for ER lumen of HEK-293 cells from Foyouzi-Youssefi et al.
er_ca_low = np.log10(400*microM_to_num_per_pL*vol_er)
er_ca_high = np.log10(600*microM_to_num_per_pL*vol_er)
dream_param(model.parameters['Ca_0'], loc=er_ca_low, width=er_ca_high-er_ca_low)
# Uses SwarmParam default 2-orders of magnitude on either side of the nominal values.
#dream_param.add_all_kinetic_params(model)
v_cell = model.parameters['Vcell'].value
kd_low = np.log10(38*nM_to_num_per_pL*v_cell)
kd_high = np.log10(430*nM_to_num_per_pL*v_cell)
dream_param(model.parameters['Kd_PAR2_bind_TAT'], loc=kd_low, width=kd_high-kd_low)
# Add the kinetic parameters -- the try/except guards against Expressions which can't be calibrated.
# Set the search width to 6 orders of magnitude centered on nominal values.
kf_difflim = 1.66e-3 # 1/s*number diffusion limit of 1x10^-9 1/M*s for Vcell
kf_difflim_log10 = np.log10(kf_difflim)
diff_limited = ['kf_PAR2_bind_TAT', 'kf_PAR2_bind_Gaq', 'kf_rgs_bind_gaq', 'kf_PLC_bind_Gaq', 'kf_PLC_bind_PIP2']
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
kf_upper = np.log10(1e-4)
kf_lower = np.log10(1e-9)
dream_param(model.parameters['kf_PAR2_bind_TAT'], loc=kf_lower, width=(kf_upper-kf_lower))
dream_param['kcat_tranport_erCa'] = (uniform,2, 6)
# Remove these parameters from the calibration
dream_param -= model.parameters['k_gtp_to_gdp_auto']
dream_param -= model.parameters['k_gtp_to_gdp_rgs']
dream_param -= model.parameters['kdeg_ip3']
# Get the fancy index for the parameters being calibrated
calibrate_mask = dream_param.fancy_index(model.parameters)

param_values = np.array([parm.value for parm in model.parameters])

# Mask for initial concentration of 2AT
twoat_mask = [parm.name=='TAT_0' for parm in model.parameters]


def cost_chi2(position):
    """Variance weighted chi-squared statistic.
    This function returns a goodness of fit metric based on the
    chi-squared statistic for use with chi-squared minimization.
    The cost returned by this function is given by
        chi**2/Nd = (1/Nd) sum_{i=1 to Nd} (y_exp,i - y_sim,i)**2 / (sigma_exp,i)**2 ,
    where Nd is number of experimental data points being fitted.
    The function has a minimum value of zero corresponding to a perfect fit.
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
        y_sim = sim['FRET'][fidx]
        # Get experimental measurement and variance
        y_exp = exp_data[sig_name].values
        sigma_exp = exp_data[err_name].values
        Nd += len(y_exp)
        # Compute error between simulation and experiment
        chi2 = np.sum((y_exp - y_sim)**2 / sigma_exp**2)
        if np.isnan(chi2):
            return np.inf
        chi2_vals.append(chi2)
    total_chi2 = np.sum(chi2_vals)
    return total_chi2 / Nd


def cost_chi2_peakbias(position):
    """Variance weighted chi-squared statistic with extra peak weighting.
    This function returns a goodness of fit metric based on the
    chi-squared statistic for use with chi-squared minimization.
    The cost returned by this function is given by
        chi**2/Nd = (1/Nd) sum_{i=1 to Nd} (y_exp,i - y_sim,i)**2 / (sigma_exp,i)**2 ,
    where Nd is number of experimental data points being fitted. However, the values corresponding
    to the peaks in the experimental data sets along with the two surrounding points are weighted more heavily
    in their chi-squared contributions, which is intended to slightly bias cost mimimization
    towards recovery of the experimental activation peaks.
    The function has a minimum value of zero corresponding to a perfect fit.
    However, a value of order 1 is typically considered to be a reasonably good fit while a
    value >> 1 indicates a poor fit.
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
        # Normal chi-squared
        chi2 = np.sum((y_exp - y_sim)**2 / sigma_exp**2)
        # Now increase the relative weight of fitting the
        # peak value and the two points on either side of it.
        # This is to help improve recovery of the peak activation.
        peak_idx = np.argmax(y_exp)
        slice_lower = peak_idx - 1
        if slice_lower < 0:
            slice_lower = 0
        slice_upper = peak_idx + 2
        bias_factor = 9
        if slice_lower == peak_idx:
            # The peak is the first timepoint so there
            # are only 2 peak points to reweight.
            # With bias_factor of 9 the 2 peak points,
            # max and max+1, contribute roughly 36% of total chi2 weight.
            Nd += 2*bias_factor - 2
        else:
            # The peak is not the first timepoint so there
            # are three peak points to reweight.
            # With bias_factor of 6 the 3 peak points,
            # max, max-1, and max+1, contribute roughly 35% of total chi2 weight.
            bias_factor = 6
            Nd += 3*bias_factor - 3
        chi2 += np.sum((bias_factor-1)*(y_exp[slice_lower:slice_upper] - y_sim[slice_lower:slice_upper])**2 / sigma_exp[slice_lower:slice_upper]**2)

        if np.isnan(chi2):
            return np.inf
        chi2_vals.append(chi2)
    total_chi2 = np.sum(chi2_vals)
    return total_chi2 / Nd



def loglikelihood(position):
    return -cost_chi2_peakbias(position)



nchains = 10
niterations = 10000
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
print('At iteration: ',total_iterations,' GR = ',GR)
np.savetxt('pydream_results/'+model_name+'_GelmanRubin_iteration_'+str(total_iterations)+'.txt', GR)

old_samples = sampled_params
if np.any(GR>1.2):
    starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
    while not converged:
        total_iterations += niterations
        sampled_params, log_ps = run_dream(parameters=params, likelihood=loglikelihood,
                                           niterations=niterations, nchains=nchains, start=starts, multitry=False, gamma_levels=4,
                                           adapt_gamma=True, history_thin=1, model_name=model_name,
                                           verbose=False, restart=True)


        # Save sampling output (sampled parameter values and their corresponding logps).
        for chain in range(len(sampled_params)):
            np.save('pydream_results/'+model_name+'_sampled_params_chain_' + str(chain)+'_'+str(total_iterations), sampled_params[chain])
            np.save('pydream_results/'+model_name+'_logps_chain_' + str(chain)+'_'+str(total_iterations), log_ps[chain])

        old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
        GR = Gelman_Rubin(old_samples)
        print('At iteration: ',total_iterations,' GR = ',GR)
        np.savetxt('pydream_results/'+model_name+'_GelmanRubin_iteration_' + str(total_iterations)+'.txt', GR)

        if np.all(GR<1.2):
            converged = True
