import random

# NumPy, Pandas, and Matplotlib
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import seaborn as sns

sns.set(rc={"figure.dpi": 100, "savefig.dpi": 300})
sns.set_palette("bright")
sns.set_context("notebook")
sns.set_style("ticks")

from scipy.stats import norm, halfnorm

from simplepso.pso import PSO

from pysb.simulator import ScipyOdeSimulator  # http://pysb.org/download
from swarm_it import SwarmIt, SwarmParam  #  https://github.com/LoLab-VU/swarm_it
from parm.parm import (
    model,
)  # path_to/PARM needs to be added to your PYTHONPATH
from parm import units
from parm import data
from parm import util


# Experimental data is in PARM/exp_data
# Load the Experimental Data, which is in PARM/exp_data
exp_data = (
    data.training_data()
)  # dynamic FRET response at different 2AT doses (Fig. S3c)
times = exp_data["Time"].values  # the time values
exp_data_2 = (
    data.kang2019_fig_s3d()
)  # the peak FRET dose-response curve data (Fig. S3d)

times = exp_data["Time"].values
print("Exp. time points:")
print(times)
tspan, fidx = util.expand_times(times)
print("Simulation time span:")
print(tspan)
# PySB solver options.
integrator_options = {"rtol": 1e-6, "atol": 1e-6}
solver_args = {"integrator": "lsoda", "integrator_options": integrator_options}
solver = ScipyOdeSimulator(model, tspan=tspan, **solver_args)

data_sets_all = ["10", "31.6", "100", "316", "1000", "3160"]
# data_sets = ["10", "316", "1000", "3160"]
# Reverse so it from high to low 2AT
data_sets = list(reversed(data_sets_all))


def fancy_index(swarm_param, parameters):
    idxs = list()
    for name in swarm_param.names():
        i = 0
        for param in parameters:
            if param.name == name:
                idxs.append(i)
                break
            i += 1
    return idxs


# Get the nominal parameter values.
param_values = np.array([param.value for param in model.parameters])

# Use SwarmParam instance to log parameters for sampling.
swarm_param = SwarmParam()
# Cell volume
v_cell = model.parameters["Vcyto"].value
# cell membrane area factor
sa_cm = model.parameters["SAcell"].value
# PAR2, which is in the cell membrane
# 1/micron^2 from Falkenburger et al. 2010 (muscarinic receptor)
# 500/micron^2 from Falkenburger et al. 2010 (muscarinic receptor)
# 2e3/cell low value from Brinkerhoff et al. 2008
# 2e4/cell high value from Brinkerhoff et al. 2008
# AuNP/cell in plate saturates at 160/cell (w/ estimated 20-30% cell loss)
# AuNP/cell in suspension gives ~1000/cell (w/ precise cell counting)
# 14364 = PAR2/cell for 647 antibody + 3*sigma
# 2000 ~ AuNP/cell adjusted for ~8% PAR2 positive cells.
par2_numbers = [2000, 3000, 14364]
par2_low = np.log10(min(par2_numbers))
par2_high = np.log10(max(par2_numbers))
print("PAR2: ", par2_low, 10 ** par2_low, par2_high, 10 ** par2_high)
swarm_param(model.parameters["PAR2_0"], loc=par2_low, width=(par2_high - par2_low))
# G-protein
# 1e4/cell of Brinkerhoff et al. 2008, also Yi et al. PNAS 2003
# 40/micron^2 endogenous G-protein density from Falkenburger et al. 2010
gp_numbers = [1000, 1e4, 40 * sa_cm]
gp_low = np.log10(min(gp_numbers))
gp_high = np.log10(max(gp_numbers))
print("G-protein: ", gp_low, 10 ** gp_low, gp_high, 10 ** gp_high)
swarm_param(model.parameters["Gprotein_0"], loc=gp_low, width=gp_high - gp_low)
# PLC and PIP2, which are also in the cell membrane
plc_low = np.log10(
    3 * sa_cm
)  # 3/micron^2 endogenous PLCB1 expression from Falkenburger et al. 2010
plc_high = np.log10(
    10 * sa_cm
)  # 10/micron^2 endogenous total-PLC from Falkenburger et al. 2010
# swarm_param(model.parameters["PLC_0"], loc=plc_low, width=plc_high - plc_low)
# print("PLC: ", plc_low, 10 ** plc_low, plc_high, 10 ** plc_high)
pip_low = np.log10(49997)  # basal level of PIP2 as per Lemon et al. 2003
pip_high = np.log10(
    5000 * sa_cm
)  # free PIP2 of 5000 per micrometer^2 used by Falkenburger et al. 2013
# swarm_param(model.parameters["PIP2_0"], loc=pip_low, width=pip_high - pip_low)
# print("PIP2: ", pip_low, 10 ** pip_low, pip_high, 10 ** pip_high)
# IP3R, which is in the ER membrane
v_erm = model.parameters["Verm"].value
sa_erm = model.parameters["SAer"].value
ip3r_low = np.log10(100)  # 1 nM low-end for signaling molecule range from Albeck et al.
ip3r_high = np.log10(
    1e4
)  # 1 microM high-end for signaling molecule range from Albeck et al.
swarm_param(model.parameters["IP3R_0"], loc=ip3r_low, width=ip3r_high - ip3r_low)
print("IP3R: ", ip3r_low, 10 ** ip3r_low, ip3r_high, 10 ** ip3r_high)
# Ca2+ in the ER lumen
# 400-600 microM range for ER lumen of HEK-293 cells from Foyouzi-Youssefi et al.
# with average around 530 +- 70 uM.
# Broader, less specific, range of 100-1000 microM
vol_er = model.parameters["Ver"].value
er_ca_low = np.log10(
    110 * units.microM_to_molec_per_pL * vol_er
)  # ~ -3sigma Foyouzi-Youssefi
er_ca_high = np.log10(
    950 * units.microM_to_molec_per_pL * vol_er
)  # ~ +3sigma Foyouzi-Youssefi
# swarm_param(model.parameters["Ca_E_0"], loc=er_ca_low, width=er_ca_high - er_ca_low)
# print("Ca_E: ", er_ca_low, 10 ** er_ca_low, er_ca_high, 10 ** er_ca_high)
# Add the baseline cytosolic concentration for range 10-150 nM
# 97 +- 5 nM from Tong et al. JBC 1999
# cyt_ca_low = np.log10(82*units.nM_to_molec_per_pL*v_cell) # -3sigma
# cyt_ca_high = np.log10(113*units.nM_to_molec_per_pL*v_cell) # +3sigma
# swarm_param(model.parameters['Ca_C_0'], loc=cyt_ca_low, width=cyt_ca_high-cyt_ca_low)

# Set the Kd for agonist binding to PAR2
# PAR2 agonists in HEK 293T cells - LeSarge et al. https://doi.org/10.1021/acsmedchemlett.9b00094
v_extra = model.parameters["Vextra"].value
kd_low = np.log10(
    38 * units.nM_to_molec_per_pL * v_extra
)  # Isox-Cha-Chg-ARK(Sulfo-Cy5)-NH2 has Kd = 38 nM with EC50 = 16 nM
kd_high = np.log10(
    430 * units.nM_to_molec_per_pL * v_extra
)  # 2f-LIGRLO(Sulfo-Cy5)-NH2 has Kd = 430 nM with EC50 = 296 nM
# Since 2AT has EC50 = 101.8 nM in Hek 293 cells we'll
# assume that the Kd for 2AT is somewhere between those two compounds.
swarm_param(model.parameters["Kd_PAR2_bind_TAT"], loc=kd_low, width=kd_high - kd_low)
print("Kd_PAR2_bind_TAT: ", kd_low, 10 ** kd_low, kd_high, 10 ** kd_high)
# Protein-protein association rate constants -- including 2AT+PAR2 also
ppf_params = [
    "kf_PAR2_bind_TAT",
    "kf_PAR2_activate_Gprotein",
    "kf_irrev_heterotrimer_reassociation",
    "kf_PLC_bind_Gaq",
]
# Other forward binding rate constants -- we'll assume these are constrained to the
# same ranges as the protein-protein interactions.
bkf_params = ["kf_PLC_bind_PIP2"]

kf_low = model.parameters["kf_PAR2_bind_TAT"].value / 100
kf_high = model.parameters["kf_PAR2_bind_TAT"].value * 100
swarm_param(
    model.parameters["kf_PAR2_bind_TAT"],
    loc=np.log10(kf_low),
    width=(np.log10(kf_high) - np.log10(kf_low)),
)
print("kf_PAR2_bind_TAT: ", np.log10(kf_low), kf_low, np.log10(kf_high), kf_high)
kf_low = model.parameters["kf_PAR2_activate_Gprotein"].value / 100
kf_high = model.parameters["kf_PAR2_activate_Gprotein"].value * 100
swarm_param(
    model.parameters["kf_PAR2_activate_Gprotein"],
    loc=np.log10(kf_low),
    width=(np.log10(kf_high) - np.log10(kf_low)),
)
print(
    "kf_PAR2_activate_Gprotein: ", np.log10(kf_low), kf_low, np.log10(kf_high), kf_high
)
kf_low = model.parameters["kf_irrev_heterotrimer_reassociation"].value / 100
kf_high = model.parameters["kf_irrev_heterotrimer_reassociation"].value * 100
swarm_param(
    model.parameters["kf_irrev_heterotrimer_reassociation"],
    loc=np.log10(kf_low),
    width=(np.log10(kf_high) - np.log10(kf_low)),
)
print(
    "kf_irrev_heterotrimer_reassociation: ",
    np.log10(kf_low),
    kf_low,
    np.log10(kf_high),
    kf_high,
)

# Low boundary for diffusion limited protein associations
# kf_low = 1e5  # (M*s)^-1
# # Upper limit for attractive associations
# kf_high = 1e10  # (M*s)^-1
# kf_low /= units.M_to_molec_per_pL * v_cell  # convert to 1/(s*molec/cell)
# kf_high /= units.M_to_molec_per_pL * v_cell  # convert to 1/(s*molec/cell)
# for param in ppf_params+bkf_params:
#     swarm_param(
#         model.parameters[param],
#         loc=np.log10(kf_low),
#         width=(np.log10(kf_high) - np.log10(kf_low)),
#     )

# for param in bkf_params:
#    swarm_param(model.parameters[param], loc=np.log10(kf_low), width=(np.log10(kf_high) - np.log10(kf_low)))

# Fix all the selected forward rate parameters to a
# the suggested plausible value 1e-6 (s*molec/pL)^-1 value from Aldridge et al.
# This is 6.02x10^5 1/(M*s) which falls in the range
# of diffusion limited association rates.
# kf = 1e-6 / v_cell # 1/(s*molec/pL) / v_cell
# broad range 10^5 - 10^10 1/M*s for diffusion limited
# to attactive associations
kf = 1e6  # (M*s)^-1
kf /= units.M_to_molec_per_pL * v_cell  # convert to 1/(s*molec/cell)
for param_name in ppf_params[3:] + bkf_params:
    print("Fixing forward rate: ", param_name)
    mask = [param.name == param_name for param in model.parameters]
    param_values[mask] = kf

# protein-protein dissociation rate constants
ppr_params = ["kr_PLC_bind_Gaq"]
# Other dissociations
bkr_params = ["kr_PLC_bind_PIP2"]
# Note plausible starting range from
# Aldridge et al. https://doi.org/10.1038/ncb1497
# is kr ~ 10^-3 - 10^-2 1/s.
kr_low = 1e-4  # s^-1 - 1 nM Kd for kf = 1x10^6 1/M*s
# Upper limit
kr_high = 1e4  # s^-1 - 1 mM Kd for kf = 1x10^6 1/M*s
for param in ppr_params:
    swarm_param(
        model.parameters[param],
        loc=np.log10(kr_low),
        width=(np.log10(kr_high) - np.log10(kr_low)),
    )

for param in bkr_params:
    swarm_param(
        model.parameters[param],
        loc=np.log10(kr_low),
        width=(np.log10(kr_high) - np.log10(kr_low)),
    )

# These have differnt diffusion limited forward rates:
# IP3
kf = model.parameters["kf_IP3_bind_IP3R"].value
kf_molar = (kf * v_cell) * units.M_to_molec_per_pL
# Note: Allow sub nM Kd since we don't model any initial IP3 in the cytosol.
# [IP3] ranged from about 40 nM up to 1.8 uM in Xenopus oocytes:
#    Luzzi et al. JBC 1998 https://doi.org/10.1074/jbc.273.44.28657
# High affinity pM IP3 binding has also been reported in Xenopus oocytes:
#    Demuro et al. Cell Calcium 2015 https://doi.org/10.1016/j.ceca.2015.08.003
kr_low = kf_molar * 1e-12  # kr for 1 pM Kd
kr_high = kf_molar * 1e-6  # kr for 1 uM Kd
swarm_param(
    model.parameters["kr_IP3_bind_IP3R"],
    loc=np.log10(kr_low),
    width=(np.log10(kr_high) - np.log10(kr_low)),
)
# Ca2+
kf = model.parameters["kf_erCa_bind_IP3R"].value
kf_molar = (kf * v_cell) * units.M_to_molec_per_pL
kr_low = kf_molar * 1e-9  # kr for 1 nM Kd
kr_high = kf_molar * 1e-3  # kr for 1 mM Kd
swarm_param(
    model.parameters["kr_erCa_bind_IP3R"],
    loc=np.log10(kr_low),
    width=(np.log10(kr_high) - np.log10(kr_low)),
)

# PAR2 conformational change to active state
# during two-state activation process:
#     LR <--> LR*
# Conformational change rate range
# Assume the forward rate is high (with ns to ms speeds)
kf_high = 1e2  # 1/s (1 per 1 ms, kf/kr=10000 when kr=1e-3)
kf_low = 1e-2  # 1/s (1 per 100 s, kf/kr=1 when kr=1e-3)
# Assume the reverse rate is slower (ms to s range speeds)
# kr_high = 1e2 # 1/s (1 per 10 ms)
# kr_low = 1e-2 # 1/s (1 per 100 s)
# swarm_param(model.parameters['k_activate_PAR2'], loc=np.log10(kf_low), width=(np.log10(kf_high) - np.log10(kf_low)))
kr_high = 1e5  # 1 / 10 us (kf/kr = 0.01 when kf=1e3)
kr_low = 0.1  # 10 / s (kf/kr = 10000 when kf=1e3)
# swarm_param(model.parameters['k_inactivate_PAR2'], loc=np.log10(kr_low), width=(np.log10(kr_high) - np.log10(kr_low)))

# Calcium transport through IP3R
k_low = (
    525 / 10
)  # 1/s Effective IP3R channel permeability as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
k_high = 525 * 10  # 1/s Lemon value x 100
# swarm_param(
#     model.parameters["k_tranport_erCa"],
#     loc=np.log10(k_low),
#     width=(np.log10(k_high) - np.log10(k_low)),
# )  # 1e2-1e8 1/s

# Dissociation constant for cytosolic Ca2+ binding and inhibiting IP3R
# We'll allow a broad range from 100 nM to 10 mM.
kd_low = np.log10(100 * units.nM_to_molec_per_pL * v_cell)  # 1 nM
kd_high = np.log10(10000 * units.microM_to_molec_per_pL * v_cell)  # 1 mM
# swarm_param(
#     model.parameters["Kd_cytCa_bind_IP3R"], loc=kd_low, width=(kd_high - kd_low)
# )

# IP3 degradation
kdeg_low = np.log10(1.25e-2)
kdeg_high = np.log10(1.25 * 10)
swarm_param(model.parameters["kdeg_ip3"], loc=kdeg_low, width=(kdeg_high - kdeg_low))

# PAR2 degradation
kdeg_low = np.log10(1e-5)
kdeg_high = np.log10(4e-3)
swarm_param(
    model.parameters["k_PAR2_bound_degradation"],
    loc=kdeg_low,
    width=(kdeg_high - kdeg_low),
)

# Calcium homeostasis reactions
# Cytosol to extracellular space:
# With previous model fits using just a single 1st order reaction to clear calcium we got around 4 per second.
# Since the net rate constant for Ca2+ clearance should be similar we'll
# assume this rate constant should be somewhere within an order of magnitude to that previous value.
k_low = 4e-2  # 1/s
k_high = 4e2  # 1/s
swarm_param(
    model.parameters["k_Ca_cyt_to_extra"],
    loc=np.log10(k_low),
    width=(np.log10(k_high) - np.log10(k_low)),
)
# Note that we'll couple the reverse rate constant k_Ca_extra_to_cyt to this
# one by enforcing that the rates should cancel for the initial calcium concentrations.
# Cytosol to ER lumen:
# Use the same range as
# swarm_param(model.parameters['k_Ca_cyt_to_er'], loc=np.log10(k_low), width=(np.log10(k_high)-np.log10(k_low)))
# Note that we'll couple the reverse rate constant k_Ca_er_to_cyt to this
# one by enforcing that the rates should cancel for the initial calcium concentrations.


kdeg_low = np.log10(
    1.25e-2
)  # 1/s -- 1.25/s from Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
kdeg_high = np.log10(1.25 * 10)
swarm_param(model.parameters["kdeg_ip3"], loc=kdeg_low, width=(kdeg_high - kdeg_low))


# Catalytic conversion of GTP to GDP when Gaq bound  to PLC
k_low = 1.33e-2  # 1/s - the value for Gaq alone
k_high = 1.33  # 1/s - assume it could be up to 100x faster.
swarm_param(
    model.parameters["k_gtp_to_gdp_plc"],
    loc=np.log10(k_low),
    width=(np.log10(k_high) - np.log10(k_low)),
)


# Catalytic conversion of PIP2 to IP3 by PLC when Ca is bound to PLC:
# Average value of PLCb3 and PLCb4 kcat's from Flaherty et al. 2008:
k_cat_flaherty = (22.85 + 27.89) / 2  # 1/s
# k_cat_flaherty = 27.89 # 1/s for PLCb4
k_high = 10 * k_cat_flaherty  # 1/s - up to 10x higher
k_low = k_cat_flaherty  # 1/s - up to 10x slower.
# swarm_param(
#     model.parameters["kcat_PIP2_to_IP3_Ca"],
#     loc=np.log10(k_low),
#     width=(np.log10(k_high) - np.log10(k_low)),
# )
# Catalytic conversion of PIP2 to IP3 by PLC when no Ca bound to PLC:
k_low = k_low * 0.01  # 1/s - up to 100x slower.
swarm_param(
    model.parameters["kcat_PIP2_to_IP3"],
    loc=np.log10(k_low),
    width=(np.log10(k_cat_flaherty) - np.log10(k_low)),
)
# Note: we are not explicitly enforcing kcat_PIP2_to_IP3_Ca >= kcat_PIP2_to_IP3 despite the
# implicit assumption that it probably should be in this case.

# Other parameter notes:
#  We'll leave the following at their nominal values based on literature or as
#  adapted from previous models:
#    * k_gtp_to_gdp_auto
#    * k_gtp_to_gdp_rgs
#    * kf_PLC_bind_Ca
#    * kr_PLC_bind_Ca
#    * kf_Ca_bind_buffer
#    * kr_Ca_bind_buffer
# We'll leave the following at their nominal diffusion limited
# values derived from literature reported diffusion coefficients:
#   * kf_cytCa_bind_IP3R
#   * kf_IP3_bind_IP3R
#   * kf_erCa_bind_IP3R

# Get the mask for the parameters being calibrated.
calibrate_mask = fancy_index(swarm_param, model.parameters)
print(
    "Model has {} parameters with {} being calibrated".format(
        len(model.parameters), len(calibrate_mask)
    )
)
print(calibrate_mask)
print("Parameters being calibrated:")
cal_param_names = swarm_param.names()
for name in cal_param_names:
    print("    " + name)
print(" ")
print("Parameters NOT being calibrated:")
for param in model.parameters:
    if param.name not in cal_param_names:
        print("    " + param.name)

# Mask for initial concentration of 2AT
twoat_mask = [param.name == "TAT_0" for param in model.parameters]
# Mask for the resting cytosolic Ca2+ amount for the FRET estimation.
cacrest_mask = [param.name == "Ca_C_resting" for param in model.parameters]
# Mask for the initial amount of Gprotein.
gprot_mask = [param.name == "Gprotein_0" for param in model.parameters]
# Mask for the initial amount of PAR2.
par2_mask = [param.name == "PAR2_0" for param in model.parameters]

# Masks for the calcium homeostasis reaction rate constants -- use for coupling
# the forward and reverse rates in the loglikelihood function.
kcaex2cyt_mask = [param.name == "k_Ca_extra_to_cyt" for param in model.parameters]
kcacyt2ex_mask = [param.name == "k_Ca_cyt_to_extra" for param in model.parameters]
kcacyt2er_mask = [param.name == "k_Ca_cyt_to_er" for param in model.parameters]
kcaer2cyt_mask = [param.name == "k_Ca_er_to_cyt" for param in model.parameters]

# Masks for the initial calcium amounts in different compartments.
cac0_mask = [param.name == "Ca_C_0" for param in model.parameters]
caex0_mask = [param.name == "Ca_extra_0" for param in model.parameters]
caer0_mask = [param.name == "Ca_E_0" for param in model.parameters]

# Mask for the PAR2 conformational change reverse (inactivation) rate.
kinactpar2_mask = [param.name == "k_inactivate_PAR2" for param in model.parameters]
kactpar2_mask = [param.name == "k_activate_PAR2" for param in model.parameters]
# Fix the rate:
param_values[kactpar2_mask] = 1e3  # 1 per ms
# param_values[kinactpar2_mask] = 1e-2 # 1 per 100 s

# Set the binding rate constant for cytosolic Ca2+ binding
# to IP3R to zero to effectively turn off the inhibition feedback.
# kfcacip3rbind_mask = [param.name=='kf_cytCa_bind_IP3R' for param in model.parameters]
# param_values[kfcacip3rbind_mask] = 0.

# set the amount of G-protein to equal that of PAR2
# param_values[gprot_mask] = 1000#param_values[par2_mask]
# Fix the amount of PAR2:
# param_values[par2_mask] = 160

# Get the peak value and location at 3160 nM 2AT for normalization.
norm_peak_idx = np.argmax(exp_data["3160_sig"].values[1:])
norm_peak_val = exp_data["3160_sig"].values[1:][norm_peak_idx]

# For fitting the FRET ratios at each 2AT concentration.
like_fret = dict()
logl_max = 0.0
peak_idxs = dict()
logl_maxs = dict()
logl_peaks = dict()
for data_set in data_sets:
    sig_name = data_set + "_sig"
    err_name = data_set + "_err"
    # Get experimental measurement and std. dev.
    y_exp = exp_data[sig_name].values[1:]
    sigma_exp = exp_data[err_name].values[1:]
    # Normalize:
    # y_exp /= norm_peak_val
    # sigma_exp /= norm_peak_val
    print(data_set)
    print(y_exp)
    print(sigma_exp)
    like_fret[data_set] = norm(loc=y_exp, scale=sigma_exp)
    mask = y_exp < 0.0
    y_exp[mask] = 0.0
    logl_max_d = np.sum(like_fret[data_set].logpdf(y_exp))
    logl_max += logl_max_d
    logl_maxs[data_set] = logl_max_d
    peak_idxs[data_set] = np.argmax(y_exp)
    logl_peaks[data_set] = like_fret[data_set].logpdf(y_exp)[peak_idxs[data_set]]

print("maximum loglikelihood: ", logl_max)

# We'll use the weights in the loglikelihood function to scale
# for importance of fitting each data point is propotional to
# how close it is to the peak value, which is used for the
# dose-response curve.
like_fret_weights = dict()
for data_set in data_sets:
    sig_name = data_set + "_sig"
    err_name = data_set + "_err"
    # Get experimental measurement and std. dev.
    y_exp = exp_data[sig_name].values[1:]
    sigma_exp = exp_data[err_name].values[1:]
    # using absolute values to prevent negative weights.
    weights = np.abs(y_exp) / np.max(np.abs(y_exp))
    weights_norm = weights / np.sum(weights)
    like_fret_weights[data_set] = weights_norm

# For comparing the dose-response data
# which uses the maximum FRET ratio at each 2AT
# concentration.
like_doseresonse = norm(
    loc=exp_data_2["Max.dRR"].values, scale=exp_data_2["Max.dRR_err"].values
)
logl_dr_max = np.sum(like_doseresonse.logpdf(exp_data_2["Max.dRR"].values))


# 600 s or 10 minutes
time_pre_equil = np.linspace(0, 1, 20, endpoint=True)


def loglikelihood(position):
    """Cost function.
    Based on determing the loglikelihood of simulated data and using
    the relationship:
        cost = - log(Likelihood)
    So, minimizing cost corresponds to maximizing the Likelihood.
    This cost includes components for fitting the training data at each 2AT
    concentration, as well as added weight for fitting the peak in the experimental
    curve, additional contribution from fitting the slope of the curve, and
    a contribution accouting for upward shifting in time of the peak at 330 nM 2AT
    after 50% of the PAR2 is inactivated.
    """
    Y = 10 ** position
    # Start with first soft constraint (data fit) data set so that if
    # the parameter set is really bad (ODE integrator problems)
    # we can exit having only tried to run one simulation.
    param_values[calibrate_mask] = Y

    # Get the initial amounts of calcium.
    cac0 = param_values[cac0_mask]
    caex0 = param_values[caex0_mask]
    caer0 = param_values[caer0_mask]
    # Now couple the reverse rate constant to the forward (being sampled) for
    # calcium homeostasis reactions.
    # kcacyt2ex = param_values[kcacyt2ex_mask]
    # kcacyt2er = param_values[kcacyt2er_mask]
    # kcaex2cyt = kcacyt2ex * cac0 / caex0
    # kcaer2cyt = kcacyt2er * cac0 / caer0

    # param_values[kcaex2cyt_mask] = kcaex2cyt
    # param_values[kcaer2cyt_mask] = kcaer2cyt
    conc_eq = None
    # Try the preequilibration.
    # try:
    #     param_values[twoat_mask] = 0
    #     #t_idx, conc_eq, cytoCa_eq, erCa_eq = util.pre_equilibration(model, time_pre_equil, parameters=param_values, tolerance=100)
    #     #conc_eq = conc_eq[0]
    #     #cytoCa_eq = cytoCa_eq[0]
    #     #erCa_eq = erCa_eq[0]
    #     param_values_eq, conc_eq = util.pre_equilibrate(model, time_pre_equil,
    #                                                     param_values=param_values,
    #                                                     tolerance=1000)
    #     param_values[:] = param_values_eq[:]
    # except:
    #     return -np.inf
    logp_vals = list()
    peak_vals = list()
    # Now simulate each 2AT concentration and get loglikelihood contribution for
    # fitting the FRET time series.
    sim_norm_val = 1.0
    for i, data_set in enumerate(data_sets):
        sig_name = data_set + "_sig"
        err_name = data_set + "_err"
        tat_conc = float(data_set) * units.nM_to_molec_per_pL * v_extra
        param_values[twoat_mask] = tat_conc
        # param_values[cacrest_mask] = cytoCa_eq
        # conc_eq[0] = tat_conc
        sim = solver.run(param_values=param_values, initials=conc_eq).all
        # Get the simulated FRET response
        y_sim = sim["FRET"][fidx][1:]
        y_peak = y_sim[peak_idxs[data_set]]
        peak_vals.append(y_peak)
        # Normalize
        if i == 0:
            sim_norm_val = y_sim[norm_peak_idx]
        # y_sim /= sim_norm_val
        logp_s = like_fret[data_set].logpdf(y_sim)
        logp_peak = logp_s[peak_idxs[data_set]] / logl_peaks[data_set]
        # The additional weights scale the importance of fitting each
        # data point relative to how close it is to peak value, which
        # is used in determining the dose-response.
        # logp_w = logp_s * like_fret_weights[data_set]
        # logp = np.sum(logp_w)
        # Compute error between simulation and experiment
        # if np.isnan(logp):
        #    return -np.inf
        # logp_vals.append(logp)
        logp = np.sum(logp_s) / logl_maxs[data_set]
        if np.isnan(logp):
            return -np.inf
        logp_vals.append(logp)
        if len(data_sets) < 6:
            # logp_vals.append(logp_peak / 6)
            logp_vals[-1] += logp_peak / 6
    # Add the dose-response component, but only if all
    # six concentrations are included.
    if len(data_sets) == 6:
        logp_dr = np.sum(like_doseresonse.logpdf(peak_vals)) / logl_dr_max
        logp_vals.append(logp_dr)

    # log(Likelihood)
    return np.sum(logp_vals)


def chi2(position):
    """Cost function.
    Based on determing the loglikelihood of simulated data and using
    the relationship:
        cost = - log(Likelihood)
    So, minimizing cost corresponds to maximizing the Likelihood.
    This cost includes components for fitting the training data at each 2AT
    concentration, as well as added weight for fitting the peak in the experimental
    curve, additional contribution from fitting the slope of the curve, and
    a contribution accouting for upward shifting in time of the peak at 330 nM 2AT
    after 50% of the PAR2 is inactivated.
    """
    Y = 10 ** position
    # Start with first soft constraint (data fit) data set so that if
    # the parameter set is really bad (ODE integrator problems)
    # we can exit having only tried to run one simulation.
    param_values[calibrate_mask] = Y

    # Try the preequilibration.
    try:
        param_values[twoat_mask] = 0
        # t_idx, conc_eq, cytoCa_eq, erCa_eq = util.pre_equilibration(model, time_pre_equil, parameters=param_values, tolerance=100)
        # conc_eq = conc_eq[0]
        # cytoCa_eq = cytoCa_eq[0]
        # erCa_eq = erCa_eq[0]
        param_values_eq, conc_eq = util.pre_equilibrate(
            model, time_pre_equil, param_values=param_values, tolerance=1000
        )
        param_values[:] = param_values_eq[:]
    except:
        return -np.inf
    chi2_vals = list()
    # Now simulate each 2AT concentration and get loglikelihood contribution for
    # fitting the FRET time series.
    sim_norm_val = 1.0
    for i, data_set in enumerate(data_sets):
        sig_name = data_set + "_sig"
        err_name = data_set + "_err"
        tat_conc = float(data_set) * units.nM_to_molec_per_pL * v_extra
        param_values[twoat_mask] = tat_conc
        # param_values[cacrest_mask] = cytoCa_eq
        conc_eq[0] = tat_conc
        sim = solver.run(param_values=param_values, initials=conc_eq).all
        # Get the simulated FRET response
        y_sim = sim["FRET"][fidx][1:]
        # Normalize
        if i == 0:
            sim_norm_val = y_sim[norm_peak_idx]
        # y_sim /= sim_norm_val
        sig_name = data_set + "_sig"
        err_name = data_set + "_err"
        # Get experimental measurement and std. dev.
        y_exp = exp_data[sig_name].values[1:]
        sigma_exp = exp_data[err_name].values[1:]
        chi2_s = (y_sim - y_exp) ** 2 / sigma_exp ** 2
        chi2_peak = chi2_s[peak_idxs[data_set]]
        chi2_avg = np.mean(chi2_s)
        # The additional weights scale the importance of fitting each
        # data point relative to how close it is to peak value, which
        # is used in determining the dose-response.
        # logp_w = logp_s * like_fret_weights[data_set]
        # logp = np.sum(logp_w)
        # Compute error between simulation and experiment
        # if np.isnan(logp):
        #    return -np.inf
        # logp_vals.append(logp)
        # logp = np.sum(logp_s)
        chi2p = max([chi2_avg, chi2_peak])
        if np.isnan(chi2p):
            return -np.inf
        chi2_vals.append(chi2p)
    # log(Likelihood)
    return np.sum(chi2_vals)


def cost(theta):
    # return chi2(theta)
    return -loglikelihood(theta)


# Function to display comparisons of Exp results and simulated results.
def display(positions, include="all", id=0):
    # Y = 10**position
    # Y[0] *=0.5
    new_positions = list()
    eq_concs = list()
    for position in positions:
        Y = 10 ** position
        # Start with first soft constraint (data fit) data set so that if
        # the parameter set is really bad (ODE integrator problems)
        # we can exit having only tried to run one simulation.
        param_values[calibrate_mask] = Y
        # param_values[twoat_mask] = 0
        # # t_idx, conc_eq, cytoCa_eq, erCa_eq = util.pre_equilibration(model, time_pre_equil, parameters=param_values, tolerance=100)
        # # conc_eq = conc_eq[0]
        # # cytoCa_eq = cytoCa_eq[0]
        # # erCa_eq = erCa_eq[0]
        # param_values_eq, conc_eq = util.pre_equilibrate(
        #     model, time_pre_equil, param_values=param_values, tolerance=1000
        # )
        # param_values[:] = param_values_eq[:]
        new_positions.append(param_values.copy())
        # eq_concs.append(conc_eq.copy())
        eq_concs = None
    times = exp_data["Time"]
    if (include == "all") or ("10" in include):
        # 10 nM 2AT
        tat_conc = 10 * units.nM_to_molec_per_pL * v_extra
        for i, position in enumerate(new_positions):
            position[twoat_mask] = tat_conc
            # eq_concs[i][0] = tat_conc
        sim = solver.run(
            param_values=new_positions, initials=eq_concs, num_processors=4
        ).all
        # print(sim[0]['FRET'])
        fret_sims = list()
        for s in sim:
            fret_sims.append(s["FRET"])
        # fret_sims = np.array(fret_sims)
        # fret_mean = np.mean(fret_sims, axis=0)
        # fret_std = np.std(fret_sims, axis=0)
        # print(fret_mean)
        # print(sim['FRET'][0])
        # return
        # Get the simulated FRET response
        # ysim_norm = sim['FRET']
        plt.errorbar(
            times,
            exp_data["10_sig"],
            exp_data["10_err"],
            capsize=2.5,
            label="10",
            color="k",
            marker="s",
            linestyle="",
        )
        for traj in fret_sims:
            plt.plot(tspan, traj, label=None, color="k")
        # plt.fill_between(
        #     tspan,
        #     fret_mean - 1.96 * fret_std,
        #     fret_mean + 1.96 * fret_std,
        #     alpha=0.25,
        #     color="k",
        # )
        # plt.show()

        # return
    if (include == "all") or ("31.6" in include):
        # 31.6 nM 2AT
        tat_conc = 31.6 * units.nM_to_molec_per_pL * v_extra
        for i, position in enumerate(new_positions):
            position[twoat_mask] = tat_conc
            # eq_concs[i][0] = tat_conc
        sim = solver.run(
            param_values=new_positions, initials=eq_concs, num_processors=4
        ).all
        # Get the simulated FRET response
        # ysim_norm = sim['FRET']
        fret_sims = list()
        for s in sim:
            fret_sims.append(s["FRET"])
        # fret_sims = np.array(fret_sims)
        # fret_mean = np.mean(fret_sims, axis=0)
        # fret_std = np.std(fret_sims, axis=0)
        # print(fret_std)
        plt.errorbar(
            times,
            exp_data["31.6_sig"],
            exp_data["31.6_err"],
            capsize=2.5,
            label="31.6",
            color="r",
            marker="o",
            linestyle="",
        )
        for fret_mean in fret_sims:
            plt.plot(tspan, fret_mean, label=None, color="r")
        # plt.fill_between(
        #     tspan,
        #     fret_mean - 1.96 * fret_std,
        #     fret_mean + 1.96 * fret_std,
        #     alpha=0.25,
        #     color="r",
        # )

    # 100 nM 2AT
    if (include == "all") or ("100" in include):
        tat_conc = 100 * units.nM_to_molec_per_pL * v_extra
        for i, position in enumerate(new_positions):
            position[twoat_mask] = tat_conc
            # eq_concs[i][0] = tat_conc
        sim = solver.run(
            param_values=new_positions, initials=eq_concs, num_processors=4
        ).all
        # Get the simulated FRET response
        # ysim_norm = sim['FRET']
        fret_sims = list()
        for s in sim:
            fret_sims.append(s["FRET"])
        fret_sims = np.array(fret_sims)
        fret_mean = np.mean(fret_sims, axis=0)
        fret_std = np.std(fret_sims, axis=0)
        plt.errorbar(
            times,
            exp_data["100_sig"],
            exp_data["100_err"],
            capsize=2.5,
            label="100",
            color="b",
            marker="^",
            linestyle="",
        )
        plt.plot(tspan, fret_mean, label=None, color="b")
        plt.fill_between(
            tspan,
            fret_mean - 1.96 * fret_std,
            fret_mean + 1.96 * fret_std,
            alpha=0.25,
            color="b",
        )
        # plt.show()
        for fret_mean in fret_sims:
            plt.plot(tspan, fret_mean, label=None, color="b", linestyle="--")

        # return
    # 316 nM 2AT
    if (include == "all") or ("316" in include):
        tat_conc = 316 * units.nM_to_molec_per_pL * v_extra
        for i, position in enumerate(new_positions):
            position[twoat_mask] = tat_conc
            # eq_concs[i][0] = tat_conc
        sim = solver.run(
            param_values=new_positions, initials=eq_concs, num_processors=4
        ).all
        # Get the simulated FRET response
        # ysim_norm = sim['FRET']
        fret_sims = list()
        for s in sim:
            fret_sims.append(s["FRET"])
        fret_sims = np.array(fret_sims)
        fret_mean = np.mean(fret_sims, axis=0)
        fret_std = np.std(fret_sims, axis=0)
        plt.errorbar(
            times,
            exp_data["316_sig"],
            exp_data["316_err"],
            capsize=2.5,
            label="316",
            color="m",
            marker="v",
            linestyle="",
        )
        plt.plot(tspan, fret_mean, label=None, color="m", linewidth=2)

        plt.fill_between(
            tspan,
            fret_mean - 1.96 * fret_std,
            fret_mean + 1.96 * fret_std,
            alpha=0.25,
            color="m",
        )
        for fret_mean in fret_sims:
            plt.plot(
                tspan, fret_mean, label=None, color="m", linewidth=2, linestyle="--"
            )
    # 1000 nM 2AT
    if (include == "all") or ("1000" in include):
        tat_conc = 1000 * units.nM_to_molec_per_pL * v_extra
        for i, position in enumerate(new_positions):
            new_positions[i][twoat_mask] = tat_conc
            # eq_concs[i][0] = tat_conc
        sim = solver.run(
            param_values=new_positions, initials=eq_concs, num_processors=4
        ).all
        # Get the simulated FRET response
        # ysim_norm = sim['FRET']
        fret_sims = list()
        for s in sim:
            fret_sims.append(s["FRET"])
        fret_sims = np.array(fret_sims)
        fret_mean = np.mean(fret_sims, axis=0)
        fret_std = np.std(fret_sims, axis=0)
        plt.errorbar(
            times,
            exp_data["1000_sig"],
            exp_data["1000_err"],
            capsize=2.5,
            label="1000",
            color="g",
            marker="D",
            linestyle="",
        )
        plt.plot(tspan, fret_mean, label=None, color="g", linewidth=2)
        plt.fill_between(
            tspan,
            fret_mean - 1.96 * fret_std,
            fret_mean + 1.96 * fret_std,
            alpha=0.25,
            color="g",
        )
        for fret_mean in fret_sims:
            plt.plot(
                tspan, fret_mean, label=None, color="g", linewidth=2, linestyle="--"
            )
    # 3160 nM 2AT
    if (include == "all") or ("3160" in include):
        tat_conc = 3160 * units.nM_to_molec_per_pL * v_extra
        for i, position in enumerate(new_positions):
            new_positions[i][twoat_mask] = tat_conc
            # eq_concs[i][0] = tat_conc
        sim = solver.run(
            param_values=new_positions, initials=eq_concs, num_processors=4
        ).all
        # Get the simulated FRET response
        # ysim_norm = sim['FRET']
        fret_sims = list()
        for s in sim:
            fret_sims.append(s["FRET"])
        # fret_sims = np.array(fret_sims)
        # fret_mean = np.mean(fret_sims, axis=0)
        # fret_std = np.std(fret_sims, axis=0)
        plt.errorbar(
            times,
            exp_data["3160_sig"],
            exp_data["3160_err"],
            capsize=2.5,
            label="3160",
            color="darkblue",
            marker="<",
            linestyle="",
        )
        for fret_mean in fret_sims:
            plt.plot(tspan, fret_mean, label=None, color="darkblue", linewidth=2)
        # plt.fill_between(
        #     tspan,
        #     fret_mean - 1.96 * fret_std,
        #     fret_mean + 1.96 * fret_std,
        #     alpha=0.25,
        #     color="darkblue",
        # )
    plt.xlabel("Time (s)")
    plt.ylabel("FRET Signal")
    plt.legend(loc=0)
    savename = "fit-comparison_exp-plot_id{}.png".format(id)
    plt.savefig(savename)
    # plt.show()
    plt.close()


# Set the number of particles in the swarm.
num_particles = 60
# Set the number of iterations for PSO run.
num_iterations = 40
# Number of processors for PSO multiprocessing.
nproc = 6
costs = list()
parms = list()
# Number of different PSO runs to do.
ntries = 10
print(
    "Doing PSO {} times with {} particles and {} iterations".format(
        ntries, num_particles, num_iterations
    )
)
print("Using {} processors".format(nproc))
for i in range(ntries):
    print("Running PSO try {}".format(i))
    # Construct the optimizer
    pso = PSO(save_sampled=False, verbose=True)
    starting_position = swarm_param.centers()
    pso.set_start_position(starting_position)
    pso.set_bounds(lower=swarm_param.lower(), upper=swarm_param.upper())
    # sets maximum speed that a particle can travel
    pso.set_speed(-0.5, 0.5)
    # run it
    # data_sets = ["10", "3160"]
    best = None
    data_sets = ["31.6", "1000"]
    # ds = data_sets_all[1:]
    # random.shuffle(ds)
    # while len(ds) > 0:
    # data_sets.append(ds.pop())
    # data_sets = sorted(data_sets, key= lambda x: float(x))
    # print(data_sets)
    pso.run(
        num_particles,
        num_iterations,
        cost_function=cost,
        stop_threshold=1e-3,
        num_processors=nproc,
    )
    best = pso.best.pos
    pso.set_start_position(best)
    if pso.best.fitness < 0.0:
        costs.append(pso.best.fitness)
        parms.append(best)
    print("best at try {} with cost {} :".format(i, pso.best.fitness))
    print(list(best))
    for i, idx in enumerate(calibrate_mask):
        line = "{} {} {} {}".format(
            idx, model.parameters[int(idx)].name, best[i], 10.0 ** best[i]
        )
        print(line)

costs = np.array(costs)
parms = np.array(parms)
sort_idx = np.argsort(costs)
for i, parm in enumerate(parms[sort_idx]):
    np.save("pso_best_chain_{}".format(i), parm)

parms_1 = parms.copy()
np.savetxt("pso_best_costs.csv", costs[sort_idx], delimiter=",")

display(parms, include=data_sets, id=2)

for i, idx in enumerate(calibrate_mask):
    line = "{} {}".format(
        idx,
        model.parameters[int(idx)].name,
    )
    print(line)
    values_1 = parms[:, i]
    # values_2 = parms_2[:, i]
    plt.title(line)
    plt.hist(values_1, label="31.6")
    # plt.hist(values_2, label="1000")
    plt.legend(loc=0)
    plt.show()
