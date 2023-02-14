import os
import numpy as np
import pysb
from pysb.simulator import ScipyOdeSimulator
from itertools import compress
import typing
from parm import parameter_masks, default_param_values
from parm import units


def expand_times(times, expand_by=100):
    # tsp = np.linspace(0, times.max(), expand_by)
    # tsp = np.concatenate((tsp, times[1:-1]))
    # tsp.sort()
    tmax = np.max(times) + 5
    tsp = np.arange(tmax)
    findex = list()
    for time in times:
        idx = np.where(tsp == time)[0][0]
        findex.append(idx)
    return tsp, findex


def run_model(
    model: pysb.Model,
    tspan: np.ndarray,
    param_values: typing.Union[None, np.ndarray, typing.List[np.ndarray]] = None,
    initials: typing.Union[None, np.ndarray, typing.List[np.ndarray]] = None,
    nprocs: int = 1,
) -> np.ndarray:
    """Run the given model using ScipyOdeSimulator.

    This is a wrapper function that simulates the model using the ScipyOdeSimulator
    with the lsoda integrator and returns the corresponding model trajectory.

    Args:
        model: The input PySB model to simulate.
        tspan: The time span to simulate the model over.
        param_values: Optional specification of parameters to use when
            simulating the model. If None, the nominal/default model parameters
            will be used. The input can be None, a single parameter vector, or
            a list of parameter vectors that will each be simulated.
        initials: Optional specification of initial concentrations to use when
            simulating the model. If None, the nominal/default model values
            will be used. The input can be None, a single vector, or
            a list of vectors that will each be simulated.

    Returns:
        The PySB model simulation trajectory as a structured NumPy array.
    """

    solver = ScipyOdeSimulator(
        model,
        tspan=tspan,
        integrator="lsoda",
    )
    m_run = solver.run(
        param_values=param_values, initials=initials, num_processors=nprocs
    )
    yout = m_run.all
    return yout


# Adapted from: https://github.com/LoLab-VU/JARM/blob/master/model_analysis/equilibration_function.py
def _pre_equilibration(model, time_search, parameters=None, tolerance=1e-6):
    """
    Parameters
    ----------
    model : pysb.Model
        pysb model
    time_search : np.array
        Time span array used to find the equilibrium
    parameters :  dict or np.array
        Model parameters used to find the equilibrium, it can be an array with all model parameters
        (this array must have the same order as model.parameters) or it can be a dictionary where the
        keys are the parameter names thatpara want to be changed and the values are the new parameter
        values.
    tolerance : float
        Tolerance to define when the equilibrium has been reached
    Returns
    -------
    """
    # Solve system for the time span provided
    solver = ScipyOdeSimulator(
        model, tspan=time_search, param_values=parameters, integrator="lsoda"
    ).run()
    if solver.nsims == 1:
        simulations = [solver.species]
    else:
        simulations = solver.species

    dt = time_search[1] - time_search[0]

    all_times_eq = [0] * solver.nsims
    all_conc_eq = [0] * solver.nsims
    for n, y in enumerate(simulations):
        time_to_equilibration = [0, 0]
        for idx in range(y.shape[1]):
            sp_eq = False
            derivative = np.diff(y[:, idx]) / dt
            derivative_range = (derivative < tolerance) & (derivative > -tolerance)
            # Indexes of values less than tolerance and greater than -tolerance
            derivative_range_idxs = list(
                compress(range(len(derivative_range)), derivative_range)
            )
            for i in derivative_range_idxs:
                # Check if derivative is close to zero in the time points ahead
                if i + 3 > len(time_search):
                    raise Exception(
                        "Equilibrium can not be reached within the time_search input"
                    )
                if (derivative[i + 3] < tolerance) | (derivative[i + 3] > -tolerance):
                    sp_eq = True
                    if time_search[i] > time_to_equilibration[0]:
                        time_to_equilibration[0] = time_search[i]
                        time_to_equilibration[1] = i
                if not sp_eq:
                    raise Exception(
                        "Equilibrium can not be reached within the time_search input"
                    )
                if sp_eq:
                    break
            else:
                raise Exception(
                    "Species s{0} ({1}) has not reached equilibrium".format(
                        idx, model.species[idx]
                    )
                )

        conc_eq = y[time_to_equilibration[1]]
        all_times_eq[n] = time_to_equilibration
        all_conc_eq[n] = conc_eq
    if solver.nsims == 1:
        cytoCa_eq = [solver.observables["cytoCa"][time_to_equilibration[1]]]
        erCa_eq = [solver.observables["erCa"][time_to_equilibration[1]]]
    else:
        cytoCa_eq = [
            obs["cytoCa"][time_to_equilibration[1]] for obs in solver.observables
        ]
        erCa_eq = [obs["erCa"][time_to_equilibration[1]] for obs in solver.observables]
    return all_times_eq, all_conc_eq, cytoCa_eq, erCa_eq


def pre_equilibrate(
    model: pysb.Model,
    time_search: np.ndarray,
    param_values: typing.Optional[np.ndarray] = None,
    tolerance: float = 10.0,
) -> typing.Tuple[np.ndarray, np.ndarray]:

    # Get a vector of the parameter values.
    if param_values is None:
        param_values = np.copy(default_param_values)
    # Make a mask for the initial concentration of the agonist, 2AT, and set
    # it to zero.
    two_at_mask = parameter_masks["TAT_0"]
    two_at_initial = param_values[two_at_mask]
    param_values[two_at_mask] = 0
    t_eq, conc_eq, cytoCa_eq, erCa_eq = _pre_equilibration(
        model, time_search, parameters=param_values, tolerance=tolerance
    )

    # Mask for the resting cytosolic Ca2+ amount for the FRET estimation.
    cacrest_mask = parameter_masks["Ca_C_resting"]
    # Adjust the FRET parameter affected by the initial cytosol calcium concentration.
    param_values[cacrest_mask] = cytoCa_eq[0]
    # Reset the intial 2AT concentration.
    param_values[two_at_mask] = two_at_initial
    initials = conc_eq[0]

    return param_values, initials


def calcium_homeostasis_reverse_rate_coupling(
    param_values: np.array,
) -> np.array:
    """Adjust parameters to couple reverse rates of calcium homeostasis to the forward rates."""

    # Masks for the calcium homeostasis reaction rate constants -- use for coupling
    # the forward and reverse rates in the loglikelihood function.
    kcaex2cyt_mask = parameter_masks["k_Ca_extra_to_cyt"]
    kcacyt2ex_mask = parameter_masks["k_Ca_cyt_to_extra"]
    kcacyt2er_mask = parameter_masks["k_Ca_cyt_to_er"]
    kcaer2cyt_mask = parameter_masks["k_Ca_er_to_cyt"]
    # Masks for the initial calcium amounts in different compartments.
    cac0_mask = parameter_masks["Ca_C_0"]
    caex0_mask = parameter_masks["Ca_extra_0"]
    caer0_mask = parameter_masks["Ca_E_0"]
    # Get the initial amounts of calcium.
    cac0 = param_values[cac0_mask]
    caex0 = param_values[caex0_mask]
    caer0 = param_values[caer0_mask]
    # Now couple the reverse rate constant to the forward for
    # calcium homeostasis reactions.
    kcacyt2ex = param_values[kcacyt2ex_mask]
    kcacyt2er = param_values[kcacyt2er_mask]
    param_values[kcaex2cyt_mask] = kcacyt2ex * cac0 / caex0
    param_values[kcaer2cyt_mask] = kcacyt2er * cac0 / caer0
    return param_values


def get_parameter_vector(
    model: pysb.Model,
) -> np.array:
    """Returns a parameter value vector corresponding to a PySB model's parameters."""
    return np.array([param.value for param in model.parameters])


def get_parameter_mask(
    model: pysb.Model,
    param_name: str,
) -> typing.List[bool]:
    """Returns a boolean mask of the parameter value vector for given PySB model and parameter name."""
    return [par.name == param_name for par in model.parameters]


def set_tat_initial_nM(
    param_values: np.array,
    tat_conc_nM: float,
) -> np.array:
    """Returns an updated parameter vector with the desired 2AT initial concentration."""

    # Mask for initial concentration of 2AT
    twoat_mask = parameter_masks["TAT_0"]
    vextra_mask = parameter_masks["Vextra"]
    v_extra = param_values[vextra_mask]
    param_values[twoat_mask] = tat_conc_nM * units.nM_to_molec_per_pL * v_extra
    return param_values

def reduce_par2(
    param_values: np.array,
    fraction: float,
    idx: int = 0,
) -> np.array:
    """Returns an updated parameter vector with the desired 2AT initial concentration."""


    par2_0 = param_values[idx]
    param_values[idx] = par2_0 * fraction 
    return param_values
  

def load_pydream_chains(
    niter: int,
    nchains: int,
    fpath="./pydream_results",
    model="parm",
    dream_type="dreamzs",
    burnin=None,
) -> np.array:
    """Loads parameter samples from PyDREAM chains and concatenates them into a single 2D array."""

    chains = list()
    if burnin is None:
        burnin = int(niter / 2)
    for i in range(nchains):
        fname = "{}_{}_{}chains_sampled_params_chain_{}_{}.npy".format(
            model, dream_type, nchains, i, niter
        )
        fname = os.path.join(fpath, fname)
        chain = np.load(os.path.abspath(fname))
        chains.append(chain[burnin:])

    return np.concatenate(chains)

def load_pydream_likelihoods(
    niter: int,
    nchains: int,
    fpath="./pydream_results",
    model="parm",
    dream_type="dreamzs",
    burnin=None,
) -> np.array:
    """Loads parameter samples from PyDREAM chains and concatenates them into a single 2D array."""

    chains = list()
    if burnin is None:
        burnin = int(niter / 2)
    for i in range(nchains):
        fname = "{}_{}_{}chains_logps_chain_{}_{}.npy".format(
            model, dream_type, nchains, i, niter
        )
        fname = os.path.join(fpath, fname)
        chain = np.load(os.path.abspath(fname))
        chains.append(chain[burnin:])

    return np.concatenate(chains)