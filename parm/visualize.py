import typing
# NumPy, Pandas, and Matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysb

from . import data
from . import util
from . import units

EXP_2AT_CONCS = {"10":10.,
                  "31.6":31.6,
                  "100":100.,
                  "316":316.,
                  "1000":1000.,
                  "3160":3160.}

EXP_PLT_COLORS = {"10":'black',
                  "31.6":'red',
                  "100":'blue',
                  "316":'magenta',
                  "1000":'green',
                  "3160":'darkblue'}

EXP_PLT_MARKERS = {"10":'s',
                  "31.6":'o',
                  "100":'^',
                  "316":'v',
                  "1000":'D',
                  "3160":'v'}

def kang2019_fig_2d():
    fig2d = data.kang2019_fig_2d()
    plt.errorbar(fig2d['Time'],fig2d['(a)_sig'], fig2d['(a)_err'], marker='s', markersize=5, color='black', label='(a)', capsize=3)
    plt.errorbar(fig2d['Time'],fig2d['(b)_sig'], fig2d['(b)_err'], marker='o', markersize=5, color='red', label='(b)', capsize=3)
    plt.errorbar(fig2d['Time'],fig2d['(c)_sig'], fig2d['(c)_err'], marker='^', markersize=5, color='blue', label='(c)', capsize=3)
    plt.errorbar(fig2d['Time'],fig2d['(d)_sig'], fig2d['(d)_err'], marker='v', markersize=5, color='magenta', label='(d)', capsize=3)
    plt.ylabel(r'$\Delta$R/R')
    plt.xlabel('Time (s)')
    plt.legend(loc='upper left')
    return

def kang2019_fig_s3c():
    figs3c = data.kang2019_fig_s3c()
    times = figs3c['Time']
    plt.errorbar(times, figs3c['10_sig'], figs3c['10_err'],
                 capsize=2.5, label='10', color='black', marker='s')
    plt.errorbar(times, figs3c['31.6_sig'], figs3c['31.6_err'],
                 capsize=2.5, label='31.6', color='red', marker='o')
    plt.errorbar(times, figs3c['100_sig'], figs3c['100_err'],
                 capsize=2.5, label='100', color='blue', marker='^')
    plt.errorbar(times, figs3c['316_sig'], figs3c['316_err'],
             capsize=2.5, label='316', color='magenta', marker='v')
    plt.errorbar(times, figs3c['1000_sig'], figs3c['1000_err'],
                 capsize=2.5, label='1000', color='green', marker='D')
    plt.errorbar(times, figs3c['3160_sig'], figs3c['3160_err'],
                 capsize=2.5, label='3160', color='darkblue', marker='<')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$\Delta$R/R')
    plt.legend(loc=0, title="[2AT] (nM)")

def kang2019_fig_s3d():
    figs3d = data.kang2019_fig_s3d()
    plt.errorbar(figs3d['[2AT]'], figs3d['Max.dRR'], figs3d['Max.dRR_err'], marker='s', linestyle="", color='black', capsize=2.5, label='Exp.')
    plt.xscale('log')
    plt.ylabel(r'Max. $\Delta$R/R')
    plt.xlabel('[2AT] (nM)')
    plt.xlim((1, 5000))

# Function to display comparisons of Exp results and simulated results.
def display_expcomp_single(
    model: pysb.Model,
    position: np.array,
    idxs_mask: list,
    nprocs: int = 1,
    save: bool = False,
    ):
    """Runs the model and plots a comparison with the experimental data for 2AT response."""

    Y = 10 ** position
    param_values = np.array([param.value for param in model.parameters])
    twoat_mask = [param.name == "TAT_0" for param in model.parameters]
    param_values[idxs_mask] = Y
    # Experimental data is in PARM/exp_data
    exp_data = data.training_data()
    times = exp_data["Time"]
    tspan, fidx = util.expand_times(times)
    kang2019_fig_s3c()

    positions = list()
    for tat_conc in EXP_2AT_CONCS.values():
        param_values[twoat_mask] = tat_conc * units.nM_to_molec_per_pL * model.parameters["Vextra"].value
        positions.append(param_values.copy())
    sim_out = util.run_model(model, tspan, param_values=positions, nprocs=nprocs)
    for i, key in enumerate(EXP_2AT_CONCS.keys()):
        ysim = sim_out[i]["FRET"]
        plt.plot(tspan, ysim, label=None, color=EXP_PLT_COLORS[key])
    if save:
        plt.savefig("display_expcomp_single.png")

def display_expcomp_multi_grid_mean_ci(
    model: pysb.Model,
    positions: typing.Union[np.ndarray, typing.List[np.ndarray]],
    idxs_mask: list,
    counts: typing.Union[None, np.ndarray] = None,
    nprocs: int = 1,
    save: bool = False,
    ):
    """Runs the model and plots a comparison with the experimental data for 2AT response."""


    param_values = np.array([param.value for param in model.parameters])
    twoat_mask = [param.name == "TAT_0" for param in model.parameters]

    # Experimental data is in PARM/exp_data
    exp_data = data.training_data()
    times = exp_data["Time"]
    tspan, fidx = util.expand_times(times)
    if counts is None:
        counts = np.array([1]*len(positions))
    f, axes = plt.subplots(2, 3, figsize=(7, 4), sharex=True, sharey=True)
    new_positions = list()
    for position in positions:
        Y = 10 ** position
        param_values[idxs_mask] = Y
        new_positions.append(param_values.copy())
    for i, key in enumerate(EXP_2AT_CONCS.keys()):
        tat_number = EXP_2AT_CONCS[key] * units.nM_to_molec_per_pL * model.parameters["Vextra"].value
        for position in new_positions:
            position[twoat_mask] = tat_number
        sim_out = util.run_model(model, tspan, param_values=new_positions, nprocs=nprocs)
        fret_sims = list()
        for s in sim_out:
            fret_sims.append(s["FRET"])
        fret_sims = np.array(fret_sims)
        fret_mean = np.average(fret_sims, axis=0, weights=counts)
        fret_std = np.sqrt(
            np.average(np.subtract(fret_sims, fret_mean) ** 2, axis=0, weights=counts)
        )
        col_sig = "{}_sig".format(key)
        col_err = "{}_err".format(key)
        # 0: [0,0]
        # 1: [0,1]
        # 2: [0, 2]
        # 3: [1, 0]
        # 4: [1, 1]
        # 5: [1, 2]
        ax_row = int(i / 3)
        ax_col = int((i - 3*ax_row))
        print(i, ax_row, ax_col)
        axes[ax_row, ax_col].errorbar(times, exp_data[col_sig], exp_data[col_err],
                     capsize=2.5, label="Exp.", color=EXP_PLT_COLORS[key],
                     marker=EXP_PLT_MARKERS[key], linestyle=' ')

        axes[ax_row, ax_col].plot(tspan, fret_mean, color=EXP_PLT_COLORS[key],
                        linewidth=2, label="Sim.")
        axes[ax_row, ax_col].fill_between(
            tspan,
            fret_mean - 1.96 * fret_std,
            fret_mean + 1.96 * fret_std,
            alpha=0.25,
            color=EXP_PLT_COLORS[key],
        )
        axes[ax_row, ax_col].set_title("[2AT] = {} nM".format(key), fontdict={"fontsize": 10})
        axes[ax_row, ax_col].legend(loc=0, frameon=False)           
    for ax in axes.flat:
        ax.set(xlabel="Time (s)", ylabel=r"$\Delta$R/R ")
    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axes.flat:
        ax.label_outer()
    plt.tight_layout()    
    if save:
        plt.savefig("display_expcomp_multi_grid_mean_ci.png")

def plot_marginal_dists(
    model: pysb.Model,
    samples: typing.Union[np.ndarray, typing.List[np.ndarray]],
    idxs_mask: list,
    nbins: int = 50,
    nprocs: int = 1,
    save: bool = False,
    ):
    """Plot the marginal distributions for the sampled parameters.
    Adapted from:
    https://github.com/LoLab-VU/JARM/blob/master/model_analysis/pars_dists_plot.py
    """
    
    ndims = len(idxs_mask)
    colors = sns.color_palette(n_colors=ndims)
    columns = 4
    rows = int(ndims / columns)
    if (ndims % columns) > 0:
        rows += 1
    counter = 0

    f, axes = plt.subplots(rows, columns, figsize=(7, 8), sharex=True)
    for r in range(rows):
        for c in range(columns):
            weights = np.ones_like(samples[:, counter])/float(len(samples[:, counter]))
            axes[r, c].hist(samples[:, counter], bins=nbins, color=colors[counter], weights=weights)
            axes[r, c].set_title(model.parameters[idxs_mask[counter]].name, fontdict={'fontsize':8})
            # axes[r, c].set_xlim(-6, 6)
            counter += 1

            if counter >= len(idxs_mask):
                break
    f.add_subplot(111, frameon=False)
    f.subplots_adjust(wspace=0.4)
    f.subplots_adjust(hspace=0.5)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("Log10(Parameter value)", fontsize=14)
    plt.ylabel("Probability", fontsize=14, labelpad=15)
    plt.tight_layout()
    # plt.show()
    if save:
        plt.savefig('pars_dist_plot_{}bins.png'.format(nbins), format='png', bbox_inches="tight")
    
    return

def display_dose_response_comparison(
    model: pysb.Model,
    positions: typing.Union[np.ndarray, typing.List[np.ndarray]],
    idxs_mask: list,
    counts: typing.Union[None, np.ndarray] = None,
    nprocs: int = 1,
    save: bool = False,
    ):
    """Runs the model and plots a comparison with the experimental data for 2AT response."""


    param_values = np.array([param.value for param in model.parameters])
    twoat_mask = [param.name == "TAT_0" for param in model.parameters]

    # Experimental data is in PARM/exp_data
    exp_data = data.training_data()
    times = exp_data["Time"]
    tspan, fidx = util.expand_times(times)
    if counts is None:
        counts = np.array([1]*len(positions))
    new_positions = list()
    for position in positions:
        Y = 10 ** position
        param_values[idxs_mask] = Y
        new_positions.append(param_values.copy())
    #dr_exp = data.kang2019_fig_s3d()
    kang2019_fig_s3d()
    dr_sim = list()    
    for i, key in enumerate(EXP_2AT_CONCS.keys()):
        tat_number = EXP_2AT_CONCS[key] * units.nM_to_molec_per_pL * model.parameters["Vextra"].value
        for position in new_positions:
            position[twoat_mask] = tat_number
        sim_out = util.run_model(model, tspan, param_values=new_positions, nprocs=nprocs)
        fret_sims = list()
        for s in sim_out:
            fret_sims.append(s["FRET"])
        fret_sims = np.array(fret_sims)
        fret_mean = np.average(fret_sims, axis=0, weights=counts)
        fret_std = np.sqrt(
            np.average(np.subtract(fret_sims, fret_mean) ** 2, axis=0, weights=counts)
        )
        idx_max = np.argmax(fret_mean)
        dr_sim.append([EXP_2AT_CONCS[key], fret_mean[idx_max], fret_std[idx_max]])
    dr_sim = np.array(dr_sim)
    plt.errorbar(dr_sim[:,0], dr_sim[:,1], yerr=dr_sim[:,2], color='r', label='Sim.', marker='o', linestyle=' ')
    plt.legend(loc=0, frameon=False)
    plt.tight_layout()    
    if save:
        plt.savefig("display_dose_response_comparison.png")