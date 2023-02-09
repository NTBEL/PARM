import typing
# NumPy, Pandas, and Matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysb

from . import data

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
    plt.errorbar(figs3d['[2AT]'], figs3d['Max.dRR'], figs3d['Max.dRR_err'], marker='s', linestyle="", color='black', capsize=2.5)
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
    tspan, fidx = expand_times(times)
    kang2019_fig_s3c()

    positions = list()
    for tat_conc in EXP_2AT_CONCS.values():
        param_values[twoat_mask] = tat_conc * units.nM_to_molec_per_pL * model.parameters["Vextra"].value
        positions.append(param_values.copy())
    sim_out = util.run_model(model, tspan, param_values=positions, nprocs=nprocs)
    for i, key in enumerate(EXP_2AT_CONCS.keys()):
        ysim = sim[i]["FRET"]
        plt.plot(tspan, ysim, label=None, color=EXP_PLT_COLORS[key])
    if save:
        plt.savefig("display_expcomp_single.png")

def display_expcomp_multi_grid_mean_ci(
    model: pysb.Model,
    positions: : typing.Union[np.ndarray, typing.List[np.ndarray]],
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
    tspan, fidx = expand_times(times)
    if counts is None:
        counts = np.array([1]*len(positions))
    f, axes = plt.subplots(2, 3, figsize=(4, 7), sharex=True, sharey=True)
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
        for i, s in enumerate(sim_out):
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
        axes[ax_row, ax_col].errorbar(times, exp_data[col_sig], exp_data[col_err],
                     capsize=2.5, label=key, color=EXP_PLT_COLORS[key],
                     marker=EXP_PLT_MARKERS[key])

        axes[ax_row, ax_col].plot(tspan, fret_mean, color=EXP_PLT_COLORS[key],
                        linewidth=2, label="Sim.")
        axes[ax_row, ax_col].fill_between(
            tspan,
            fret_mean - fret_std,
            fret_mean + fret_std,
            alpha=0.25,
            color=EXP_PLT_COLORS[key],
        )
    if save:
        plt.savefig("display_expcomp_multi_grid_mean_ci.png")
