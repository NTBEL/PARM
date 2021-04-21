# NumPy, Pandas, and Matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import data

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
