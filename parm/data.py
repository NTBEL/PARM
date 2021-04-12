"""This module defines functions to read/load experimental data used for training and analysis.
"""
import os

import numpy as np
import pandas as pd

PATH = os.path.dirname(os.path.abspath(__file__))

def training_data():
    """Loads the FRET signals at different 2AT concentrations used for model training.

    Returns:
        pandas.DataFrame
    """
    return kang2019_fig_s3c()

def kang2019_fig_s3c():
    """Loads the FRET signals at different 2AT concentrations from Fig. S3C of Kang et al. 2019.

    Note that the data has been clipped and the times shifted so that the point
    where 2AT was added (180 s) is now the zero timepoint.

    Returns:
        pandas.DataFrame
    """
    fpath = os.path.abspath(os.path.join(PATH, '../exp_data/FRET_data_Kang2019_FigS3C.csv'))
    df = pd.read_csv(fpath)
    return df

def kang2019_fig_2d():
    """Loads the FRET signals from Fig. 2D of Kang et al. 2019.

    Returns:
        pandas.DataFrame
    """
    fpath = os.path.abspath(os.path.join(PATH, '../exp_data/FRET_data_Kang2019_Fig2D.csv'))
    df = pd.read_csv(fpath)
    return df

def kang2019_fig_2f_raw():
    """Loads the raw FRET signals corresponding to Fig. 2F of Kang et al. 2019.

    Returns:
        pandas.DataFrame
    """
    fpath = os.path.abspath(os.path.join(PATH, '../exp_data/raw_FRETratio_data_Kang2019_Fig2F.csv'))
    df = pd.read_csv(fpath, dtype={'LaserEnergy':np.int, 'PulseNumber':np.int})
    # Go ahead and compute Average and Standard Deviation over the samples:
    #   Sample_0 to Sample_7.
    df = df.assign(SampleAverage=df.iloc[:,3:].mean(axis=1))
    df = df.assign(SampleStd=df.iloc[:,3:].std(axis=1))
    return df

def kang2019_fig_2f():
    """Recreates the data used in Fig. 2F of Kang et al. 2019.

    Returns:
        pandas.DataFrame
    """
    df_raw = kang2019_fig_2f_raw()
    pulse_ns = [1,5,10]
    fluences = [20,60,100]
    # Control, which is used for normalization
    norm_val = None
    df_raw_0 = df_raw[(df_raw.LaserEnergy == 0) & (df_raw.PulseNumber == 0)]
    df_raw_0_t = df_raw_0[df_raw_0.Time > 180]
    max_vals = list()
    for s in range(8):
        sample = "Sample_{}".format(s)
        s_avg = df_raw_0_t[sample]
        s_avg_max = s_avg.max()
        max_vals.append(s_avg_max)
    norm_val = np.mean(np.array(max_vals))
    norm_val_err = np.std(np.array(max_vals))
    # Laser cases
    norm_drr = list()
    for pn in pulse_ns:
        for f in fluences:
            df_raw_fn = df_raw[(df_raw.LaserEnergy == f) & (df_raw.PulseNumber == pn)]
            df_raw_fn_t = df_raw_fn[df_raw_fn.Time > 180]
            max_vals = list()
            for s in range(8):
                sample = "Sample_{}".format(s)
                s_avg = df_raw_fn_t[sample]
                s_avg_max = s_avg.max()
                max_vals.append(s_avg_max)
            n_drr = np.mean(np.array(max_vals)/norm_val)
            # Propagate from deviation of the peak values and deviation of the control.
            n_drr_err = np.sqrt((np.std(np.array(max_vals)/norm_val)/n_drr)**2 + (norm_val_err/norm_val)**2) * n_drr
            norm_drr.append(dict({'LaserEnergy':f, 'PulseNumber':pn, 'NormPeakdRR':n_drr, 'NormPeakdRR_err':n_drr_err}))

    return pd.DataFrame(norm_drr)

def kang2019_fig_s3d():
    df = training_data()
    concs = ['10', '31.6', '100', '316', '1000', '3160']
    dat = list()
    for c in concs:
        idx_max = df[c+'_sig'].argmax()
        max_dr = df[c+'_sig'][idx_max]
        max_dr_err = df[c+'_err'][idx_max]
        #print(idx_max, max_dr, max_dr_err)
        dat.append(dict({'[2AT]':float(c), 'Max.dRR':max_dr, 'Max.dRR_err':max_dr_err}))
    dat = pd.DataFrame(dat)
    return dat
