"""These are statistical tests for the Infrequent sampling results."""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp
from scipy import stats
import pandas as pd


def perform_ks_analysis(dataframe):
    """
    Perform the KS Test and determines statistics.

    Parameters:
    -----------
    dataframe   : panda dataframe
                  Contains the effective times of rare events.

    Returns:
    ----------
    statistics  : string
                  Displays the statistical information of the KS test.
    """
    data = dataframe['Teff']
    min = np.min(data)
    max = np.max(data)
    bins = 10*np.size(data)

    time = np.logspace(np.log10(min), np.log10(max), num=bins)
    mu = np.mean(data)

    time_centers = np.r_[0.5 * (time[:-1] + time[1:])]

    hist, bins2 = np.histogram(data, bins=time, density=False)
    cdf = np.cumsum(hist)*1.0/data.size

    # Fit the CDF
    taufit, pcov = curve_fit(analyticalCDF, time_centers, cdf, mu)
    # print "mu (ns)\t\t", mu
    # print "taufit (ns)\t", taufit[0]

    points = 1e5
    randdata = np.random.gamma(1, taufit, np.size(data)*points)

    # perfrom the KS test to see if data points from MetaD are statistically
    # the same as the data points from the analytical fit
    stat, p = ks_2samp(data, randdata)

    statistics = ("mu:" + str(np.mean(data)) + "\n" + "mu_sem:" +
                  str(stats.sem(data)) + "\n" + "sigma:" +
                  str(np.std(data, ddof=1)) + "\n" + "t_m:" +
                  str(np.median(data)) + "\n" + "tau:" + str(taufit) + "\n" +
                  "mu_sigma_ratio:" + str(np.mean(data)/np.std(data, ddof=1)) +
                  "\n" + "log2mu_median_ratio:" +
                  str(np.log(2)*np.mean(data)/np.median(data)) + "\n" +
                  "tau_mu_ratio:" + str(taufit/np.mean(data)) + "\n" +
                  "p-value:" + str(p) + "\n" + "ks-stat:" + str(stat) + "\n" +
                  "events recorded:"+str(np.size(data)))

    return statistics

    # random sampling on data set


def sampling(dataframe, num_iters, sampsize):
    """
    Perform boostrapping procedure for error analysis.

    Parameters:
    -----------
    dataframe   : panda dataframe
                  Contains the effective times of rare events.

    num_iters   : int
                  Number of iterations for bootstrapping

    sampsize    : int
                  Size of bootstrapped samples

    Returns:
    ----------
    means       : array of floats
                  Average effective times of each accepted bootstrapped sample
    pvals       : array of floats
                  P values of each accepted bootstrapped sample
    rejects     : int
                  Number of samples rejected
    """
    # if sampsize > 100
    # sampsize = 100
    data = dataframe['Teff']
    means = 0.0
    pvals = 0.0
    points = 1e4  # number of sampling points for p-val
    alpha = 0.05
    # for i in range((num_iters)):
    # while np.size(means) <= num_iters:
    smalldata = np.random.choice(data, sampsize, replace=True)
    # hist / CDF fit / etc
    min = np.min(smalldata)
    max = np.max(smalldata)
    bins = 10*np.size(smalldata)
    time = np.logspace(np.log10(min), np.log10(max), num=bins)
    mu = np.mean(smalldata)
    time_centers = np.r_[0.5 * (time[:-1] + time[1:])]
    hist, bins2 = np.histogram(smalldata, bins=time, density=False)
    cdf = np.cumsum(hist)*1.0/smalldata.size
    taufit, pcov = curve_fit(analyticalCDF, time_centers, cdf, mu)
    # analysis
    randdata = np.random.gamma(1, taufit, np.size(smalldata)*points)
    stat, p = ks_2samp(smalldata, randdata)
    if p > alpha:
        means = mu
        pvals = p
        reject = 'No'
        # debugprint p, mu
        # means.resize(means.size+1)
        # pvals.resize(pvals.size+1)
    if p < alpha:
        reject = 'Yes'
    # this is just book keeping to remove the last 0 element
    # means = means[:(means.size-1)]
    # pvals = pvals[:(pvals.size-1)]
    return means, pvals, reject


def analyticalCDF(times, tau):
    """Return analytical CDF for a set of data and tau."""
    return 1-np.exp(-times/tau)
    # lets make some plots
    # fig = plt.figure(figsize=(6, 6))
    # fig.subplots_adjust(left=None, bottom=None, right=None, top=None,
    #                      wspace=None,
    #                     hspace=0.2)
    #
    # axes = fig.add_subplot(111)
    # axes.plot(bins2[1:bins], cdf, label='$CDF$')
    # axes.set_xscale('log')
    # axes.plot(time_centers, analyticalCDF(time_centers, taufit),
    #           label='$analytical\ CDF$')
    # first_legend = plt.legend(loc=0)
    # axes.set_xlabel('$log\ time\ (ns)$')
    # axes.set_ylabel('$P_{n\geq1}$')
    # plt.show()
