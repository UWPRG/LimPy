"""This is a script to run the langevin integrator."""

import langevin_functions as lf

from simulate1D import simulate_1Dsystem
from simulate2D import simulate_2Dsystem
from statistical_functions import perform_ks_analysis, sampling
import sys
import os
import numpy as np
import pandas as pd
import csv

inputsfile = sys.argv[1]
received = lf.get_parameters(inputsfile)
inps = received[0]
mdps = received[1]
dimension = received[2]
method = received[3]
potfunc = received[4]
filetitle = received[5]
makeplot = received[6]
plot_freq = received[7]
make_movie = received[8]
np.set_printoptions(threshold=np.nan)
trials = mdps[-1]
checkprogress = 0
while checkprogress < trials+1:
    if dimension == '1-D Potential':
        trial = simulate_1Dsystem(inps, mdps, dimension, method, potfunc,
                                  filetitle, makeplot, plot_freq, make_movie)
    else:
        trial = simulate_2Dsystem(inps, mdps, dimension, method, potfunc,
                                  filetitle, makeplot, plot_freq, make_movie)
    if method == 'Infrequent WT MetaD':
        with open(filetitle + '_Allevents.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow([trial[3], trial[0], trial[1],])
        data = pd.read_csv(filetitle + '_Allevents.csv')
        checkprogress = len(data) + 2  # (1 to discount header, other the +1)
    else:

        if dimension == '1-D Potential':
            colvar = pd.DataFrame({'CV': trial[0], 'E': trial[1]})
            colvar.reset_index('CV')
        else:
            colvar = pd.DataFrame({'CV1': trial[0][:, 0],
                                   'CV2': trial[0][:, 1],
                                   'E': trial[1]})
            colvar.reset_index('CV1')
        colvar.index.name = 'Step'
        colvar.to_csv(filetitle+'_COLVAR.csv')
        with open(filetitle + '_info.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow(['RMSD', 'RMSDkld', 'RMSD alignerr'])
                writer.writerow([trial[2]])
                writer.writerow([trial[3]])
        break
if os.path.isfile(filetitle + '_info.csv') is False:
    with open(filetitle + '_info.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow([trial[-1]])

if method == 'Infrequent WT MetaD':
    ks_results = perform_ks_analysis(filetitle + '_Allevents.csv')
    monitor = 0
    if os.path.isfile('bootstrapped.csv') is False:
        with open('bootstrapped.csv', "ab") as f:
                writer = csv.writer(f)
                writer = writer.writerow(['Means', 'Pvals', 'Rejected'])
    while monitor <= 1000:
        (means, pvals, reject) = sampling(filetitle + '_Allevents.csv', 1000,
                                          round(trials/2))
        with open('bootstrapped.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow([means, pvals, reject])
        checkprogress = pd.read_csv('bootstrapped.csv')
        checkaccept = checkprogress[checkprogress['Rejected'] == 'No']
        if ((len(checkprogress) - len(checkaccept))/len(checkprogress) >
           0.90 and len(checkprogress) > 100):
            break
        monitor = len(checkaccept)

    finisheddata = pd.read_csv('bootstrapped.csv')
    validdata = finisheddata[finisheddata['Rejected'] == 'No']
    rejectedtrials = (len(finisheddata) - len(validdata))
    if os.path.isfile(filetitle + '_statistics.csv') is False:
        with open(filetitle + '_statistics.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow(['Bootstrapped: Mean Escape Time',
                                 'Mean p-value', '# of Trials Rejected'])
                writer.writerow([validdata['Means'].mean(),
                                 validdata['Pvals'].mean(),
                                 rejectedtrials])
                writer.writerow([ks_results])
