"""This is a script to run the langevin integrator."""

import langevin_functions as lf

from simulate1D import simulate_1Dsystem
from simulate2D import simulate_2Dsystem
from statistical_functions import perform_ks_analysis, sampling
import sys
import os
import numpy as np
import pandas as pd
import pdb
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

while checkprogress < trials:
    if dimension == '1-D Potential':
        trial = simulate_1Dsystem(inps, mdps, dimension, method, potfunc,
                                  filetitle, makeplot, plot_freq, make_movie)
    else:
        trial = simulate_2Dsystem(inps, mdps, dimension, method, potfunc,
                                  filetitle, makeplot, plot_freq, make_movie)

    if method == 'Infrequent WT MetaD':
        if checkprogress == 0:
            timedata = pd.DataFrame({'Time': [trial[0]],
                                     'Teff': [trial[1]],
                                     'Event': [trial[3]]})
            checkprogress = len(timedata)
        else:
            newdata = pd.DataFrame({'Time': [trial[0]],
                                    'Teff': [trial[1]],
                                    'Event': [trial[3]]})
            timedata = timedata.append(newdata, ignore_index=True)
            checkprogress = len(timedata)
    else:

        if dimension == '1-D Potential':
            colvar = pd.DataFrame({'CV': trial[1][0], 'E': trial[1][1]})
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
                writer.writerow([trial[-2]])

if method == 'Infrequent WT MetaD':
    timedata.to_csv(filetitle + '_Allevents.csv', delimiter=',')
    ks_results = perform_ks_analysis(timedata)
    if len(timedata[timedata['Event'] == 'A']) > 0:
        ks_resultsA = perform_ks_analysis(timedata[timedata['Event'] == 'A'])
    if len(timedata[timedata['Event'] == 'B']) > 0:
        ks_resultsB = perform_ks_analysis(timedata[timedata['Event'] == 'B'])
    monitor = 0
    # if os.path.isfile('bootstrapped.csv') is False:
    #     with open('bootstrapped.csv', "ab") as f:
    #             writer = csv.writer(f)
    #             writer = writer.writerow(['Means', 'Pvals', 'Rejected'])
    # while monitor <= 1000:
    #     (means, pvals, reject) = sampling(filetitle + '_Allevents.csv', 1000,
    #                                       round(len(timedata)/2))
    #     with open('bootstrapped.csv', "ab") as f:
    #             writer = csv.writer(f)
    #             writer.writerow([means, pvals, reject])
    #     checkprogress = pd.read_csv('bootstrapped.csv')
    #     checkaccept = checkprogress[checkprogress['Rejected'] == 'No']
    #     if ((len(checkprogress) - len(checkaccept))/len(checkprogress) >
    #        0.90 and len(checkprogress) > 100):
    #         break
    #     monitor = len(checkaccept)
    #
    # finisheddata = pd.read_csv('bootstrapped.csv')
    # validdata = finisheddata[finisheddata['Rejected'] == 'No']
    # rejectedtrials = (len(finisheddata) - len(validdata))
    if os.path.isfile(filetitle + '_statistics.csv') is False:
        with open(filetitle + '_statistics.csv', "ab") as f:
                writer = csv.writer(f)
                # writer.writerow(['Bootstrapped: Mean Escape Time',
                #                  'Mean p-value', '# of Trials Rejected'])
                # writer.writerow([validdata['Means'].mean(),
                #                  validdata['Pvals'].mean(),
                #                  rejectedtrials])
                writer.writerow(['All events'])
                writer.writerow([ks_results])
                if len(timedata[timedata['Event'] == 'A']) > 0:
                    writer.writerow(['A events'])
                    writer.writerow([ks_resultsA])
                if len(timedata[timedata['Event'] == 'B']) > 0:
                    writer.writerow(['B events'])
                    writer.writerow([ks_resultsB])
