import langevin_functions as lf
import potential_functions as pf
from statistical_functions import perform_ks_analysis, sampling
import sys
import numpy as np
import pandas as pd
import csv
import pdb

inputsfile = sys.argv[1]
received = lf.get_parameters(inputsfile)
inps = received[0]
mdps = received[1]
dimension = received[2]
method = received[3]
potfunc = received[4]
filetitle = received[5]
makeplot = received[6]

np.set_printoptions(threshold=np.nan)
trials = mdps[-1]
checkprogress = 0
while checkprogress < trials+1:
    if dimension == '1-D Potential':
        trial = lf.simulate_1Dsystem(inps, mdps, dimension, method, potfunc,
                                     filetitle, makeplot)
    else:
        trial = lf.simulate_2Dsystem(inps, mdps, dimension, method, potfunc,
                                     filetitle, makeplot)
    if method == 'Infrequent WT MetaD':
        with open(filetitle + '_Allevents.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow([trial[0], trial[1]])
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
        colvar.to_csv(filetitle+'_COLVAR.csv')
        with open(filetitle + '_info.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow(['RMSD', 'RMSDkld', 'RMSD alignerr'])
                writer.writerow([trial[2]])
        break
with open(filetitle + '_info.csv', "ab") as f:
            writer = csv.writer(f)
            writer.writerow([trial[-1]])

if method == 'Infrequent WT MetaD':
    ks_results = perform_ks_analysis(filetitle + '_Allevents.csv')
    boot_strapped = sampling(filetitle + '_Allevents.csv', 1000,
                             round(trials/2))

    with open(filetitle + '_statistics.csv', "ab") as f:
            writer = csv.writer(f)
            writer.writerow(['Mean Escape Time', 'Mean p-value',
                             '# of Trials Rejected'])
            writer.writerow([boot_strapped])
            writer.writerow([ks_results])
