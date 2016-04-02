"""This is a script to run the langevin integrator in parallel."""

import langevin_functions as lf

from simulate1D import simulate_1Dsystem
from simulate2D import simulate_2Dsystem
from statistical_functions import perform_ks_analysis

from mpi4py import MPI
import sys
# import os
import numpy as np
import csv

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
    inputsfile = sys.argv[1]
    received = lf.get_parameters(inputsfile)

else:
    received = None
inputdata = comm.bcast(received, root=0)

inps = inputdata[0]
mdps = inputdata[1]
dimension = inputdata[2]
method = inputdata[3]
potfunc = inputdata[4]
filetitle = inputdata[5]
makeplot = inputdata[6]
plot_freq = inputdata[7]
make_movie = inputdata[8]
trials = mdps[-1]
num_iter = 3
if potfunc == 'Infrequent WT MetaD':
    for i in range(0, num_iter):
        if dimension == '1-D Potential':
            trial = simulate_1Dsystem(inps, mdps, dimension, method, potfunc,
                                      filetitle, makeplot, plot_freq, make_movie)
        else:
            trial = simulate_2Dsystem(inps, mdps, dimension, method, potfunc,
                                      filetitle, makeplot, plot_freq, make_movie)
        if i == 0:
            timedata = pd.DataFrame({'Time': [trial[0]],
                                     'Teff': [trial[1]],
                                     'Event': [trial[3]]})
        else:
            newdata = pd.DataFrame({'Time': [trial[0]],
                                     'Teff': [trial[1]],
                                     'Event': [trial[3]]})
    # print timedata
    collected_time_data = comm.gather(timedata, root=0)

    if rank == 0:
        collect = np.asarray(collected_time_data)
        collect = np.reshape(collect, (num_iter*size, 3))
        np.savetxt(filetitle+'_Allevents.csv', collect, delimiter=',')
        ks_results = perform_ks_analysis(filetitle + '_Allevents.csv')

        with open(filetitle + '_statistics.csv', "ab") as f:
                writer = csv.writer(f)
                writer.writerow([ks_results])
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
        if rank == 0:
            colvar.to_csv(filetitle+'_COLVAR.csv')
            with open(filetitle + '_info.csv', "ab") as f:
                    writer = csv.writer(f)
                    writer.writerow(['RMSD', 'RMSDkld', 'RMSD alignerr'])
                    writer.writerow([trial[2]])
                    writer.writerow([trial[3]])
