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
# print 'rank:', rank
# print 'inputdata:', inputdata
# print 'type:', type(inputdata)
inps = inputdata[0]
mdps = inputdata[1]
dimension = inputdata[2]
method = inputdata[3]
potfunc = inputdata[4]
filetitle = inputdata[5]
makeplot = inputdata[6]
trials = mdps[-1]
num_iter = 3
for i in range(0, num_iter):
    if dimension == '1-D Potential':
        trial = simulate_1Dsystem(inps, mdps, dimension, method, potfunc,
                                  filetitle, makeplot)
    else:
        trial = simulate_2Dsystem(inps, mdps, dimension, method, potfunc,
                                  filetitle, makeplot)
    if i == 0:
        timedata = np.array([trial[0], trial[1]])
    else:
        timedata = np.append(timedata, np.array([trial[0], trial[1]]))
# print timedata
collected_time_data = comm.gather(timedata, root=0)

if rank == 0:
    print collected_time_data
    collect = np.asarray(collected_time_data)
    collect = np.reshape(collect, (num_iter*size, 2))
    np.savetxt(filetitle+'_Allevents.csv', collect, delimiter=',')
    ks_results = perform_ks_analysis(filetitle + '_Allevents.csv')

    with open(filetitle + '_statistics.csv', "ab") as f:
            writer = csv.writer(f)
            writer.writerow([ks_results])
