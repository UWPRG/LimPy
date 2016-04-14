"""This is a script to run the langevin integrator in parallel."""

import langevin_functions as lf

from simulate1D import simulate_1Dsystem
from simulate2D import simulate_2Dsystem
from statistical_functions import perform_ks_analysis

from mpi4py import MPI
import sys
import pandas as pd
import pdb
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

if method == 'Infrequent WT MetaD':
    for i in range(0, 30):
        if dimension == '1-D Potential':
            trial = simulate_1Dsystem(inps, mdps, dimension, method, potfunc,
                                      filetitle, makeplot, plot_freq,
                                      make_movie)
        else:
            trial = simulate_2Dsystem(inps, mdps, dimension, method, potfunc,
                                      filetitle, makeplot, plot_freq,
                                      make_movie)
        if i == 0:
            timedata = pd.DataFrame({'Time': [trial[0]],
                                     'Teff': [trial[1]],
                                     'Event': [trial[3]]})
        else:
            newdata = pd.DataFrame({'Time': [trial[0]],
                                    'Teff': [trial[1]],
                                    'Event': [trial[3]]})
            timedata = timedata.append(newdata, ignore_index=True)

    collected_time_data = comm.gather(timedata, root=0)

    if rank == 0:
        collect = pd.concat(collected_time_data, ignore_index=True)
        collect.reset_index('Time')
        collect.index.name = 'Trial'
	# pdb.set_trace()
	collect = collect[collect['Event'] != 'NULL']
        ks_results = perform_ks_analysis(collect)
	if len(collect[collect['Event'] == 'A']) > 0:
		ks_resultsA = perform_ks_analysis(collect[collect['Event'] == 'A'])
	if len(collect[collect['Event'] == 'B']) > 0:
		ks_resultsB = perform_ks_analysis(collect[collect['Event'] == 'B'])
        with open(filetitle + '_statistics.csv', "ab") as f:
                writer = csv.writer(f)
		writer.writerow(['All Events'])
                writer.writerow([ks_results])
		if len(collect[collect['Event'] == 'A']) > 0:
                	writer.writerow(['A Events'])
                	writer.writerow([ks_resultsA])
                if len(collect[collect['Event'] == 'B']) > 0:
			writer.writerow(['B Events'])
                	writer.writerow([ks_resultsB])
	collect.to_csv(filetitle+'_Allevents.csv', delimiter=',')
else:
        if dimension == '1-D Potential':
            trial = simulate_1Dsystem(inps, mdps,
                                      dimension, method,
                                      potential, filetitle,
                                      makeplot, plot_freq,
                                      make_movie)
            colvar = pd.DataFrame({'CV': trial[0], 'E': trial[1]})
            colvar.reset_index('CV')
        else:
            trial = simulate_2Dsystem(inps, mdps,
                                      dimension, method,
                                      potential, filetitle,
                                      makeplot, plot_freq,
                                      make_movie)
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
