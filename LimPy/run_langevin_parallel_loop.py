"""This is a script to run the langevin integrator in parallel."""

import langevin_functions as lf

from simulate1D import simulate_1Dsystem
from simulate2D import simulate_2Dsystem
from statistical_functions import perform_ks_analysis
import os
from mpi4py import MPI
import sys
import pandas as pd
import pdb
import numpy as np
import csv

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# if rank == 0:
#    inputsfile = sys.argv[1]
#    inputdata = lf.get_parameters(inputsfile)

# else:
#    inputdata = None
# received = comm.bcast(inputdata, root=0)

inputsfile = sys.argv[1]
received = lf.get_parameters(inputsfile)
inps = received[0]
mdps = received[1]
method = received[2]
potfunc = received[3]
bcs = received[4]
filetitle = received[5]
makeplot = received[6]
plot_freq = received[7]
make_movie = received[8]
ebound = received[9]
np.set_printoptions(threshold=np.nan)
trials = mdps[-1]
checkprogress = 0
mon = 0
print rank
print make_movie
if method == 'Infrequent WT MetaD':
    for j in range(0, 5):
        for i in range(0, 1):
            if potfunc.dimension == '1-D Potential':
                trial = simulate_1Dsystem(inps, mdps, method, potfunc, bcs,
                                          filetitle, makeplot, plot_freq,
                                          make_movie, ebound)
            else:
                trial = simulate_2Dsystem(inps, mdps, method, potfunc, bcs,
                                          filetitle, makeplot, plot_freq,
                                          make_movie, ebound)
            print rank
            if i == 0:
                timedata = pd.DataFrame({'Time': [trial[0]],
                                         'Teff': [trial[1]],
                                         'Event': [trial[3]],
                                         'First Pass': [trial[4]]})
            else:
                newdata = pd.DataFrame({'Time': [trial[0]],
                                        'Teff': [trial[1]],
                                        'Event': [trial[3]],
                                        'First Pass': [trial[4]]})
                timedata = timedata.append(newdata, ignore_index=True)

        collected_time_data = comm.gather(timedata, root=0)

        if rank == 0:
            collect = pd.concat(collected_time_data, ignore_index=True)
            collect.reset_index('Time')
            collect.index.name = 'Trial'
            if os.path.isfile(filetitle+'_Allevents.csv') == False:
                collect.to_csv(filetitle+'_Allevents.csv')
            else:
                old_data = pd.read_csv(filetitle+'_Allevents.csv')
                with open(filetitle+'_Allevents.csv', 'a') as f:
                    collect.to_csv(f, header=False)

    if rank == 0:
        final_collect = pd.read_csv(filetitle+'_Allevents.csv')
        final_collect = final_collect[final_collect['Event'] != 'NULL']
        ks_results = perform_ks_analysis(final_collect)
        bpresent = (final_collect['Event'] == 'B').any()
        apresent = (final_collect['Event'] == 'A').any()
        if apresent == True:
            adata = final_collect[final_collect['Event'] == 'A']
            ks_resultsA = perform_ks_analysis(adata)
        if bpresent == True:
            bdata = final_collect[final_collect['Event'] == 'B']
            ks_resultsB = perform_ks_analysis(bdata)

        with open(filetitle + '_statistics.csv', "ab") as f:
            writer = csv.writer(f)
            writer.writerow(['All Events'])
            writer.writerow([ks_results])
            if apresent == True:
                writer.writerow(['A Events'])
                writer.writerow([ks_resultsA])
            if bpresent == True:
                writer.writerow(['B Events'])
                writer.writerow([ks_resultsB])
else:
    print rank
