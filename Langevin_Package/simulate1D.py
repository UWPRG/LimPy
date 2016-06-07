"""Simulate a 1D system."""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import pdb
import math

import langevin_functions as lf


def simulate_1Dsystem(inps, mdps, method, potfunc, bcs, filetitle,
                      makeplot, plot_freq, make_movie, ebound):
    """
    Simulatesa walker in a 1D potential.

    Parameters:
    -----------
        inps       : Numpy Array
                     Input parameters for intializing system
                     (Steps, step size, X0, Temp, Mass, X min, X max,
                     X increment, Y0,Y min, Y max, Y increment, Gamma)

        mdps       : Numpy Array
                     Metadynamics Parameters
                     (Gaussian Height, Gaussian Width, Deposition
                     Frequency (step), Well Temperature, Trials)

        dimension  : string
                     Dimensionality of system ('1-D Potential' or
                     '2-D Potential')

        method     : string
                     Defines the sampling method (Molecular Dynamics,
                     Metadyamics, Well-Tempered Metadynamics, or
                     Infrequent WT MetaD)

        potfunc    : string
                     Defines the potential function to integrate

        filetitle  : string
                     Name of files for output

        makeplot   : Boolean
                     If True, make plots, else don't make plots

        plot_freq  : integer
                     Defines how often the plot generated is updated

        make_movie : Boolean
                     If True, save plot images as pngs. If false, then no
                     images are saved.
    Returns:
    --------
        sim_time   : float
                     Simulation time

        teff       : float
                     Effective time from bias

        info       : string
                     Simulation information
        rmsds      : array of floats
                     contains rmsd, rmsd aligned to average, and rmsd with
                     weighting
        coords     : array of floats
                     Coordinates of walker
    """

    steps = inps[0]
    dt = inps[1]
    x0 = inps[2]
    T = inps[3]
    m = inps[4]
    xmin = inps[5]
    xmax = inps[6]
    xinc = inps[7]
    kb = inps[13]
    winit = mdps[0]
    delta = mdps[1]
    hfreq = mdps[2]
    DT = mdps[3]
    w = np.array([0.0])
    bcs = bcs[0]
    dimension = potfunc.dimension
    if (make_movie == 'True'):
        if os.path.exists(filetitle+"_movies"):
            os.rmdir(filetitle+"_movies")
        os.mkdir(filetitle+"_movies")
        frame = 0
        os.chdir(filetitle+"_movies")

    gamma = inps[12]  # Friction factor
    beta = 1 / T / (kb)  # units of 1/kcal

    xlong = np.arange(xmin, xmax+xinc, xinc)
    coords = np.empty(int(steps))
    E = np.empty(int(steps))
    history = np.array([0.0])
    dep_count = 0
    time = np.array([0.0])
    walkerpot = np.array([0.0])

    pot_base = potfunc.get_potential(xlong)
    bias = np.zeros_like(pot_base)
    iv = potfunc.get_potential(x0)
    coords[0] = x0
    if makeplot == 'True':
        plt.plot(xlong, pot_base, '-b')
        plt.plot(x0, iv, 'ro', markersize=10)
        plt.axis([xmin, xmax, ebound[0],
                 ebound[1]])
        plt.xlabel("CV(s)")
        plt.ylabel("F")
        plt.draw()
        plt.pause(0.0001)

    v0 = np.random.normal(0, 1, 1)
    T1 = m*v0**2/kb
    vscaled = v0 * (T/T1)**(0.5)
    p = vscaled * m
    # is_periodic = selected_bc(iv, selected_pot(x0)[1], coords[0])[4]

    FES = np.zeros_like(xlong)
    icount = np.zeros_like(xlong)

    E[0] = 0.5*p**2/m + iv

    info = ('Parameters: \n' + 'Number of steps: ' + str(steps) + '\n' +
            'Initial x coordinate ' + str(x0) + '\n' + 'Temperature ' +
            str(T) + '\n' + 'Timestep ' + str(dt) + '\n' +
            ' Initial Hill Height ' + str(winit) + '\n' +
            'Hill Width ' + str(delta) + '\n' +
            'Deposition Frequency (steps)' + str(hfreq) + '\n' +
            'Well Temperature ' + str(DT) + '\n' + 'Gaussian ' + str(gamma) +
            '\n' + 'Potential ' + str(potfunc))
    i = 0

    while i < steps - 1:

        if (method == "Infrequent WT MetaD"):
            (triggered, path) = potfunc.get_triggered(coords[i])

            if triggered is True:
                totaltime = time[i]
                teff = lf.calc_teff(walkerpot, beta, dt)
                for xc in range(0, xlong.size):
                    bias[xc] = bias[xc]+lf.calc_biased_pot(xlong[xc],
                                                           history[dep_count:],
                                                           w[dep_count:],
                                                           delta, dimension)
                [bins, FES] = lf.calc_FES_1D(coords, bias, xlong, method,
                                             beta, T, DT)

                init_center_point = int((x0 - np.min(bins))/(bins[1]-bins[0]))
                end_center_point = int((coords[i] - np.min(bins)) /
                                       (bins[1]-bins[0]))

                rare_E = FES[end_center_point]
                initial_E = FES[init_center_point]
                for a in range(-100, 100):
                    if (end_center_point + a <= FES.size and
                       end_center_point + a > 0):
                        other_rare_E = FES[end_center_point + a]
                        other_initial_E = FES[init_center_point + a]
                        if other_rare_E > rare_E:
                            rare_E = other_rare_E
                        if other_initial_E < initial_E:
                            initial_E = other_initial_E
                barrier = rare_E - initial_E

                return (totaltime, teff, info, path, barrier)

        if sp.mod(i, hfreq) == 0 and i > 0:
            if(i == hfreq):
                history[0] = coords[i]
                w[0] = winit
                if bcs.type == 'Periodic':
                    history = bcs.add_depositions(coords[i], history)
                    w = np.append(w, winit)
                    w = np.append(w, winit)
            else:

                if method == 'Metadynamics':
                    history = np.append(history, coords[i])
                    w = np.append(w, winit)
                    if bcs.type == 'Periodic':
                        (history) = bcs.add_depositions(coords[i], history)
                        w = np.append(w, winit)
                        w = np.append(w, winit)
                elif (method == 'Well-Tempered Metadynamics' or
                        method == "Infrequent WT MetaD"):
                    VR = lf.calc_biased_pot(coords[i], history, w, delta,
                                            dimension)
                    history = np.append(history, coords[i])
                    w = np.append(w, winit * np.exp(-VR / (kb*DT)))
                    if bcs.type == 'Periodic':
                        (history) = bcs.add_depositions(coords[i], history)
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))
        [pnew, vnew, newcoord, bcbias] = lf.integrate_step(coords[i], history,
                                                           w, delta, DT,
                                                           potfunc, bcs, p, m,
                                                           dt, gamma, beta,
                                                           dimension)
        p = pnew

        coords[i+1] = newcoord
        E[i+1] = 0.5 * p**2/m + vnew

        time = np.append(time, dt*(i+1))

        walkerpot = np.append(walkerpot, lf.calc_biased_pot(coords[i+1],
                                                            history, w, delta,
                                                            dimension)+bcbias)
        if method != "Infrequent WT MetaD":
            [FES, icount] = recreate_1DFES(FES, icount, coords[i+1],
                                           xinc, xmin, xmax, E[i+1])
        if makeplot == 'True' and sp.mod(i, plot_freq) == 0:

            for xc in range(0, xlong.size):

                bias[xc] = bias[xc] + lf.calc_biased_pot(xlong[xc],
                                                         history[dep_count:],
                                                         w[dep_count:],
                                                         delta, dimension)
            walkv = vnew + lf.calc_biased_pot(coords[i+1], history, w, delta,
                                              dimension)
            # pdb.set_trace()

            dep_count = len(history)
            if len(history) == 1:
                dep_count = dep_count-1

            plt.clf()
            plt.plot(xlong, bias+pot_base, '-r')
            plt.plot(xlong, pot_base, '-b')
            plt.plot(coords[i+1], walkv, 'ro', markersize=10)
            plt.axis([xmin, xmax, ebound[0],
                     ebound[1]])
            plt.xlabel("CV(s)")
            plt.ylabel("F")
            plt.draw()
            plt.pause(0.0001)
            if (make_movie == 'True'):
                filename = "movieframe" + str(frame)
                plt.savefig(filename + '.png', bbox_inches='tight')
                frame = frame + 1
        i = i + 1

    if(method != "Infrequent WT MetaD"):
        for xc in range(0, xlong.size):
            bias[xc] = bias[xc] + lf.calc_biased_pot(xlong[xc],
                                                     history[dep_count:],
                                                     w, delta, dimension)

        FES = lf.calc_FES_1D(coords, bias,
                             xlong, method, beta, T, DT)
        rmsds = lf.calc_rmsd(FES[1], beta, pot_base)

        return (coords, E, FES, rmsds, info)

    elif(method == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        path = 'NULL'
        barrier = 0
        return (totaltime, teff, info, path, barrier)


def recreate_1DFES(FES, icount, coord, xinc, xmin, xmax, E):
    """
    Receive and returns an array that recreates the FES.

    Parameters:
    -----------
        FES     : Array of floats
                  Energy values corresponding to x location on x dimension

        icount  : Array of integers
                  Stores number of counts sampled at each location

        coord   : float
                  location of walker

        xinc    : float
                  increment of grid

        xmin    : float
                  minimum value in grid

        xmax    : float
                  maximum value in grid

        E       : float
                  Energy value to be stored

    Returns:
    --------
        FES     : Array of floats
                  Energy values corresponding to x location on x dimension
                  (updated)

        icount  : Array of integers
                  Number of counts sampled at each location (updated)

    """
    index = int(round((round(coord, int(abs(math.log10(xinc)))) +
                      (0-xmin))/xinc))
    if coord > xmin and coord < xmax:
        FES[index] = ((FES[index] * (icount[index]) + E) /
                      (icount[index] + 1))
        icount[index] = icount[index] + 1

    return (FES, icount)
