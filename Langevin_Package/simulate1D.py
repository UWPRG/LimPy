"""Simulate a 1D system."""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import pdb
import math

from potential_functions import get_potential_dict, get_boundary_condition_dict
import langevin_functions as lf


def simulate_1Dsystem(inps, mdps, dimension, method, potfunc, filetitle,
                      makeplot, plot_freq, make_movie):
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

    time = np.array([0.0])
    walkerpot = np.array([0.0])

    (pot_dict,_) = get_potential_dict()
    bc_dict = get_boundary_condition_dict()
    try:
        selected_pot = pot_dict[potfunc]
    except KeyError:
        print 'That potential function has not been loaded into the dictionary'
    try:
        selected_bc = bc_dict[potfunc]
    except KeyError:
        print 'That boundary condition has not been loaded into the dictionary'

    baseline = selected_pot(xlong)
    iv = selected_pot(x0)[0]
    pot_base = baseline[0]
    coords[0] = x0
    if makeplot == 'True':
        plt.plot(xlong, pot_base, '-b')
        plt.plot(x0, iv, 'ro', markersize=10)
        plt.axis([xmin, xmax, min(pot_base)-5,
                 max(pot_base)+5])
        plt.xlabel("CV(s)")
        plt.ylabel("F")
        plt.draw()
        plt.pause(0.0001)

    v0 = np.random.normal(0, 1, 1)
    T1 = m*v0**2/kb
    vscaled = v0 *(T/T1)**(0.5)
    p = vscaled * m
    is_periodic = selected_bc(iv, selected_pot(x0)[1], coords[0])[4]

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
            (_,_,triggered,path) = selected_pot(coords[i])

            if triggered is True:

                totaltime = time[i]
                teff = lf.calc_teff(walkerpot, beta, dt)
                return (totaltime, teff, info, path)

        if sp.mod(i, hfreq) == 0 and i > 0:
            if(i == hfreq):
                history[0] = coords[i]
                w[0] = winit
                if is_periodic is True:
                    history = np.append(history, coords[i]+(xmax-xmin))
                    history = np.append(history, coords[i]-(xmax-xmin))
                    w = np.append(w, winit)
                    w = np.append(w, winit)
            else:
                if method == 'Metadynamics':
                    history = np.append(history, coords[i])
                    w = np.append(w, winit)
                    if is_periodic is True:
                        history = np.append(history, coords[i]+(xmax-xmin))
                        history = np.append(history, coords[i]-(xmax-xmin))
                        w = np.append(w, winit)
                        w = np.append(w, winit)
                elif (method == 'Well-Tempered Metadynamics' or
                        method == "Infrequent WT MetaD"):
                    VR = lf.calc_biased_pot(coords[i], history, w, delta,
                                            dimension)
                    history = np.append(history, coords[i])
                    w = np.append(w, winit * np.exp(-VR / (kb*DT)))
                    if is_periodic is True:
                        history = np.append(history, coords[i]+(xmax-xmin))
                        history = np.append(history, coords[i]-(xmax-xmin))
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))
        [pnew, vnew, newcoord, bcbias] = lf.integrate_step(coords[i], history,
                                                           w, delta, DT,
                                                           potfunc, p, m, dt,
                                                           gamma, beta,
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
            bias = np.copy(pot_base)
            for xc in range(0, xlong.size):
                bias[xc] = bias[xc] + lf.calc_biased_pot(xlong[xc], history,
                                                         w, delta, dimension)
            walkv = vnew + lf.calc_biased_pot(coords[i+1], history, w, delta,
                                              dimension)

            plt.clf()
            plt.plot(xlong, bias, '-r')
            plt.plot(xlong, pot_base, '-b')
            plt.plot(coords[i+1], walkv, 'ro', markersize=10)
            plt.axis([xmin, xmax, min(pot_base)-5,
                     max(pot_base)+5])
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
        colvar100 = lf.calc_colvar(coords, history, w, delta, dimension,
                                   xlong, method, beta, T, DT):
        rmsds = lf.calc_rmsd(colvar100[1], beta, pot_base)

        return (coords, colvar100, rmsds, info)

    elif(method == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        path = 'NULL'
        return (totaltime, teff, info, path)


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
