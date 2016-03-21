"""Simulate a 1D system."""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


import math


from potential_functions import get_potential_dict
import langevin_functions as lf


def simulate_1Dsystem(inps, mdps, dimension, method, potfunc, filetitle,
                      makeplot):
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

        filetitle  :  string
                    Name of files for output

        makeplot   : Boolean
                     If True, make plots, else don't make plots
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
    if (method == 'Metadynamics'):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        w = np.array([0.0])
        DT = float("inf")

    if (method == 'Well-Tempered Metadynamics'):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        DT = mdps[3]
        w = np.array([0.0])
    if (method == 'MD'):
        winit = 0
        delta = 1
        hfreq = steps*2
        DT = 10000
        w = np.array([0.0])
    if (method == "Infrequent WT MetaD"):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        DT = mdps[3]
        w = np.array([0.0])

    gamma = inps[12]  # Friction factor
    beta = 1 / T / (kb)  # units of 1/kcal

    xlong = np.arange(xmin, xmax+xinc, xinc)
    coords = np.empty(int(steps))
    E = np.empty(int(steps))
    history = np.array([0.0])

    time = np.array([0.0])
    walkerpot = np.array([0.0])

    pot_dict = get_potential_dict()

    try:
        selected_pot = pot_dict[potfunc]
    except KeyError:
        print 'That potential function has not been loaded into the dictionary'

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

    v0 = sp.rand(1)-0.5
    p = v0 * m

    FES = np.zeros_like(xlong)
    icount = np.zeros_like(xlong)

    E[0] = 0.5*p**2 + iv

    info = ('Parameters: \n' + 'Number of steps: ' + str(steps) + '\n' +
            'Initial x coordinate ' + str(x0) + '\n' + 'Temperature ' +
            str(T) + '\n' + 'Timestep ' + str(dt) + '\n' +
            ' Initial Hill Height ' + str(winit) + '\n' +
            'Hill Width ' + str(delta) + '\n' +
            'Deposition Frequency (steps)' + str(hfreq) + '\n' +
            'Well Temperature ' + '\n' + str(DT)+'Gaussian ' + str(gamma) +
            '\n' + str(DT)+'Potential ' + str(potfunc))
    i = 0
    while i < steps - 1:

        if (method == "Infrequent WT MetaD"):
            triggered = selected_pot(coords[i])[2]

            if triggered is True:
                totaltime = time[i]
                teff = lf.calc_teff(walkerpot, beta, dt)
                return (totaltime, teff, info)

        if sp.mod(i, hfreq) == 0 and i > 0:

            if(i == hfreq):

                history[0] = coords[i]
                w[0] = winit
            else:
                if method == 'Metadynamics':
                    history = np.append(history, coords[i])
                    w = np.append(w, winit)
                elif (method == 'Well-Tempered Metadynamics' or
                        method == "Infrequent WT MetaD"):
                    VR = lf.calc_biased_pot(coords[i], history, w, delta,
                                            dimension)
                    history = np.append(history, coords[i])
                    w = np.append(w, winit * np.exp(-VR / (kb*DT)))

        [pnew, vnew, newcoord, bcbias] = lf.integrate_step(coords[i], history,
                                                           w, delta, DT,
                                                           potfunc, p, m, dt,
                                                           gamma, beta,
                                                           dimension)
        p = pnew

        coords[i+1] = newcoord
        E[i+1] = 0.5 * p**2 + vnew

        time = np.append(time, dt*(i+1))

        walkerpot = np.append(walkerpot, lf.calc_biased_pot(coords[i+1],
                                                            history, w, delta,
                                                            dimension)+bcbias)
        if method != "Infrequent WT MetaD":
            [FES, icount] = recreate_1DFES(FES, icount, coords[i+1],
                                           xinc, xmin, xmax, E[i+1])
        if makeplot == 'True' and sp.mod(i, 10000) == 0:
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

        i = i + 1

    if(method != "Infrequent WT MetaD"):
        rmsds = lf.calc_rmsd(FES, beta, pot_base)
        return (coords, E, rmsds, info)

    elif(method == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        return (totaltime, teff, info)


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
