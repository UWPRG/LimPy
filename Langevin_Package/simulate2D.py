import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import pdb
import math
import os

from potential_functions import get_potential_dict, get_boundary_condition_dict
import langevin_functions as lf


def simulate_2Dsystem(inps, mdps, dimension, method, potfunc, filetitle,
                      makeplot):
    """Simulates a walker in a 1D potential

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
    y0 = inps[8]
    ymin = inps[9]
    ymax = inps[10]
    yinc = inps[11]
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
    ylong = np.arange(ymin, ymax+yinc, yinc)

    history = np.array([0.0, 0.0])
    coords = np.empty([int(steps), 2])
    E = np.empty(int(steps))

    time = np.array([0.0])
    walkerpot = np.array([0.0])

    pot_dict = get_potential_dict()
    try:
        force = pot_dict[potfunc]
    except KeyError:
        print 'That potential function has not been loaded into the dictionary'

    baseline = force(xlong, ylong)
    iv = force(x0, y0)[0]
    pot_base = baseline[0]
    FES = np.zeros_like(pot_base)
    icount = np.zeros_like(FES)
    coords[0, 0] = x0
    coords[0, 1] = y0
    cmap = plt.cm.PRGn
    levels = np.arange(np.min(pot_base), np.max(pot_base)+0.25, 0.25)
    v0x = sp.rand(1) - 0.5
    v0y = sp.rand(1) - 0.5
    px = v0x * m
    py = v0y * m
    p = np.array([px, py])
    E[0] = 0.5 * (px**2 + py**2) + iv

    if makeplot == 'True':
        plt.clf()
        plt.ion()
        plt.subplot(221)

        cset1 = plt.contourf(xlong, ylong, pot_base, levels,
                             cmap=plt.cm.get_cmap(cmap, levels.size - 1))
        plt.colorbar(cset1)
        plt.title('2-D Potential')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.scatter(x0, y0, marker='o', color='r', zorder=10)

        plt.subplot(222)
        cset2 = plt.contourf(xlong, ylong, pot_base, levels,
                             cmap=plt.cm.get_cmap(cmap, levels.size - 1))
        plt.colorbar(cset2)
        plt.draw()
        plt.pause(0.0001)

    info = ('Parameters: \n' + 'Number of steps: ' + str(steps) + '\n' +
            'Initial x coordinate ' + str(x0) + 'Initial y coordinate ' +
            str(y0) + '\n' + 'Temperature ' +
            str(T) + '\n' + 'Timestep ' + str(dt) + '\n' +
            ' Initial Hill Height ' + str(winit) + '\n' +
            'Hill Width ' + str(delta) + '\n' +
            'Deposition Frequency (steps)' + str(hfreq) + '\n' +
            'Well Temperature ' + '\n' + str(DT)+'Gaussian ' + str(gamma) +
            '\n' + str(DT)+'Potential ' + str(potfunc))
    i = 0
    while i < steps - 1:

        if (method == "Infrequent WT MetaD"):
            triggered = force(coords[i, 0], coords[i, 1])[2]
            if triggered is True:
                totaltime = time[i]
                teff = lf.calc_teff(walkerpot, beta, dt)
                return (totaltime, teff, info)

        if sp.mod(i, hfreq) == 0 and i > 0:

            if(i == hfreq):

                history[0] = coords[i, 0]
                history[1] = coords[i, 1]
                w[0] = winit

            else:

                if method == 'Metadynamics':
                    w = np.append(w, winit)
                    history = np.vstack((history, np.array([coords[i, 0],
                                                           coords[i, 1]])))

                elif (method == 'Well-Tempered Metadynamics' or
                        method == "Infrequent WT MetaD"):
                    VR = lf.calc_biased_pot(np.array([coords[i, 0],
                                                      coords[i, 1]]), history,
                                            w, delta, dimension)
                    history = np.vstack((history, np.array([coords[i, 0],
                                                           coords[i, 1]])))
                    w = np.append(w, winit * np.exp(-VR / (kb*DT)))

        [pnew, vnew, newcoord, bcbias] = lf.integrate_step(coords[i], history,
                                                           w, delta, DT,
                                                           potfunc, p, m, dt,
                                                           gamma, beta,
                                                           dimension)
        p = pnew

        coords[i+1, 0] = newcoord[0]
        coords[i+1, 1] = newcoord[1]
        vnew = force(coords[i+1, 0], coords[i+1, 1])[0]

        E[i+1] = 0.5/m * (p[0]**2 + p[1]**2) + vnew

        time = np.append(time, dt*(i+1))

        walkerpot = np.append(walkerpot,
                              (lf.calc_biased_pot(np.array([coords[i+1, 0],
                                                            coords[i+1, 1]]),
                                                  history, w, delta,
                                                  dimension) + bcbias))
        if method != "Infrequent WT MetaD":
            [FES, icount] = recreate_2DFES(FES, icount,
                                           np.array([coords[i+1, 0],
                                                     coords[i+1, 1]]),
                                           xinc, xmin, xmax,
                                           yinc, ymin, ymax, E[i+1])
        if makeplot == 'True' and sp.mod(i, 1000) == 0:
            bias = np.copy(pot_base)
            for yc in range(0, ylong.size):
                for xc in range(0, xlong.size):
                    bias[yc, xc] = (bias[yc, xc] +
                                    lf.calc_biased_pot(np.array([xlong[xc],
                                                                 ylong[yc]]),
                                                       history, w, delta,
                                                       dimension))
            walkv = vnew + lf.calc_biased_pot(np.array([coords[i+1, 0],
                                                        coords[i+1, 1]]),
                                              history, w, delta, dimension)
            plt.clf()
            plt.subplot(221)
            cset2 = plt.contourf(xlong, ylong, pot_base, levels,
                                 cmap=plt.cm.get_cmap(cmap, levels.size - 1))
            plt.colorbar(cset2)
            plt.scatter(coords[i+1, 0], coords[i+1, 1], marker='o',
                        color='r', zorder=10)
            plt.xlabel("CV1")
            plt.ylabel("CV2")
            plt.subplot(222)
            cset3 = plt.contourf(xlong, ylong, bias, levels,
                                 cmap=plt.cm.get_cmap(cmap, levels.size - 1))
            plt.colorbar(cset3)
            plt.scatter(coords[i+1, 0], coords[i+1, 1], marker='o',
                        color='r', zorder=10)
            plt.xlabel("CV1")
            plt.ylabel("CV2")
            plt.subplot(223)
            cset4 = plt.contourf(xlong, ylong, bias-pot_base, levels,
                                 cmap=plt.cm.get_cmap(cmap, levels.size - 1))
            plt.colorbar(cset4)
            plt.scatter(coords[i+1, 0], coords[i+1, 1], marker='o',
                        color='r', zorder=10)
            plt.xlabel("CV1")
            plt.ylabel("CV2")
            plt.draw()
            plt.pause(0.0001)
        i = i + 1

    if(method != "Infrequent WT MetaD"):
        editFES = np.array(0)
        editvcalc = np.array(0)
        for k in range(0, ylong.size):
            for j in range(0, xlong.size):
                if(icount[k, j] > 0):
                    editFES = np.append(editFES, FES[k, j])
                    editvcalc = np.append(editvcalc, pot_base[k, j])
        rmsds = lf.calc_rmsd(FES, beta, pot_base)
        return (coords, E, rmsds, info)

    elif(method == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        return (coords, E, teff, totaltime, info)


def recreate_2DFES(FES, icount, coords, xinc, xmin, xmax, yinc, ymin, ymax, E):
    """Receives and returns an array that recreates the FES

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

            yinc    : float
                      increment of grid

            ymin    : float
                      minimum value in grid

            ymax    : float
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
    xindex = int(round((round(coords[0],
                        int(abs(math.log10(xinc)))) +
                        (0 - xmin)) / xinc))

    yindex = int(round((round(coords[1],
                        int(abs(math.log10(yinc)))) +
                        (0 - ymin)) / yinc))
    if (coords[0] > xmin and coords[0] < xmax and
            coords[1] > ymin and coords[1] < ymax):

        FES[yindex, xindex] = ((FES[yindex, xindex] *
                               (icount[yindex, xindex]) + E) /
                               (icount[yindex, xindex] + 1))
        icount[yindex, xindex] = icount[yindex, xindex] + 1

    return (FES, icount)