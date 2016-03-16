"""
This is a Langevin Integrator capable of implementing Molecular Dynamics,
Metadynamics, Well-Tempered Metadynamics (WTMD), and Infrequen Metadynamics.
It has the capability of operating on 1-D or 2-D potentials, if the
potentials are supplied by the user

"""


import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

import pdb
import math
import os

from potential_functions import get_potential_dict


def get_parameters(input_file):
    """Organizes input parameters for Langevin Integrator

    Parameters
    --------------
    input_file : string
                 Filename which contains the desired input parameters


    Returns
    --------
    inps       : Numpy Array
                 Input parameters for intializing system
                 (Steps, step size, X0, Temp, Mass, X min, X max, X increment,
                 Y0,Y min, Y max, Y increment, Gamma)

    mdps       : Numpy Array
                 Metadynamics Parameters
                 (Gaussian Height, Gaussian Width, Deposition Frequency (step),
                 Well Temperature, Trials)

    dimension  : string
                 Dimensionality of system ('1-D Potential' or '2-D Potential')

    method     : string
                 Defines the sampling method (Molecular Dynamics, Metadyamics,
                 Well-Tempered Metadynamics, or Infrequent WT MetaD)

    potfunc    : string
                 Defines the potential function to integrate

    filetitle  :  string
                Name of files for output

    makeplot   : String
                 True or False as to whether make plots

    """

    inputsfile = input_file
    inputs = pd.read_csv(inputsfile+'.csv')

    dimension = str(inputs['Dimension'][0])
    method = str(inputs['Method'][0])
    filetitle = str(inputs['Data Filename'][0])
    inps = np.zeros(16)
    inps[0] = float(inputs['Steps'][0])
    inps[1] = float(inputs['Step size'][0])
    inps[2] = float(inputs['X0'][0])
    inps[3] = float(inputs['Temperature'][0])
    inps[4] = float(inputs['Mass'][0])
    inps[5] = float(inputs['Xmin'][0])
    inps[6] = float(inputs['Xmax'][0])
    inps[7] = float(inputs['Xincrement'][0])
    inps[8] = float(inputs['Y0'][0])
    inps[9] = float(inputs['Ymin'][0])
    inps[10] = float(inputs['Ymax'][0])
    inps[11] = float(inputs['Yincrement'][0])
    inps[12] = float(inputs['Gamma'][0])
    inps[13] = float(inputs['Movie'][0])
    makeplot = float(inputs['Plotting'][0])
    potfunc = str(inputs['Potential_Function'][0])
    mdps = np.zeros(5)
    mdps[0] = float(inputs['Gaussian Height'][0])
    mdps[1] = float(inputs['Gaussian Width'][0])
    mdps[2] = float(inputs['Deposition Frequency'][0])
    mdps[3] = float(inputs['Well Temperature'][0])
    mdps[4] = float(inputs['Trials'][0])

    return (inps, mdps, dimension, method, potfunc, filetitle)


def calc_biased_pot(coords, history, w, delta):
    """
        Calculate the biased force and biased potential.

        Parameters:
        -----------

        coords      :  array of floats
                       X and Y coordinates

        history     : array of floats
                      Locations of where gaussians are deposited

        w           : array of floats
                      Gaussian height

        delta       : float
                      Gaussian width

        Returns:
        --------

        VR          : float (or array of floats)
                      Bias potential

    """

    VR = sum(w*np.exp(-(coords-history)**2 / 2 / delta**2))

    return VR


def calc_biased_force(coords, history, w, delta, DT, base_force):
    """
        Calculate the biased force and biased potential.

        Parameters:
        -----------

        coords      :  array of floats
                       X and Y coordinates

        history     : array of floats
                      Locations of where gaussians are deposited

        w           : array of floats
                      Gaussian height

        delta       : float
                      Gaussian width

        DT          : float
                      Well-Temperature

        base_force  : float
                       Underlying potential energy

        Returns:
        --------

        Fbias       : float (or array of floats)
                      Biased force

    """

    F = sum(w * (coords-history) / delta**2 *
            np.exp(-(coords-history)**2 / 2 / delta**2))
    Fbias = base_force + F

    return Fbias


def calc_teff(walkerpot, beta, dt):
    """Calculates the effective time

    Parameters:
    -----------
    walkerpot: array of floats
               Bias potential of walker at each point in time

    beta:      float
               1/kT

    dt:        float
               timestep

    Returns:
    --------
    teff:      float
               Effective time
    """

    walkercalc = np.delete(walkerpot, 0)
    teff = sum(dt * np.exp(walkercalc * beta))

    return teff


def integrate_1D_step(coords, history, w,  delta, DT, potfunc, p0, m):
    """Move walker by one step in 1D via Langevin Integrator

        Parameters:
        -----------
        coords      :  array of floats
                       X and Y coordinates

        history     : array of floats
                      Locations of where gaussians are deposited

        w           : array of floats
                      Gaussian height

        delta       : float
                      Gaussian width

        DT          : float
                      Well-Temperature

        potfunc     : string
                      potential energy function

        p0           : float
                      momentum

        m           : float
                      Mass

        Returns:
        --------
        pnew        : float
                      new momentum

        new_coords  : float
                      new coordinates

    """

    pot_dict = get_potential_dict
    force = pot_dict[potfunc]

    f = force(coords)[1]
    fbiased = calc_biased_force(coords, history, w, delta, DT, f)

    R1 = sp.rand(1) - 0.5
    R2 = sp.rand(1) - 0.5

    pplus = c1*p0 + c2*R1
    newcoords = coords + (pplus/m) * dt + fbiased/m * ((dt**2) / 2)

    f2 = force(newcoords)[1]
    f2biased = calc_biased_force(coords, history, w, delta, DT, f2)

    pminus = pplus + (fbiased/2 + f2biased/2)*dt
    pnew = c1*pminus + c2*R2

    return [pnew, new_coords]


def recreate_FES(FES, icount, coord, xinc, xmin, xmax, E):
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

            E       : float
                      Energy value to be stored

    """
    index = int(round((round(coord, int(abs(math.log10(xinc)))) +
                      (0-xmin))/xinc))
    if coord > xmin and coord < xmax:
        FES[index] = ((FES[index] * (icount[index]) + E) /
                      (icount[index] + 1))
        icount[index] = icount[index] + 1

    return (FES, icount)


def calc_rmsd(FES, beta, baseline):
    """Calculates 3 differed RMSD of the calculated FES and the actual FES

        Parameters:
        -----------
            FES      : Array of floats
                       Representation of the free energy surface

            beta     : float
                       1/kT

            baseline : Array of floats
                       Underlying potential saved as a grid
    """

    rmsd = np.sqrt(np.sum((((FES-baseline)) * beta)**2) / FES.shape[0])
    rmskld = np.sqrt(np.sum(np.multiply(np.power(((FES - FES.mean()) -
                     (baseline-baseline.mean())) * beta, 2),
                     np.exp((-baseline*beta)))) / np.sum(np.exp((-baseline) *
                                                                beta)))
    rmsalignerr = np.sqrt(np.sum(np.power(((FES-FES.mean()) -
                          (baseline-baseline.mean()))*beta, 2)) /
                          np.size(baseline))
    return np.array([rmsd, rmskld, rmsalignerr])


def simulate_1Dsystem(inps, mdps, dimension, method, potfunc, filetitle,
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

    if(potdim == '1-D Potential'):
        steps = inps[0]
        dt = inps[1]
        x0 = inps[2]
        T = inps[3]
        m = inps[4]
        xmin = inps[5]
        xmax = inps[6]
        xinc = inps[7]

    if (sm == 'Metadynamics'):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        w = np.array([0])
        DT = float("inf")

    if (sm == 'Well-Tempered Metadynamics'):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        DT = mdps[3]
        w = np.array([0])
    if (sm == 'MD'):
        winit = 0
        delta = 1
        hfreq = steps*2
        DT = 10000
        w = np.array([0])
    if (sm == "Infrequent WT MetaD"):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        DT = mdps[3]
        w = np.array([0])

    gamma = inps[12]  # Friction factor
    beta = 1 / T / (1.987E-3)  # units of 1/kcal
    c1 = np.exp(-gamma * dt / 2)  # (Eq.13a)
    c2 = np.sqrt((1 - c1**2) * m / beta)  # (Eq.13b)

    xlong = np.arange(xmin, xmax+xinc, xinc)
    coords = np.empty(int(steps))
    E = np.empty(int(steps))
    history = np.array([0.0])

    time = np.array([0.0])
    walkerpot = np.array([0.0])

    pot_dict = get_potential_dict
    force = pot_dict[potfunc]

    baseline = force(xlong)
    iv = force(x0)[0]
    pot_base = baseline[0]
    coords[0] = x0
    if makeplot == 'True':
        plt.plot(xlong, baseline, '-b')
        plt.plot(x0, v1, 'ro', markersize=10)
        plt.axis([xmin, xmax, xmin-xmin/5, xmax+xmax/5])
        plt.xlabel("CV(s)")
        plt.ylabel("F")
        plt.draw()
        plt.pause(0.0001)

    v0 = sp.rand(1)-0.5
    p = v0 * m

    initial_vals = force(coords[0])

    pot0 = initial_vals[0]
    f0 = initial_vals[1]

    FES = np.zeros_like(xlong)
    icount = np.zeros_like(xlong)

    E[0] = 0.5*p**2 + v0

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

        if (sm == "Infrequent WT MetaD"):
            triggered = force(coords[i])[2]

            if triggered is True:
                totaltime = time[i]
                teff = calc_teff(walkerpot, beta, dt)
                return (totaltime, teff, info)

        if sp.mod(i, hfreq) == 0 and i > 0:

            if(i == hfreq):
                # pdb.set_trace()
                sx[0] = q[i, 0]
                sy[0] = q[i, 1]
                w[0] = winit
            else:
                if sm == 'Metadynamics':
                    s = np.append(s, q[i])
                    w = np.append(w, winit)
                elif (sm == 'Well-Tempered Metadynamics' or
                        sm == "Infrequent WT MetaD"):
                    VR = calc_biased_pot(coords, history, w, delta, DT)
                    s = np.append(s, q[i])
                    w = np.append(w, winit * np.exp(-vr /
                                                    (1.987E-3*DT)))
        [pnew, newcoord] = integrate_1D_step(coords, history, w,  delta, DT,
                                             potfunc, p0, m)
        p = pnew
        coords[i+1] = newcoord

        E[i+1] + 0.5 * p**2 + force(coords[i+1])[0]

        time = np.append(time, dt*(i+1))

        walkerpot = np.append(walkerpot, calc_biased_pot(coords, history,
                                                         w, delta, DT))
        if sm != "Infrequent WT MetaD":
            [FES, icount] = recreate_FES(FES, icount, coord,
                                         xinc, xmin, xmax, E[i+1])
        if plotit == 'True':
            bias = baseline
            for xc in xlong.size:
                bias[xc] = bias[xc] + calc_biased_force(xlong[xc], history,
                                                        w, delta)
            walkv = coords[i+1] + calc_biased_force(coords[i+1], history,
                                                    w, delta)
            plt.clf()
            plt.plot(xlong, bias, '-r')
            plt.plot(xlong, baseline, '-b')
            plt.plot(coords[i+1], walkv, 'ro', markersize=10)
            plt.axis([xmin, xmax, xmin-xmin/5, xmax+xmax/5])
            plt.xlabel("CV(s)")
            plt.ylabel("F")
            plt.draw()
            plt.pause(0.0001)
        i = i + 1

    if(sm != "Infrequent WT MetaD"):
        rmsds = calc_rmsd(FES, beta, baseline)
        return (coords, E, rmsds, info)

    elif(sm == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        return (coords, E, teff, totaltime, info)
