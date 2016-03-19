"""
This is a Langevin Integrator capable of implementing Molecular Dynamics,
Metadynamics, Well-Tempered Metadynamics (WTMD), and Infrequen Metadynamics.
It has the capability of operating on 1-D or 2-D potentials, if the
potentials are supplied by the user

"""


import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import pdb
import math
import os

from potential_functions import get_potential_dict, get_boundary_condition_dict


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
    inputs = pd.read_csv(inputsfile)
    inputs.index = inputs['Parameter']
    inputs = inputs.transpose()
    inputs = inputs.ix[1:]
    dimension = str(inputs['Dimension'][0])
    method = str(inputs['Method'][0])
    filetitle = str(inputs['Data Filename'][0])
    inps = np.zeros(13)
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
    makeplot = str((inputs['Plotting'][0]))
    potfunc = str(inputs['Potential_Function'][0])
    mdps = np.zeros(5)
    mdps[0] = float(inputs['Gaussian Height'][0])
    mdps[1] = float(inputs['Gaussian Width'][0])
    mdps[2] = float(inputs['Deposition Frequency'][0])
    mdps[3] = float(inputs['Well Temperature'][0])
    mdps[4] = float(inputs['Trials'][0])

    return (inps, mdps, dimension, method, potfunc, filetitle, makeplot)


def calc_biased_pot(coords, history, w, delta, dimension):
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

        dimension   : string
                      indicates if 1D or 2D potential

        Returns:
        --------

        VR          : float (or array of floats)
                      Bias potential

    """
    if dimension == '1-D Potential':
        VR = sum(w*np.exp(-(coords-history)**2 / 2 / delta**2))
    elif dimension == '2-D Potential':
        if history.size > 2:
            VR = sum(w*np.exp(-(coords[0] - history[:, 0])**2 / 2 / delta**2) *
                     np.exp(-(coords[1] - history[:, 1])**2 / 2 / delta**2))
        else:
            VR = sum(w*np.exp(-(coords[0] - history[0])**2 / 2 / delta**2) *
                     np.exp(-(coords[1] - history[1])**2 / 2 / delta**2))
    return VR


def calc_biased_force(coords, history, w, delta, base_force, dimension):
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

        base_force  : float
                       Underlying potential energy

        dimension   : string
                      indicates if 1D or 2D potential

        Returns:
        --------

        Fbias       : float (or array of floats)
                      Biased force

    """
    if dimension == '1-D Potential':
        F = sum(w * (coords-history) / delta**2 *
                np.exp(-(coords-history)**2 / 2 / delta**2))
        Fbias = base_force + F

    else:
        if history.size > 2:

            Fbiasx = sum(w*(((coords[0]-history[:, 0])/delta**2) *
                         np.exp(-(coords[0]-history[:, 0])**2/2/delta**2) *
                         np.exp(-(coords[1]-history[:, 1])**2/2/delta**2)))
            Fbiasy = sum(w*(((coords[1]-history[:, 1])/delta**2) *
                         np.exp(-(coords[0]-history[:, 0])**2/2/delta**2) *
                         np.exp(-(coords[1]-history[:, 1])**2/2/delta**2)))
            Fbiasx = Fbiasx + base_force[0]
            Fbiasy = Fbiasy + base_force[1]
            Fbias = np.array([Fbiasx, Fbiasy])
        else:
            Fbiasx = sum(w*((coords[0]-history[0])/delta**2) *
                         np.exp(-(coords[0]-history[0])**2/2/delta**2) *
                         np.exp(-(coords[1]-history[1])**2/2/delta**2))
            Fbiasy = sum(w*((coords[1]-history[1])/delta**2) *
                         np.exp(-(coords[0]-history[0])**2/2/delta**2) *
                         np.exp(-(coords[1]-history[1])**2/2/delta**2))
            Fbiasx = Fbiasx + base_force[0]
            Fbiasy = Fbiasy + base_force[1]
            Fbias = np.array([Fbiasx, Fbiasy])
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


def integrate_step(coords, history, w,  delta, DT, potfunc, p0, m, dt,
                   gamma, beta, dimension):
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

        dt          : float
                      time step

        gamma       : float
                       friction factor

        beta        : float
                      1/kT

        dimension  : string
                     Dimensionality of system ('1-D Potential' or
                     '2-D Potential')

        Returns:
        --------
        pnew        : float
                      new momentum

        vnew        : float
                      new potential energy

        new_coords  : float
                      new coordinates

        bcbias      : float
                      bias from boundary condition

    """

    pot_dict = get_potential_dict()
    try:
        force = pot_dict[potfunc]
    except KeyError:
        print 'That potential function has not been loaded into the dictionary'

    bc_dict = get_boundary_condition_dict()
    try:
        apply_bc = bc_dict[potfunc]
    except KeyError:
        print 'That boundary condition has not been loaded into the dictionary'

    c1 = np.exp(-gamma * dt / 2)  # (Eq.13a)
    c2 = np.sqrt((1 - c1**2) * m / beta)  # (Eq.13b)

    if dimension == '1-D Potential':

        f = force(coords)[1]
        fbiased = calc_biased_force(coords, history, w, delta, f, dimension)

        R1 = sp.rand(1) - 0.5
        R2 = sp.rand(1) - 0.5

        pplus = c1*p0 + c2*R1
        newcoords = (coords + (pplus/m) * dt + fbiased/m * ((dt**2) / 2))[0]

        [vnew, f2, _] = force(newcoords)
        [vnew, f2, newcoords, bcbias] = apply_bc(vnew, f2, newcoords)
        f2biased = calc_biased_force(newcoords, history, w, delta, f2,
                                     dimension)
        pminus = pplus + (fbiased/2 + f2biased/2)*dt
        pnew = c1*pminus + c2*R2

    else:
        f = force(coords[0], coords[1])[1]
        [fbiasedx, fbiasedy] = calc_biased_force(coords, history,
                                                 w, delta, f, dimension)
        R1x = sp.rand(1) - 0.5
        R2x = sp.rand(1) - 0.5
        R1y = sp.rand(1) - 0.5
        R2y = sp.rand(1) - 0.5

        # x direction
        pplusx = c1*p0[0] + c2*R1x
        newcoordx = (coords[0] +
                     (pplusx/m) * dt + fbiasedx/m * ((dt**2) / 2))[0]

        # y direction
        pplusy = c1*p0[1] + c2*R1y
        newcoordy = (coords[1] +
                     (pplusy/m) * dt + fbiasedy/m * ((dt**2) / 2))[0]

        [vnew, f2, _] = force(newcoordx, newcoordy)
        newcoords = np.array([newcoordx, newcoordy])
        [vnew, f2, newcoords, bcbias] = apply_bc(vnew, f2, newcoords)
        [f2biasedx, f2biasedy] = calc_biased_force(newcoords, history, w,
                                                   delta, f2, dimension)

        pminusx = pplusx + (fbiasedx/2 + f2biasedx/2)*dt
        pnewx = c1*pminusx + c2*R2x

        pminusy = pplusy + (fbiasedy/2 + f2biasedy/2)*dt
        pnewy = c1*pminusy + c2*R2y

        pnew = np.array([pnewx, pnewy])

    return (pnew, vnew, newcoords, bcbias)


def recreate_1DFES(FES, icount, coord, xinc, xmin, xmax, E):
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

        Returns:
        --------
            rmsds    : Array of floats
                       Holds RMSDs of calculated FES and actual FES
    """

    rmsd = np.sqrt(np.sum((((FES-baseline)) * beta)**2) / FES.shape[0])
    rmskld = np.sqrt(np.sum(np.multiply(np.power(((FES - FES.mean()) -
                     (baseline-baseline.mean())) * beta, 2),
                     np.exp((-baseline*beta)))) / np.sum(np.exp((-baseline) *
                                                                beta)))
    rmsalignerr = np.sqrt(np.sum(np.power(((FES-FES.mean()) -
                          (baseline-baseline.mean()))*beta, 2)) /
                          np.size(baseline))
    rmsds = np.array([rmsd, rmskld, rmsalignerr])

    return rmsds


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
    beta = 1 / T / (1.987E-3)  # units of 1/kcal

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

    bc_dict = get_boundary_condition_dict()
    try:
        apply_bc = bc_dict[potfunc]
    except KeyError:
        print 'That boundary condition has not been loaded into the dictionary'
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
    iv = force(x0, y0)[0]
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
                teff = calc_teff(walkerpot, beta, dt)
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
                    VR = calc_biased_pot(np.array([coords[i, 0],
                                                  coords[i, 1]]), history, w,
                                         delta, dimension)
                    history = np.vstack((history, np.array([coords[i, 0],
                                                           coords[i, 1]])))
                    w = np.append(w, winit * np.exp(-VR / (1.987E-3*DT)))

        [pnew, vnew, newcoord, bcbias] = integrate_step(coords[i], history, w,
                                                        delta, DT, potfunc, p,
                                                        m, dt, gamma, beta,
                                                        dimension)
        p = pnew

        coords[i+1, 0] = newcoord[0]
        coords[i+1, 1] = newcoord[1]
        vnew = force(coords[i+1, 0], coords[i+1, 1])[0]

        E[i+1] = 0.5/m * (p[0]**2 + p[1]**2) + vnew

        time = np.append(time, dt*(i+1))

        walkerpot = np.append(walkerpot,
                              (calc_biased_pot(np.array([coords[i+1, 0],
                                                        coords[i+1, 1]]),
                                               history, w, delta, dimension) +
                               bcbias))
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
                                    calc_biased_pot(np.array([xlong[xc],
                                                             ylong[yc]]),
                                                    history, w, delta,
                                                    dimension))
            walkv = vnew + calc_biased_pot(np.array([coords[i+1, 0],
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
        rmsds = calc_rmsd(FES, beta, pot_base)
        return (coords, E, rmsds, info)

    elif(method == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        return (coords, E, teff, totaltime, info)


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

    steps = inps[0]
    dt = inps[1]
    x0 = inps[2]
    T = inps[3]
    m = inps[4]
    xmin = inps[5]
    xmax = inps[6]
    xinc = inps[7]

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
    beta = 1 / T / (1.987E-3)  # units of 1/kcal
    c1 = np.exp(-gamma * dt / 2)  # (Eq.13a)
    c2 = np.sqrt((1 - c1**2) * m / beta)  # (Eq.13b)

    xlong = np.arange(xmin, xmax+xinc, xinc)
    coords = np.empty(int(steps))
    E = np.empty(int(steps))
    history = np.array([0.0])

    time = np.array([0.0])
    walkerpot = np.array([0.0])

    pot_dict = get_potential_dict()

    try:
        force = pot_dict[potfunc]
    except KeyError:
        print 'That potential function has not been loaded into the dictionary'

    bc_dict = get_boundary_condition_dict()

    try:
        apply_bc = bc_dict[potfunc]
    except KeyError:
        print 'That potential function has not been loaded into the dictionary'

    baseline = force(xlong)
    iv = force(x0)[0]
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

        if (method == "Infrequent WT MetaD"):
            triggered = force(coords[i])[2]

            if triggered is True:
                totaltime = time[i]
                teff = calc_teff(walkerpot, beta, dt)
                return (totaltime, teff, info)

        if sp.mod(i, hfreq) == 0 and i > 0:

            if(i == hfreq):

                history[0] = coords[i]
                w[0] = winit
            else:
                if method == 'Metadynamics':
                    history = np.append(s, coords[i])
                    w = np.append(w, winit)
                elif (method == 'Well-Tempered Metadynamics' or
                        method == "Infrequent WT MetaD"):
                    VR = calc_biased_pot(coords[i], history, w, delta,
                                         dimension)
                    history = np.append(history, coords[i])
                    w = np.append(w, winit * np.exp(-VR / (1.987E-3*DT)))

        [pnew, vnew, newcoord, bcbias] = integrate_step(coords[i], history, w,
                                                        delta, DT, potfunc, p,
                                                        m, dt, gamma, beta,
                                                        dimension)
        p = pnew

        coords[i+1] = newcoord
        E[i+1] = 0.5 * p**2 + vnew

        time = np.append(time, dt*(i+1))

        walkerpot = np.append(walkerpot, calc_biased_pot(coords[i+1], history,
                                                         w, delta, dimension) +
                              bcbias)
        if method != "Infrequent WT MetaD":
            [FES, icount] = recreate_1DFES(FES, icount, coords[i+1],
                                           xinc, xmin, xmax, E[i+1])
        if makeplot == 'True' and sp.mod(i, 1000) == 0:
            bias = np.copy(pot_base)
            for xc in range(0, xlong.size):
                bias[xc] = bias[xc] + calc_biased_pot(xlong[xc], history,
                                                      w, delta, dimension)
            walkv = vnew + calc_biased_pot(coords[i+1], history, w, delta,
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
        rmsds = calc_rmsd(FES, beta, pot_base)
        return (coords, E, rmsds, info)

    elif(method == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        return (totaltime, teff, info)
