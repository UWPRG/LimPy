"""
Functions used by the simulate functions for the Langevin Integrator.

This is a Langevin Integrator capable of implementing Molecular Dynamics,
Metadynamics, Well-Tempered Metadynamics (WTMD), and Infrequen Metadynamics.
It has the capability of operating on 1-D or 2-D potentials, if the
potentials are supplied by the user

"""


import numpy as np
import scipy as sp
import pandas as pd


from potential_functions import get_potential_dict, get_boundary_condition_dict


def get_parameters(input_file):
    """Organize input parameters for Langevin Integrator.

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
    inputs = pd.read_csv(inputsfile, comment='#')
    inputs.index = inputs['Parameter']
    inputs = inputs.transpose()
    inputs = inputs.ix[1:]
    dimension = str(inputs['Dimension'][0])
    method = str(inputs['Method'][0])
    filetitle = str(inputs['Data Filename'][0])
    inps = np.zeros(14)
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
    inps[13] = float(inputs['Kb'][0])
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
    """
    Calculate the effective time.

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
    """
    Move walker by one step in 1D via Langevin Integrator.

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


def calc_rmsd(FES, beta, baseline):
    """
    Calculate 3 differed RMSD of the calculated FES and the actual FES.

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
