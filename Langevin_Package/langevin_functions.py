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
import pdb

#from potential_functions import get_potential_dict, get_boundary_condition_dict
import potential_class
import boundarycondition

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

    method = str(inputs['Method'][0])
    filetitle = str(inputs['Data Filename'][0])
    potfunc = str(inputs['Potential_Function'][0])
    pots = potential_class.get_potential_dict()
    potential= pots[potfunc]
    potfunc = potential()

    if 'Potential Parameters' in inputs.columns:
        potparams = inputs['Potential Parameters']
        potparams = np.fromstring(potparams[0], dtype=float, sep=' ')
        potfunc.set_parameters(potparams)
    bc_dict = boundarycondition.get_boundary_condition_dict()
    if 'X Boundary Condition' in inputs.columns:
        xbc = inputs['X Boundary Condition'][0]
        xbc = bc_dict[xbc]
        xbc = xbc()
        if xbc.type == 'Periodic':
            xbc.set_new_upper_location(float(inputs['Xmax'][0]))
            xbc.set_new_lower_location(float(inputs['Xmin'][0]))
        else:
            if inputs['Upper X BC'][0] == 'True':
                xbc.set_new_upper_location(float(inputs['Xmax'][0]))
            if inputs['Lower X BC'][0] == 'True':
                xbc.set_new_lower_location(float(inputs['Xmin'][0]))
    else:
        xbc = bc_dict['Empty']
        xbc = xbc()
    bcs = [xbc]
    inps = np.zeros(14)
    inps[0] = float(inputs['Steps'][0])
    inps[1] = float(inputs['Step size'][0])
    inps[3] = float(inputs['Temperature'][0])
    inps[4] = float(inputs['Mass'][0])
    inps[12] = float(inputs['Gamma'][0])
    inps[13] = float(inputs['Kb'][0])
    makeplot = str((inputs['Plotting'][0]))
    if potfunc.dimension == '1-D Potential':
            inps[2] = float(inputs['X0'][0])
            inps[5] = float(inputs['Xmin'][0])
            inps[6] = float(inputs['Xmax'][0])
            inps[7] = float(inputs['Xincrement'][0])
            inps[8] = 0.0
            inps[9] = 0.0
            inps[10] = 0.0
            inps[11] = 0.0
    elif potfunc.dimension == '2-D Potential':
            inps[2] = float(inputs['X0'][0])
            inps[5] = float(inputs['Xmin'][0])
            inps[6] = float(inputs['Xmax'][0])
            inps[7] = float(inputs['Xincrement'][0])
            inps[8] = float(inputs['Y0'][0])
            inps[9] = float(inputs['Ymin'][0])
            inps[10] = float(inputs['Ymax'][0])
            inps[11] = float(inputs['Yincrement'][0])
            if 'Y Boundary Condition' in inputs.columns:
                ybc = inputs['Y Boundary Condition'][0]
                ybc = bc_dict[ybc]
                ybc = ybc()
                if ybc.type == 'Periodic':
                    ybc.set_new_upper_location(float(inputs['Ymax'][0]))
                    ybc.set_new_lower_location(float(inputs['Ymin'][0]))
                else:
                    if inputs['Upper Y BC'][0] == 'True':
                        ybc.set_new_upper_location(float(inputs['Ymax'][0]))
                    if inputs['Lower Y BC'][0] == 'True':
                        ybc.set_new_lower_location(float(inputs['Ymin'][0]))

            else:
                ybc = bc_dict['Empty']
                ybc = ybc()
            bcs.append(ybc)
    if makeplot == 'True':
        plot_freq = int((inputs['Plot Freq'][0]))
        plot_emax = int((inputs['Plot Emax'][0]))
        plot_emin = int((inputs['Plot Emin'][0]))
        make_movie = str((inputs['Make Movie'][0]))
    else:
        plot_freq = inps[0]*2.0
        plot_emax = 0.0
        plot_emin = 0.0
        make_movie = 'False'
    ebound = np.array([plot_emin, plot_emax])
    mdps = np.zeros(5)
    if method == 'MD':
        mdps[0] = 0.0
        mdps[1] = 0.01
        mdps[2] = inps[0]*2.0
        mdps[3] = 1.0
        mdps[4] = 1.0
    elif method == 'Metadynamics':
            mdps[0] = float(inputs['Gaussian Height'][0])
            mdps[1] = float(inputs['Gaussian Width'][0])
            mdps[2] = float(inputs['Deposition Frequency'][0])
            mdps[3] = float("inf")
            mdps[4] = 1.0
    elif method =='Well-Tempered Metadynamics':
            mdps[0] = float(inputs['Gaussian Height'][0])
            mdps[1] = float(inputs['Gaussian Width'][0])
            mdps[2] = float(inputs['Deposition Frequency'][0])
            biasfactor = float(inputs['Bias Factor'][0])
            mdps[3] = (biasfactor-1)*inps[3]
            mdps[4] = 1.0
    else:
        mdps[0] = float(inputs['Gaussian Height'][0])
        mdps[1] = float(inputs['Gaussian Width'][0])
        mdps[2] = float(inputs['Deposition Frequency'][0])
        mdps[3] = float(inputs['Well Temperature'][0])
        mdps[4] = float(inputs['Trials'][0])
        if 'Lower Rare Event' in inputs.columns:
            potfunc.set_lower_rare_event(float(inputs['Lower Rare Event']))
        if 'Upper Rare Event' in inputs.columns:
            potfunc.set_upper_rare_event(float(inputs['Upper Rare Event']))

    return (inps, mdps, method, potfunc, bcs, filetitle, makeplot,
            plot_freq, make_movie, ebound)


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


def integrate_step(coords, history, w,  delta, DT, potfunc, bcs, p0, m, dt,
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
    # (pot_dict,_) = get_potential_dict()
    # try:
    #     selected_pot = pot_dict[potfunc]
    # except KeyError:
    #     print 'That potential function has not been loaded into the dictionary'
    #
    # bc_dict = get_boundary_condition_dict()
    # try:
    #     apply_bc = bc_dict[potfunc]
    # except KeyError:
    #     print 'That boundary condition has not been loaded into the dictionary'
    # pdb.set_trace()
    c1 = np.exp(-gamma * dt / 2)  # (Eq.13a)
    c2 = np.sqrt((1 - c1**2) * m / beta)  # (Eq.13b)
    if dimension == '1-D Potential':

        f = potfunc.get_force(coords)
        fbiased = calc_biased_force(coords, history, w, delta, f, dimension)

        R1 = np.random.normal(0, 1, 1)
        R2 = np.random.normal(0, 1, 1)

        pplus = c1*p0 + c2*R1
        newcoords = (coords + (pplus/m) * dt + fbiased/m * ((dt**2) / 2))[0]

        vnew_nobc = potfunc.get_potential(newcoords)
        f2_nobc = potfunc.get_force(newcoords)
        #[vnew, f2, newcoords, bcbias, _] = apply_bc(vnew, f2, newcoords)

        if newcoords < bcs.location[0]:
            vnew = bcs.get_potential(newcoords,
                                     potfunc.get_potential( bcs.location[0]))
            f2 = bcs.get_force(newcoords, f2_nobc)
            newcoords = bcs.get_new_location(newcoords)
            bcbias = bcs.get_bc_bias(coords)
        elif newcoords > bcs.location[1]:
            vnew = bcs.get_potential(newcoords,
                                     potfunc.get_potential( bcs.location[1]))
            f2 = bcs.get_force(newcoords, f2_nobc)
            newcoords = bcs.get_new_location(newcoords)
            bcbias = bcs.get_bc_bias(coords)
        else:
            bcbias = 0
            vnew = vnew_nobc
            f2 = f2_nobc

        f2biased = calc_biased_force(newcoords, history, w, delta, f2,
                                     dimension)
        pminus = pplus + (fbiased/2 + f2biased/2)*dt
        pnew = c1*pminus + c2*R2

    else:
        xbc = bcs[0]
        ybc = bcs[1]

        f = potfunc.get_force(coords)
        [fbiasedx, fbiasedy] = calc_biased_force(coords, history,
                                                 w, delta, f, dimension)
        R1x = np.random.normal(0, 1, 1)
        R2x = np.random.normal(0, 1, 1)
        R1y = np.random.normal(0, 1, 1)
        R2y = np.random.normal(0, 1, 1)

        # x direction
        pplusx = c1*p0[0] + c2*R1x
        newcoordx = (coords[0] +
                     (pplusx/m) * dt + fbiasedx/m * ((dt**2) / 2))[0]

        # y direction
        pplusy = c1*p0[1] + c2*R1y
        newcoordy = (coords[1] +
                     (pplusy/m) * dt + fbiasedy/m * ((dt**2) / 2))[0]

        # [vnew, f2, _, _] = selected_pot(newcoordx, newcoordy)
        newcoords = np.array([newcoordx, newcoordy])
        vnew_nobc = potfunc.get_potential(newcoords)
        f2_nobc = potfunc.get_force(newcoords)

        # [vnew, f2, newcoords, bcbias, _] = apply_bc(vnew, f2, newcoords)
        if newcoords[0] < xbc.location[0]:
            vnew = xbc.get_potential(newcoords[0],
                                     potfunc.get_potential(np.array([xbc.location[0],
                                                                     newcoords[1]])))
            f2x = xbc.get_force(newcoords[0], f2_nobc[0])
            newcoords[0] = xbc.get_new_location(newcoords[0])
            bcbiasx = xbc.get_bc_bias(newcoords[0])
        elif newcoords[0] > xbc.location[1]:
            vnew = bcs.get_potential(newcoords,
                                     potfunc.get_potential(np.array([xbc.location[0],
                                                                    newcoords[1]])))
            f2x = xbc.get_force(newcoords[0], f2_nobc[0])
            newcoords [0]= xbc.get_new_location(newcoords[0])
            bcbiasx = xbc.get_bc_bias(newcoords[0])
        else:
            bcbiasx = 0
            vnew = vnew_nobc
            f2x = f2_nobc[0]
        if newcoords[1] < ybc.location[0]:
            vnew = ybc.get_potential(newcoords[1],
                                     potfunc.get_potential(np.array([newcoords[0],
                                                                     ybc.location[0]])))
            f2y = ybc.get_force(newcoords[1], f2_nobc[1])
            newcoords[1] = ybc.get_new_location(newcoords[1])
            bcbiasy = ybc.get_bc_bias(newcoords[1])
        elif newcoords[1] > ybc.location[1]:
            vnew = ybc.get_potential(newcoords[1],
                                     potfunc.get_potential(np.array([newcoords[0],
                                                                     ybc.location[1]])))
            f2y = ybc.get_force(newcoords[1], f2_nobc[1])
            newcoords[1] = ybc.get_new_location(newcoords[1])
            bcbiasy = ybc.get_bc_bias(newcoords[1])
        else:
            bcbiasy = 0
            vnew = vnew_nobc
            f2y = f2_nobc[1]
        f2=np.array([f2x,f2y])
        [f2biasedx, f2biasedy] = calc_biased_force(newcoords, history, w,
                                                   delta, f2, dimension)

        pminusx = pplusx + (fbiasedx/2 + f2biasedx/2)*dt
        pnewx = c1*pminusx + c2*R2x

        pminusy = pplusy + (fbiasedy/2 + f2biasedy/2)*dt
        pnewy = c1*pminusy + c2*R2y

        pnew = np.array([pnewx, pnewy])
        bcbias = bcbiasx = bcbiasy
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


def calc_colvar_1D(coords, history, w, delta, xlong, method, beta,
                   T, DT):

    vbias = np.zeros_like(xlong)
    if method == 'MD':
        hist, bins = np.histogram(coords, bins=xlong, density=True)
        F = np.log(hist)*-1/beta
    elif method == 'Metadynamics':
        for cv in range(0, xlong.size):
            vbias[cv] = calc_biased_pot(xlong[cv], history, w, delta,
                                        dimension)
        F = -1*vbias
    else:
        for cv in range(0, xlong.size):
            vbias[cv] = calc_biased_pot(xlong[cv], history, w, delta,
                                        dimension)
        F = -1*vbias*(T+DT)/DT
    return (xlong, F)


def calc_colvar_2D(coords, bias, xlong, ylong, method,
                   beta, T, DT):

    if method == 'MD':
        hist, xbins, ybins = np.histogram2d(coords,
                                            bins=[len(xlong), len(ylong)],
                                            density=True)
        F = np.log(hist)*-1/beta
    elif method == 'Metadynamics':
        F = -1*bias
    else:
        F = -1*bias*(T+DT)/DT
    return (xlong, F)
