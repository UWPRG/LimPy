"""Simulate a 1D system."""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os
import math
import pdb

import langevin_functions as lf


def simulate_2Dsystem(inps, mdps, method, potfunc, bcs, filetitle,
                      makeplot, plot_freq, make_movie,
                      ebound):
    """
    Simulate a walker in a 2D potential.

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
    y0 = inps[8]
    ymin = inps[9]
    ymax = inps[10]
    yinc = inps[11]
    kb = inps[13]
    winit = mdps[0]
    delta = mdps[1]
    hfreq = mdps[2]
    DT = mdps[3]
    w = np.array([0.0])
    xbc = bcs[0]
    ybc = bcs[1]
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
    ylong = np.arange(ymin, ymax+yinc, yinc)

    history = np.array([0.0, 0.0])
    coords = np.empty([int(steps), 2])
    E = np.empty(int(steps))

    time = np.array([0.0])
    walkerpot = np.array([0.0])

    # (pot_dict,_) = get_potential_dict()
    # bc_dict = get_boundary_condition_dict()
    # try:
    #     selected_pot = pot_dict[potfunc]
    # except KeyError:
    #     print 'That potential function has not been loaded into the dictionary'
    # try:
    #     selected_bc = bc_dict[potfunc]
    # except KeyError:
    #     print 'That potential function has not been loaded into the dictionary'

    iv = potfunc.get_potential(np.array([x0, y0]))
    pot_base = potfunc.get_potential(np.array([xlong, ylong]))
    bias = np.zeros_like(pot_base)
    FES = np.zeros_like(pot_base)
    icount = np.zeros_like(FES)
    coords[0, 0] = x0
    coords[0, 1] = y0
    cmap = cmap = plt.cm.jet
    levels = np.arange(ebound[0],ebound[1],(ebound[1]-ebound[0])/10)
    v0x = np.random.normal(0, 1, 1)
    v0y = np.random.normal(0, 1, 1)
    T1x = m*v0x**2/kb
    T1y = m*v0y**2/kb
    vscaledx = v0x *(T/T1x)**(0.5)
    vscaledy = v0y *(T/T1y)**(0.5)
    px = vscaledx * m
    py = vscaledy * m
    p = np.array([px, py])
    E[0] = 0.5 * (px**2 + py**2) + iv

    if makeplot == 'True':
        plt.clf()
        plt.ion()
        # plt.subplot(221)

        cset1 = plt.contourf(xlong, ylong, pot_base, levels,
                             cmap=plt.cm.get_cmap(cmap, levels.size - 1))
        plt.colorbar(cset1)
        plt.title('2-D Potential')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.scatter(x0, y0, marker='o', color='r', zorder=10)
        #
        # plt.subplot(222)
        # cset2 = plt.contourf(xlong, ylong, pot_base, levels,
        #                      cmap=plt.cm.get_cmap(cmap, levels.size - 1))
        # plt.colorbar(cset2)
        plt.draw()
        plt.pause(0.000001)

    info = ('Parameters: \n' + 'Number of steps: ' + str(steps) + '\n' +
            'Initial x coordinate ' + str(x0) + 'Initial y coordinate ' +
            str(y0) + '\n' + 'Temperature ' +
            str(T) + '\n' + 'Timestep ' + str(dt) + '\n' +
            ' Initial Hill Height ' + str(winit) + '\n' +
            'Hill Width ' + str(delta) + '\n' +
            'Deposition Frequency (steps)' + str(hfreq) + '\n' +
            'Well Temperature ' + str(DT) + '\n' + 'Gaussian ' + str(gamma) +
            '\n' + 'Potential ' + str(potfunc))
    i = 0
    dep_count = 0
    while i < steps - 1:

        if (method == "Infrequent WT MetaD"):
            (triggered,path) = potfunc.get_triggered(np.array([coords[i, 0],
                                                               coords[i, 1]]))
            if triggered is True:
                totaltime = time[i]
                teff = lf.calc_teff(walkerpot, beta, dt)
                return (totaltime, teff, info, path)

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

                    if xbc.type == 'Periodic':
                        history = np.vstack((history, np.array([coords[i, 0] +
                                             (xmax-xmin), coords[i, 1]])))
                        history = np.vstack((history, np.array([coords[i, 0] -
                                             (xmax-xmin), coords[i, 1]])))
                        w = np.append(w, winit)
                        w = np.append(w, winit)
                    if ybc.type == 'Periodic':
                        history = np.vstack((history, np.array([coords[i, 0],
                                                               coords[i, 1] +
                                                               (ymax-ymin)])))
                        history = np.vstack((history, np.array([coords[i, 0],
                                                               coords[i, 1] -
                                                               (ymax-ymin)])))
                        w = np.append(w, winit)
                        w = np.append(w, winit)

                elif (method == 'Well-Tempered Metadynamics' or
                        method == "Infrequent WT MetaD"):
                    VR = lf.calc_biased_pot(np.array([coords[i, 0],
                                                      coords[i, 1]]), history,
                                            w, delta, dimension)
                    history = np.vstack((history, np.array([coords[i, 0],
                                                           coords[i, 1]])))
                    w = np.append(w, winit * np.exp(-VR / (kb*DT)))

                    if xbc.type == 'Periodic':
                        history = np.vstack((history, np.array([coords[i, 0] +
                                             (xmax-xmin), coords[i, 1]])))
                        history = np.vstack((history, np.array([coords[i, 0] -
                                             (xmax-xmin), coords[i, 1]])))
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))
                    if ybc.type == 'Periodic':
                        history = np.vstack((history, np.array([coords[i, 0],
                                                               coords[i, 1] +
                                                               (ymax-ymin)])))
                        history = np.vstack((history, np.array([coords[i, 0],
                                                               coords[i, 1] -
                                                               (ymax-ymin)])))
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))
                        w = np.append(w, winit * np.exp(-VR / (kb*DT)))

        [pnew, vnew, newcoord, bcbias] = lf.integrate_step(coords[i], history,
                                                           w, delta, DT,
                                                           potfunc, bcs, p, m,
                                                           dt, gamma, beta,
                                                           dimension)
        p = pnew

        coords[i+1, 0] = newcoord[0]
        coords[i+1, 1] = newcoord[1]

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
        if makeplot == 'True' and sp.mod(i, plot_freq) == 0 and i > 0:
            # pdb.set_trace()

            for yc in range(0, ylong.size):
                for xc in range(0, xlong.size):
                    bias[yc, xc] = (bias[yc, xc] +
                                    lf.calc_biased_pot(np.array([xlong[xc],
                                                       ylong[yc]]),
                                                       history[dep_count:],
                                                       w[dep_count:],
                                                       delta, dimension))
            dep_count = len(history)
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
            cset3 = plt.contourf(xlong, ylong, bias+pot_base, levels,
                                 cmap=plt.cm.get_cmap(cmap, levels.size - 1))
            plt.colorbar(cset3)
            plt.scatter(coords[i+1, 0], coords[i+1, 1], marker='o',
                        color='r', zorder=10)
            plt.xlabel("CV1")
            plt.ylabel("CV2")
            plt.subplot(223)
            levels2 = np.arange(0, ebound[0]*-1, (ebound[0]*-1-0)/10)
            cset4 = plt.contourf(xlong, ylong, bias, levels2,
                                 cmap=plt.cm.get_cmap(cmap, levels2.size - 1))
            plt.colorbar(cset4)
            plt.scatter(coords[i+1, 0], coords[i+1, 1], marker='o',
                        color='r', zorder=10)
            plt.xlabel("CV1")
            plt.ylabel("CV2")
            # plt.draw()
            print i
            plt.pause(0.0001)
            if (make_movie == 'True'):
                filename = "movieframe" + str(frame)
                plt.savefig(filename + '.png', bbox_inches='tight')
                frame = frame + 1

        i = i + 1

    if(method != "Infrequent WT MetaD"):
        editFES = np.array(0)
        editvcalc = np.array(0)
        for yc in range(0, ylong.size):
            for xc in range(0, xlong.size):
                bias[yc, xc] = (bias[yc, xc] +
                                lf.calc_biased_pot(np.array([xlong[xc],
                                                   ylong[yc]]),
                                                   history[dep_count:],
                                                   w[dep_count:],
                                                   delta, dimension))
        colvar100 = lf.calc_colvar_2D(coords, bias,
                                      xlong, ylong, method, beta, T, DT)
        rmsds = lf.calc_rmsd(FES, beta, pot_base)
        pdb.set_trace()
        return (coords, E, rmsds, info)

    elif(method == "Infrequent WT MetaD"):
        teff = 0
        info = info + 'NO RARE EVENT'
        totaltime = 0
        path="NULL"
        return (coords, E, teff, totaltime, info, path)


def recreate_2DFES(FES, icount, coords, xinc, xmin, xmax, yinc, ymin, ymax, E):
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
