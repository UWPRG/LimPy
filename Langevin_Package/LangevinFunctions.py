"""These are the functions for Langevin Integrator, Force, and Rare Events."""

# import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

import pdb
import math
import os
# import sys


def LIMD(inps, mdps, potdim, sm, movieflag):
    """Perform calcs for langevin integrator with enhanced sampling."""
    # If statements to receive inputs based on systems
    if(potdim == '1-D Potential'):
        steps = inps[0]
        dt = inps[1]
        x0 = inps[2]
        T = inps[3]
        m = inps[4]
        xmin = inps[5]
        xmax = inps[6]
        xinc = inps[7]
    if (potdim == '2-D Potential'):
        steps = inps[0]
        dt = inps[1]
        x0 = inps[2]
        T = inps[3]
        m = inps[4]  # units of kg
        xmin = inps[5]
        xmax = inps[6]
        xinc = inps[7]
        y0 = inps[8]
        ymin = inps[9]
        ymax = inps[10]
        yinc = inps[11]
    if (sm == 'Metadynamics'):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        w = np.array([0.0])
        DT = float("inf")

    if (sm == 'Well-Tempered Metadynamics'):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        DT = mdps[3]
        w = np.array([0.0])
    if (sm == 'MD'):
        winit = 0
        delta = 1
        hfreq = steps*2
        DT = 10000
        w = np.array([0.0])
    if (sm == "Infrequent WT MetaD"):
        winit = mdps[0]
        delta = mdps[1]
        hfreq = mdps[2]
        DT = mdps[3]
        w = np.array([0.0])

    # How often movie is created if desired
    iratio = 1000

    # print('Parameters:')
    # print('Number of steps '+ str(steps))
    # print('Initial x coordinate '+str(x0))
    # if (potdim == '2-D Potential'):
    #    print('Initial y coordinate '+str(y0))
    # print('Temperature '+str(T))
    # print('Timestep '+str(dt))

    # Parameters for integrator
    gamma = 0.5  # Friction factor
    beta = 1 / T / (1.987E-3)  # units of 1/kcal
    c1 = np.exp(-gamma * dt / 2)  # (Eq.13a)
    c2 = np.sqrt((1 - c1**2) * m / beta)  # (Eq.13b)
    # print ('c1: '+str(c1))
    # print ('c2: '+str(c2))
    # print 1/beta #kt =[kcal]
    # Determines if .png files created
    if (movieflag == 1):
        if os.path.exists("movies125"):
            os.rmdir("movies125")
        os.mkdir("movies125")
        frame = 0
        os.chdir("movies125")

    # Initial Configuration for 1-D
    if (potdim == '1-D Potential'):
        # pdb.set_trace()
        xlong = np.arange(xmin, xmax+xinc, xinc)
        q = np.empty(int(steps))
        E = np.empty(int(steps))
        s = np.array([0.0])

        time = np.array([0.0])
        walkerpot = np.array([0.0])
        vcalc = np.zeros_like(xlong)
        for xc in range(0, xlong.size):
            vcalc[xc] = force(xlong[xc], s, w, delta, DT, winit)[0]
        # Fbase = force(xlong, s, w, delta, DT, winit)
        # vcalc = Fbase[0]
        # fb = Fbase[1]
        q[0] = x0
        # Initial Velocity
        v0 = sp.rand(1)-0.5

        p = v0*m
        F1 = force(q[0], s, w, delta, DT, winit)
        v1 = F1[0]
        #  first = F1[1]

        E[0] = 0.5*p**2 + v1
        FES = np.zeros_like(xlong)

        icount = np.zeros_like(xlong)
        # Initial plotting
        # plt.ion()
        # plt.plot(xlong, vcalc, '-b')
        # plt.plot(x0, v1, 'ro', markersize=10)
        # plt.axis([-4, 4, -12, 6])
        # plt.xlabel("CV(s)")
        # plt.ylabel("F")
        # plt.draw()
        # plt.pause(0.0001)
        # pdb.set_trace()
    # Initial Configuration for 2-D
    if (potdim == '2-D Potential'):
        # pdb.set_trace()
        xlong = np.arange(xmin, xmax+xinc, xinc)
        ylong = np.arange(ymin, ymax+yinc, yinc)
        sx = np.array([0.0])
        sy = np.array([0.0])

        FTD = tdforce(xlong, ylong, sx, sy, w, delta, DT, winit)
        vcalc = FTD[0]

        FES = np.zeros_like(vcalc)
        icount = np.zeros_like(FES)
        q = np.empty([int(steps), 2])
        E = np.empty(int(steps))
        time = np.array([0.0])
        walkerpot = np.array([0.0])
        q[0, 0] = x0
        q[0, 1] = y0

        v0x = sp.rand(1) - 0.5
        v0y = sp.rand(1) - 0.5
        px = v0x * m
        py = v0y * m

        FTD1 = tdforce(q[0, 0], q[0, 1], sx, sy, w, delta, DT, winit)
        v1 = FTD1[0]
        # firstx = FTD1[1]
        # firsty = FTD1[2]
        E[0] = 0.5 * (px**2 + py**2) + v1
        # Initial Plotting
        # plt.clf()
        # plt.ion()
        # plt.subplot(221)
        # cmap = plt.cm.PRGn
        # levels = np.arange(-3, 2, 0.25)
        # cset1 = plt.contourf(xlong, ylong, vcalc, levels,
        #                     cmap=plt.cm.get_cmap(cmap, levels.size - 1))
        # plt.colorbar(cset1)
        # plt.title('2-D Potential')
        # plt.xlabel("X")
        # plt.ylabel("Y")
        # plt.scatter(x0, y0, marker='o', color='r', zorder=10)
        # plt.subplot(222)
        # cset2 = plt.contourf(xlong, ylong, vcalc,levels,cmap =
        # plt.cm.get_cmap(cmap, levels.size - 1))
        # plt.colorbar(cset2)

        # plt.draw()
        # plt.pause(0.0001)

    # Iteration across steps
    i = int(0)
    teff = 0.0
    if (potdim == '1-D Potential'):
        info = ('Parameters: \n' + 'Number of steps: ' + str(steps) + '\n' +
                'Initial x coordinate ' + str(x0) + '\n' + 'Temperature ' +
                str(T) + '\n' + 'Timestep ' + str(dt) + '\n' +
                ' Initial Hill Height ' + str(winit) + '\n' +
                'Hill Width ' + str(delta) + '\n' +
                'Deposition Frequency (steps)' + str(hfreq) + '\n' +
                'Well Temperature ' + str(DT))
    elif (potdim == '2-D Potential'):
        info = ('Parameters: \n'+'Number of steps: ' + str(steps) + '\n' +
                'Initial x coordinate ' + str(x0) + '\n' +
                'Initial y coordinate ' + str(y0) + '\n' + 'Temperature ' +
                str(T) + '\n' + 'Timestep ' + str(dt) + '\n' +
                'Initial Hill Height ' + str(winit) + '\n' + 'Hill Width ' +
                str(delta) + '\n' + 'Deposition Frequency (steps)' +
                str(hfreq) + '\n' + 'Well Temperature ' + str(DT))
    # Loop over until rareevent happens or simulation time runs out.
    while i < steps - 1:

            if (potdim == '1-D Potential'):
                triggered = rareevent(q[i])

                # Checks if Rare Event has occured to stop
                if triggered is True and sm == "Infrequent WT MetaD":
                    totaltime = time[i]
                    walkercalc = np.delete(walkerpot, 0)
                    teff = sum(dt * np.exp(walkercalc * beta))
                    hillnum = i/hfreq
                    return (totaltime, teff, info, hillnum)
        # Check if deposit of gaussian
                if sp.mod(i, hfreq) == 0 and i > 0:

                    if(i == hfreq):
                        s[0] = q[i]
                        w[0] = winit
                    else:
                        if sm == 'Metadynamics':
                            s = np.append(s, q[i])
                            w = np.append(w, winit)
                        if (sm == 'Well-Tempered Metadynamics' or
                                sm == "Infrequent WT MetaD"):

                                vr = sum(w * np.exp(-(q[i] - s)**2 /
                                                    2 / delta**2))
                                s = np.append(s, q[i])
                                w = np.append(w, winit * np.exp(-vr /
                                                                (1.987E-3*DT)))
        # Calculate The force and potential at current step
                # pdb.set_trace()
                Fnow = force(q[i], s, w, delta, DT, winit)
                # v = Fnow[0]
                f = Fnow[1]

        # Langevin integrator equations from paper
                R1 = sp.rand(1) - 0.5
                R2 = sp.rand(1) - 0.5
                pplus = c1*p + c2*R1
                q[i+1] = q[i] + (pplus/m) * dt + f/m * ((dt**2) / 2)

                Ffut = force(q[i+1], s, w, delta, DT, winit)
                v2 = Ffut[0]
                f2 = Ffut[1]

                pminus = pplus + (f/2 + f2/2)*dt
                p = c1*pminus + c2*R2
                E[i+1] = 0.5 * p**2 + v2
                # time and walker potential track for infrequent sampling
                time = np.append(time, dt*(i+1))
                if q[i+1] > 4:
                    walkerpot = np.append(walkerpot,
                                          sum(w * np.exp(-(q[i+1]-s)**2 /
                                              2 / delta**2)) +
                                          ((100 * (q[i+1] - 4)**4)))

                else:
                    walkerpot = np.append(walkerpot,
                                          sum(w * np.exp(-(q[i+1]-s)**2 /
                                              2 / delta**2)))

                # indexing for saving energy values in grid

                index = int(round((round(q[i+1], int(abs(math.log10(xinc)))) +
                            (0-xmin))/xinc))
                # print index
                # print q[i+1]
                if q[i + 1] > xmin and q[i + 1] < xmax:
                    FES[index] = ((FES[index] * (icount[index]) + E[i+1]) /
                                  (icount[index] + 1))
                    icount[index] = icount[index] + 1

                # updating plot
                if sp.mod(i, iratio) == 0 and i > 0:
                    # print('Step '+ str(i))
                    # print('X coordinate '+str(q[i]))
                    # print('Force '+str(f))
                    # print('Energy '+str(E[i+1]))
                    bias = np.copy(vcalc)

                    if s.size > 1:
                        k = 0
                        while k < xlong.size:
                            bias[k] = (vcalc[k] +
                                       sum(w * np.exp(-(xlong[k] - s)**2 /
                                                      2 / delta**2)))
                            k = k + 1
            # # pdb.set_trace()
            # # Walker

                        v2 = v2 + sum(w * np.exp(-(q[i+1] - s)**2 /
                                                 2/delta**2))

            # Plotting

                    # plt.clf()
                    # plt.plot(xlong, bias, '-r')
                    # plt.plot(xlong, vcalc, '-b')
                    # plt.plot(q[i+1], v2, 'ro', markersize=10)
                    # plt.axis([xmin, xmax, min(vcalc)-8, max(vcalc)+8])
                    # plt.xlabel("CV(s)")
                    # plt.ylabel("F")
                    # plt.draw()
                    # plt.pause(0.0001)
                    if (movieflag == 1):
                        filename = "movieframe" + str(frame)
                        plt.savefig(filename + '.png', bbox_inches='tight')
                        frame = frame + 1

    # Two Dimensional Potential
            if (potdim == '2-D Potential'):
                triggered = twodrarevent(q[i, 0], q[i, 1])
                if triggered is True and sm == "Infrequent WT MetaD":
                    totaltime = time[i]
                    print q[i, 0]
                    print q[i, 1]
                    walkercalc = np.delete(walkerpot, 0)
                    teff = sum(dt * np.exp(walkercalc * beta))
                    return (totaltime, teff, info)
                if sp.mod(i, hfreq) == 0 and i > 0:

                    if(i == hfreq):
                        # pdb.set_trace()
                        sx[0] = q[i, 0]
                        sy[0] = q[i, 1]
                        w[0] = winit
                    else:
                        # pdb.set_trace()
                        VR = sum(w*np.exp(-(q[i, 0] - sx)**2 / 2 / delta**2) *
                                 np.exp(-(q[i, 1] - sy)**2 / 2 / delta**2))
                        w = np.append(w, winit * np.exp(-VR / (1.987E-3*DT)))
                        sx = np.append(sx, q[i, 0])
                        sy = np.append(sy, q[i, 1])
        # Calculate The force and potential at current step
                # pdb.set_trace()
                Fnow = tdforce(q[i, 0], q[i, 1], sx, sy, w, delta, DT, winit)
                # v = Fnow[0]
                fx = Fnow[1]
                fy = Fnow[2]

        # integrator stuff from paper
                R1x = sp.rand(1) - 0.5
                R2x = sp.rand(1)-0.5
                pplusx = c1*px + c2*R1x
                q[i+1, 0] = q[i, 0] + (pplusx/m)*dt + fx/m*((dt**2)/2)

                R1y = sp.rand(1) - 0.5
                R2y = sp.rand(1) - 0.5
                pplusy = c1*py + c2*R1y
                q[i+1, 1] = q[i, 1] + (pplusy/m)*dt + fy/m*((dt**2)/2)

                Ffut = tdforce(q[i+1, 0], q[i+1, 1], sx, sy, w, delta, DT,
                               winit)
                v2 = Ffut[0]
                f2x = Ffut[1]
                f2y = Ffut[2]

                pminusx = pplusx + (fx/2 + f2x/2)*dt
                pminusy = pplusy+(fy/2 + f2y/2)*dt
                px = c1*pminusx + c2*R2x
                py = c1*pminusy + c2*R2y
                E[i+1] = 0.5*(px**2 + py**2) + v2
                time = np.append(time, dt*(i + 1))
                walkerpot = np.append(walkerpot,
                                      sum(w * np.exp(-((q[i+1, 0] - sx)**2) /
                                          2 / delta**2) *
                                          np.exp(-((q[i+1, 1] - sy)**2) / 2 /
                                                 delta**2)))

                xindex = int(round((round(q[i+1, 0],
                                    int(abs(math.log10(xinc)))) +
                                    (0 - xmin)) / xinc))

                yindex = int(round((round(q[i+1, 1],
                                    int(abs(math.log10(yinc)))) +
                                    (0 - ymin)) / yinc))
                if (q[i+1, 0] > xmin and q[i+1, 0] < xmax and
                        q[i+1, 1] > ymin and q[i+1, 1] < ymax):

                    FES[yindex, xindex] = ((FES[yindex, xindex] *
                                           (icount[yindex, xindex]) + E[i+1]) /
                                           (icount[yindex, xindex] + 1))
                    icount[yindex, xindex] = icount[yindex, xindex] + 1

                # if sp.mod(i, iratio) == 0 and i > 0:
                    # pdb.set_trace()
                    # bias = np.copy(vcalc)

                    # if sx.size > 1:
                    #     # pdb.set_trace()
                    #     xcoord = np.arange(-1, 1.01, .01)
                    #     ycoord = np.arange(-2, 2.01, .01)
                    #     for k in range(0, ycoord.size-1):
                    #         for j in range(0, xcoord.size-1):
                    #             bias[k, j] = (bias[k, j] +
                    #                           sum(w * np.exp(-((xcoord[j] -
                    #                                             sx)**2) / 2 /
                    #                                          delta**2) *
                    #                             np.exp(-((ycoord[k]-sy)**2) /
                    #                               2 / delta**2)))
                    # v = v2 + sum(w * np.exp(-((q[i+1, 0]-sx)**2) / 2 /
                    #             delta**2)*np.exp(-((q[i+1, 1]-sy)**2) /
                    #                              2 / delta**2))
                    # Plotting
                    # plt.clf()
                    # cmap = plt.cm.PRGn
                    # levels = np.arange(-5, 20, 0.1)

                    # plt.subplot(221)
                    # plt.title('2-D Potential')
                    # plt.xlabel("X")
                    # plt.ylabel("Y")
                    # xrare = np.array([-1, 1])
                    # xyrare = np.array([-0.5, -0.5])
                    # yrare = np.array([-2, 2])
                    # yxrare = np.array([0.1, 0.1])
                    # xpoints = np.array([q[0, 0]])
                    # ypoints = np.array([q[0, 1]])
                    # for b in range(0, i):
                    #     xpoints = np.append(xpoints, q[b, 0])
                    #     ypoints = np.append(ypoints, q[b, 1])

                    # plt.plot(xpoints, ypoints, '--r', linewidth=1.0)
                    # cset1 = plt.contourf(xlong, ylong, vcalc, levels,
                    #                       cmap=plt.cm.get_cmap(cmap,
                    #                                         levels.size - 1))
                    # plt.colorbar(cset1)
                    # plt.scatter(q[i+1, 0], q[i+1, 1], marker='o',
                    #              color='r', zorder=10)
                    # plt.xlim(0, 3)
                    # plt.ylim(-2, 2)
                    # plt.plot(xrare, xyrare)
                    # plt.plot(yxrare, yrare)
                    # plt.subplot(222)
                    # cset2 = plt.contourf(xlong, ylong, vcalc, levels,
                    #                      cmap=plt.cm.get_cmap(cmap,
                    #                                         levels.size - 1))
                    # plt.colorbar(cset2)
                    # plt.scatter(q[i+1, 0], q[i+1, 1], marker='o',
                    #             color='r', zorder=10)
                    # plt.subplot(223)
                    # cset3 = plt.contourf(xlong, ylong, bias, levels,
                    #                      cmap=plt.cm.get_cmap(cmap,
                    #                                         levels.size - 1))
                    # plt.colorbar(cset3)
                    # plt.scatter(q[i+1, 0], q[i+1, 1], marker='o',
                    #             color='r', zorder=10)
                    # plt.subplot(224)
                    # cset4 = plt.contourf(xlong, ylong, bias-vcalc, levels,
                    #                      cmap=plt.cm.get_cmap(cmap,
                    #                                         levels.size - 1))
                    # plt.colorbar(cset4)
                    # plt.scatter(q[i+1, 0], q[i+1, 1], marker='o',
                    #             color='r', zorder=10)

                    # plt.draw()
                    # plt.pause(0.0001)

                    if (movieflag == 1):
                        filename = "movieframe" + str(frame)
                        plt.savefig(filename + '.png', bbox_inches='tight')
                        frame = frame + 1
                    print('Step ' + str(i))
                    print('X coordinate ' + str(q[i+1, 0]))
                    print('Y coordinate ' + str(q[i+1, 1]))
                    print('Energy ' + str(E[i+1]))
            i = i + 1
    # if non infrequent sampling, calculate RMSD for convergence
    if(potdim == '1-D Potential' and sm != "Infrequent WT MetaD"):
        rmsd = np.sqrt(np.sum((((FES-vcalc)) * beta)**2) / FES.shape[0])
        rmskld = np.sqrt(np.sum(np.multiply(np.power(((FES - FES.mean()) -
                         (vcalc-vcalc.mean())) * beta, 2),
                         np.exp(-vcalc*beta))) / np.sum(np.exp(-vcalc * beta)))
        rmsalignerr = np.sqrt(np.sum(np.power(((FES-FES.mean()) -
                              (vcalc-vcalc.mean()))*beta, 2)) / np.size(vcalc))
    elif(potdim == '2-D Potential'and sm != "Infrequent WT MetaD"):

        editFES = np.array(0)
        editvcalc = np.array(0)
        for k in range(0, ylong.size):
            for j in range(0, xlong.size):
                if(icount[k, j] > 0):
                    editFES = np.append(editFES, FES[k, j])
                    editvcalc = np.append(editvcalc, vcalc[k, j])
        editFES = np.delete(editFES, 0)
        editvcalc = np.delete(editvcalc, 0)

        rmsd = np.sqrt(np.sum((((editFES-editvcalc)) * beta)**2) /
                       np.size(editFES))

        rmskld = np.sqrt(np.sum(np.multiply(np.power(((editFES -
                                editFES.mean()) -
                         (editvcalc-editvcalc.mean()))*beta, 2),
                         np.exp(-editvcalc*beta))) /
                         np.sum(np.exp(-editvcalc*beta)))

        rmsalignerr = np.sqrt(np.sum(np.power(((editFES-editFES.mean()) -
                              (editvcalc-editvcalc.mean()))*beta, 2)) /
                              np.size(editFES))
    # export info if infrequent sampling but rare event not triggered
    if(sm == "Infrequent WT MetaD"):
        totaltime = time[i]
        walkercalc = np.delete(walkerpot, 0)
        teff = sum(dt*np.exp(walkercalc*beta))
        return (totaltime, teff, info)
    return(rmsd, rmskld, rmsalignerr, info)

# Function for the underlying potential and force calculations


def force(r, s, w, delta, DT, winit):
    """Calculate the force and potential based on location (1-D)."""
    if (r < -4):
        V = 100 * (r+4)**4
        Fpot = 100 * 4 * (r+4)**3

    elif (r > 4):
        V = 100 * (r-4)**4
        Fpot = 100 * 4 * (r-4)**3

    else:
        V = (-5 * np.exp(-(r - 2/0.75)**2) - 10*np.exp(-(r + 2/0.75)**2))
        Fpot = (-5 * 2 * -1 * (r - 2/0.75) * np.exp(-(r - 2/0.75)**2) - 10 *
                2 * -1 * (r + 2/0.75) * np.exp(-(r + 2/0.75)**2))
    Fbias = sum(w * (r-s) / delta**2 * np.exp(-(r-s)**2 / 2 / delta**2))
    F = Fpot * -1 + Fbias
    return np.array([V, F])


def tdforce(x, y, sx, sy, w, delta, DT, winit):
    """Calculate the force and potential based on location (2-D)."""
    if type(x) is not np.float64:
        V = np.empty([y.size, x.size])
        Fx = np.empty([y.size, x.size])
        Fy = np.empty([y.size, x.size])
        for k in range(0, y.size - 1):
            for j in range(0, x.size - 1):
                V[k, j] = (np.cos(2*math.pi*x[j])*(1 + 4*y[k]) +
                           math.pi*y[k]**2 - 0.75*np.cos(2*math.pi*x[j]/3))
                Fx = (((2*math.pi/3*0.75)*np.sin(2*math.pi*x[j]/3) -
                      2*math.pi*(1+4*y[k])*np.sin(2*math.pi*x[j])))
                Fy = ((2*math.pi*y[k]+4*np.cos(2*math.pi*x[j])))
            Fpotx = Fx * -1
            Fpoty = Fy * -1

    elif type(x) is np.float64:

        V = (np.cos(2*math.pi*x)*(1+4*y) + math.pi*y**2 -
             0.75*np.cos(2*math.pi*x/3))
        Fx = (((2*math.pi/3*0.75)*np.sin(2*math.pi*x/3) -
              2*math.pi*(1+4*y)*np.sin(2*math.pi*x)))
        Fy = ((2*math.pi*y+4*np.cos(2*math.pi*x)))

        Fbiasx = sum(w*((x-sx)/delta**2) * np.exp(-(x-sx)**2/2/delta**2) *
                     np.exp(-(y-sy)**2/2/delta**2))
        Fbiasy = sum(w*((y-sy)/delta**2) * np.exp(-(x-sx)**2/2/delta**2) *
                     np.exp(-(y-sy)**2/2/delta**2))
        Fpotx = Fx*-1+Fbiasx
        Fpoty = Fy*-1+Fbiasy

        # Boundaries to keep walker in system
    return [V, Fpotx, Fpoty]

# Functions defining rare events


def rareevent(X):
    """Define Rare event for (1-D) Infrequent MetaD."""
    if X < -1:
        return True
    else:
        return False


def twodrarevent(X, Y):
    """Define Rare event for (2-D) Infrequent MetaD."""
    if X < 0.75 and Y > 0 or X > 2.25 and Y > 0:
        return True
    else:
        return False
