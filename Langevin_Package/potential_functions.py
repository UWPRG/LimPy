"""This file contains the various potential functions we have defined"""

import numpy as np
import math
import pdb
import langevin_functions


def cosine_potential(coords):
    """
        Potential Energy and Force of Cos(x)

        Parameters:
        -----------

        coords  : float (or array of floats)
                Location


        Returns:
        --------

        V       : float (or array of floats)
                  Potential Energy

        F       : float (or array of floats)
                  Force

        Trigger : Boolean
                  Has rare event occurred (True) or not (False)
    """

    V = np.cos(coords) * 2.5
    F = np.sin(coords) * 2.5

    if type(coords) is np.float64 and (coords < 0 or coords > 6.3):
        Trigger = True
    else:
        Trigger = False

    return (V, F, Trigger)


def two_gaussian_potential(coords):
    """
        Calculate the force and potential based on location (1-D).

        Parameters:
        -----------

        coords  :  array of floats
                   X and Y coordinates


        Returns:
        --------

        V       : float (or array of floats)
                  Potential Energy

        F       : float (or array of floats)
                  Force

        Trigger : Boolean
                  Has rare event occurred (True) or not (False)
    """
    V = (-5 * np.exp(-(coords - 2/0.75)**2) -
         10*np.exp(-(coords + 2/0.75)**2))
    F = ((-5 * 2 * -1 * (coords - 2/0.75) * np.exp(-(coords - 2/0.75)**2) -
         10 * 2 * -1 * (coords + 2/0.75) * np.exp(-(coords + 2/0.75)**2)) *
         (-1))

    if type(coords) is np.float64 and coords < -2.0:
        Trigger = True
    else:
        Trigger = False

    return (V, F, Trigger)


def pv_2D_potential(x, y):
    """
        Calculate the force and potential based on location (2-D).
        Parameters:
        -----------

        x       : array of floats
                  X coordinates

        y       : array of floats
                  Y coordinates

        Returns:
        --------

        V       : float (or array of floats)
                  Potential Energy

        Fpot    : float (or array of floats)
                  Force in x and y direction

        Trigger : Boolean
                  Has rare event occurred (True) or not (False)
    """

    if type(x) is not np.float64:

        V = np.zeros([y.size, x.size])
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
        Trigger = False

    else:

        V = (np.cos(2*math.pi*x)*(1+4*y) + math.pi*y**2 -
             0.75*np.cos(2*math.pi*x/3))
        Fpotx = (((2*math.pi/3*0.75)*np.sin(2*math.pi*x/3) -
                 2*math.pi*(1+4*y)*np.sin(2*math.pi*x)))*-1
        Fpoty = ((2*math.pi*y+4*np.cos(2*math.pi*x)))*-1

        if x < 0.75 and y > 0 or x > 2.25 and y > 0:
            Trigger = True
        else:
            Trigger = False
    Fpot = np.array([Fpotx, Fpoty])

    return (V, Fpot, Trigger)


def get_potential_dict():
    """Returns a dictionary of all of the available potential functions
    """

    potential_dict = {'cosine_potential': cosine_potential,
                      'two_gaussian_potential': two_gaussian_potential,
                      'pv_2D_potential': pv_2D_potential}

    return potential_dict


def two_gaussian_potential_bc(vnew, f2, coords):
    """Applies Boundary Condition to the potential, force, and coordinates

        Parameters:
        -----------
        vnew       : float (or array of floats)
                     Potential Energy

        f2         : float (or array of floats)
                     Force

        coords     : float
                     coordinates

        Returns:
        --------

        vnew       : float (or array of floats)
                     Adjusted potential energy from boundary condition

        F          : float (or array of floats)
                     Adjusted force from boundary condition

        coords     : float
                     adjusted coordinates from boundary condition

        bcbias     : float
                     bias applied strictly from the boundary condition
    """
    vold = vnew
    bcbias = 0
    if (coords < -4.3193):

        vnew = 100.0 * (coords+4.0)**4.0 - 1.690133
        f2 = -100.0 * 4.0 * (coords+4.0)**3.0
        bcbias = vnew - vold
    elif (coords > 4.25882):

        vnew = 100.0 * (coords-4.0)**4.0 - 0.845067
        f2 = -100.0 * 4.0 * (coords-4.0)**3.0
        bcbias = vnew - vold
    return (vnew, f2, coords, bcbias)


def pv_2D_potential_bc(vnew, f2, coords):
    """Applies Boundary Condition to the potential, force, and coordinates

        Parameters:
        -----------
        vnew       : float (or array of floats)
                     Potential Energy

        f2         : float (or array of floats)
                     Force

        coords     : float
                     coordinates

        Returns:
        --------

        vnew       : float (or array of floats)
                     Adjusted potential energy from boundary condition

        F          : float (or array of floats)
                     Adjusted force from boundary condition

        coords     : float
                     adjusted coordinates from boundary condition

        bcbias     : float
                     bias applied strictly from the boundary condition
    """

    if (coords[0] < 0):
        coords[0] = coords[0] + 3

    elif (coords[0] > 3.0):
        coords[0] = coords[0] - 3
    bcbias = 0
    return (vnew, f2, coords, bcbias)


def cosine_potential_bc(vnew, f2, coords):
    """Applies Boundary Condition to the potential, force, and coordinates

        Parameters:
        -----------
        vnew       : float (or array of floats)
                     Potential Energy

        f2         : float (or array of floats)
                     Force

        coords     : float
                     coordinates

        Returns:
        --------

        vnew       : float (or array of floats)
                     Adjusted potential energy from boundary condition

        F          : float (or array of floats)
                     Adjusted force from boundary condition

        coords     : float
                     adjusted coordinates from boundary condition

        bcbias     : float
                     bias applied strictly from the boundary condition
    """

    if (coords < 0):
        coords = coords + 2*np.pi

    elif (coords > 2*np.pi):
        coords = coords - 2*np.pi
    bcbias = 0
    return (vnew, f2, coords, bcbias)


def get_boundary_condition_dict():
    """Returns a dictionary of all of the available potential functions
    """

    bc_dict = {'cosine_potential': cosine_potential_bc,
               'two_gaussian_potential': two_gaussian_potential_bc,
               'pv_2D_potential': pv_2D_potential_bc}

    return bc_dict
