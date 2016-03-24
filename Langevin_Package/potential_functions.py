"""This file contains the various potential functions we have defined."""

import numpy as np
import math


def cosine_potential(coords):
    """
    Calculate the Potential Energy and Force of Cos(x).

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
    if hasattr(x, "__len__") is True:

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

#def get_distances(x,y,z):
def get_distances(coord):
  Natoms = len(coord)
  distance = []
  for atom_i in range(Natoms):
    matrix_line = []
    x_i, y_i, z_i = coord[atom_i]
    for atom_j in range(Natoms):
      x_j, y_j, z_j = coord[atom_j]
      d = sqrt( (x_i-x_j)**2 + (y_i-y_j)**2 + (z_i-z_j)**2 )
      matrix_line.append( d )
    distance.append(matrix_line)
  return distance
#def read_data	
def read_atom_coord(Natoms):
  print ("For each atom: x, y, z")
  result = []
  for i in range(Natoms):
    line = sys.stdin.readline()
    xyz = [ float(i) for i in line.split() ]
    result.append(xyz)
  return result

def LJ_potential(coords):
    """
    Calculate the Potential Energy and Force from a Lennard-Jones Potential.

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
    sigma=1.#0.3345#1.
    epsilon=10.#0.0661#1.
    
    Natoms=2.
    #distance = get_distances(coord)
    V = 0.
    F = 0.
    r=coords#coord
    #for atom_i in range(Natoms):
    #   for atom_j in range(atom_i+1,Natoms):
    V = 4.*epsilon*((sigma/r)**12 - (sigma/r)**6)
   # F = (48/r**2)*((sigma/r**12)**12-(0.5*(sigma/r**6))) # should be like that for 3D V=V(x,y,z)
    F = -(24*(r**6-2.)*epsilon)/r**13
          #update force??

    if type(coords) is np.float64 and coords < -2.0:
        Trigger = True
    else:
        Trigger = False
    return (V, F, Trigger)
 

def muller_brown_potential(x, y):
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
        A = np.array([-200.0, -100.0, -170.0, 15.0])
        a = np.array([-1.0, -1.0, -6.50, 0.7])
        b = np.array([0.0, 0.0, 11.0, 0.6])
        c = np.array([-10.0, -10.0, -6.50, 0.7])
        x0 = np.array([1.0, 0.0, -0.50, -1.0])
        y0 = np.array([0.0, 0.5, 1.50, 1.0])
        for k in range(0, y.size - 1):
            for j in range(0, x.size - 1):
                V[k, j] = sum(A * np.exp(a*(x[j]-x0)**2 +
                                         b*(x[j]-x0)*(y[k]-y0) +
                                         c*(y[k]-y0)**2))
                Fx = (-400*np.exp(-1*(-1+x[j])**2 - 10*y[k]**2)*(-1+x[j]) -
                      200*np.exp(-x[j]**2-10*(-0.5+x[j])**2)*x[j] +
                      170*np.exp(-6.5*(0.5+x[j])**2+11*(0.5+x[j])*(-1.5+y[k]) -
                      6.5*(-1.5+y[k])**2)*(-13*(0.5+x[j])+11*(-1/5+y[k])) -
                      15*np.exp(0.7*(1+x[j])**2+0.6*(1+x[j])*(y[k]-1) +
                      0.7*(y[k]-1)**2)*(1.4*(1+x[j])+0.6*(y[k]-1)))
                Fy = (170*np.exp(-6.5*(0.5+x[j])**2 +
                      11*(0.5+x[j])*(-1.5+y[k]) -
                      6.5*(y[k]-1.5)**2)*(11*(0.5+x[j])-13*(y[k]-1.5)) -
                      15*np.exp(0.7*(1+x[j])**2+0.6*(1+x[j])*(y[k]-1) +
                      0.7*(y[k]-1)**2)*(0.6*(x[j]+1)+1.4*(y[k]-1)) -
                      1000*np.exp(-x[j]**2-10*(y[k]-0.5)**2)*(y[k]-0.5) -
                      4000*np.exp(-1*(x[j]-1)**2-10*y[k]**2)*y[k])
        Fpotx = Fx * -1
        Fpoty = Fy * -1
        Trigger = False

    else:

        V = sum(A * np.exp(a * (x-x0)**2 + b * (x-x0)*(y-y0) +
                           c * (y-y0)**2))
        Fpotx = (-400*np.exp(-1*(-1+x)**2 - 10*y**2)*(-1+x) -
                 200*np.exp(-x**2-10*(-0.5+x)**2)*x +
                 170*np.exp(-6.5*(0.5+x)**2+11*(0.5+x)*(-1.5+y) -
                 6.5*(-1.5+y)**2)*(-13*(0.5+x)+11*(-1/5+y)) -
                 15*np.exp(0.7*(1+x)**2+0.6*(1+x)*(y-1) +
                 0.7*(y-1)**2)*(1.4*(1+x)+0.6*(y-1)))
        Fpoty = (170*np.exp(-6.5*(0.5+x)**2 +
                 11*(0.5+x)*(-1.5+y) -
                 6.5*(y-1.5)**2)*(11*(0.5+x)-13*(y-1.5)) -
                 15*np.exp(0.7*(1+x)**2+0.6*(1+x)*(y-1) +
                 0.7*(y-1)**2)*(0.6*(x+1)+1.4*(y-1)) -
                 1000*np.exp(-x**2-10*(y-0.5)**2)*(y-0.5) -
                 4000*np.exp(-1*(x-1)**2-10*y**2)*y)

        if x > 0.4 and y < 0.1:
            Trigger = True
        else:
            Trigger = False
    Fpot = np.array([Fpotx, Fpoty])

    return (V, Fpot, Trigger)


def get_potential_dict():
    """Return a dictionary of all of the available potential functions."""
    potential_dict = {'cosine_potential': cosine_potential,
                      'two_gaussian_potential': two_gaussian_potential,
                      'pv_2D_potential': pv_2D_potential,
                      'muller_brown_potential': muller_brown_potential,
		      'LJ_potential': LJ_potential}			

    return potential_dict


def two_gaussian_potential_bc(vnew, f2, coords):
    """
    Apply Boundary Condition to the potential, force, and coordinates.

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
    """
    Apply Boundary Condition to the potential, force, and coordinates.

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

def LJ_potential_bc(vnew,f2,coords):
    """
    Apply Boundary Condition to the potential, force, and coordinates.

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
    bcbias = 0
    return (vnew, f2, coords, bcbias)


def mb_2D_potential_bc(vnew, f2, coords):
    """
    Apply Boundary Condition to the potential, force, and coordinates.

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
    bcbias = 0
    return (vnew, f2, coords, bcbias)


def cosine_potential_bc(vnew, f2, coords):
    """
    Apply Boundary Condition to the potential, force, and coordinates.

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
    if (coords < -0.1*np.pi):
        coords = coords + 2.2*np.pi

    elif (coords > 2.1*np.pi):
        coords = coords - 2.2*np.pi
    bcbias = 0
    return (vnew, f2, coords, bcbias)


def get_boundary_condition_dict():
    """Return a dictionary of all of the available potential functions."""
    bc_dict = {'cosine_potential': cosine_potential_bc,
               'two_gaussian_potential': two_gaussian_potential_bc,
               'pv_2D_potential': pv_2D_potential_bc,
               'muller_brown_potential': mb_2D_potential_bc,
               'LJ_potential': LJ_potential_bc}

    return bc_dict
