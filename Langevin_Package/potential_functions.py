"""This file contains the various potential functions we have defined"""


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

    V = np.cos(coords)
    Fpot = (-1) * np.sin(coords)

    if coords < -2.0:
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

    """
    if (type(coords) is np.float64 and coords < -4):
        V = 100 * (coords+4)**4
        F = -100 * 4 * (coords+4)**3
    elif (type(coords) is np.float64 and coords > 4):
        V = 100 * (coords-4)**4
        F = -100 * 4 * (coords-4)**3
    else:
        V = (-5 * np.exp(-(coords - 2/0.75)**2) -
             10*np.exp(-(coords + 2/0.75)**2))
        F = (-5 * 2 * -1 * (coords - 2/0.75) * np.exp(-(coords - 2/0.75)**2) -
             10 * 2 * -1 * (coords + 2/0.75) * np.exp(-(coords + 2/0.75)**2))
        # len2 = s.size-1
        # VR = sum(w*np.exp(-(r-s)**2 / 2 / delta**2))
        # w[len2] = winit * np.exp(-VR / (1.987E-3*DT))
        # Fbias = sum(w * (r-s) / delta**2 * np.exp(-(r-s)**2 / 2 / delta**2))
        # F = Fpot * -1 + Fbias
    return (V, F, Trigger)


def pratyush_voter_2D_potential(coords):
    """
        Calculate the force and potential based on location (2-D).
        Parameters:
        -----------

        coords  :  array of floats
                   X and Y coordinates

        Returns:
        --------

        V       : float (or array of floats)
                  Potential Energy

        Fpotx   : float (or array of floats)
                  Force in the x direction

        Fpoty   : float (or array of floats)
                  Force in the y direction

    """
    x = coords[0]
    y = coords[1]
    sx = history[0]
    sy = history[1]
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
        Fpotx = (((2*math.pi/3*0.75)*np.sin(2*math.pi*x/3) -
                 2*math.pi*(1+4*y)*np.sin(2*math.pi*x)))*-1
        Fpoty = ((2*math.pi*y+4*np.cos(2*math.pi*x)))*-1

        # len2 = sx.size - 1
        # VR = sum(w*np.exp(-(x-sx)**2/2/delta**2) *
        #          np.exp(-(y-sy)**2/2/delta**2))
        # w[len2] = winit*np.exp(-VR/(1.987E-3*DT))
        #
        # Fbiasx = sum(w*((x-sx)/delta**2) * np.exp(-(x-sx)**2/2/delta**2) *
        #              np.exp(-(y-sy)**2/2/delta**2))
        # Fbiasy = sum(w*((y-sy)/delta**2) * np.exp(-(x-sx)**2/2/delta**2) *
        #              np.exp(-(y-sy)**2/2/delta**2))
        # Fpotx = Fx*-1+Fbiasx
        # Fpoty = Fy*-1+Fbiasy

    return (V, Fpotx, Fpoty)


def get_potential_dict():
    """Returns a dictionary of all of the available potential functions
    """

    potential_dict = {'cosine_potential': cosine_potential,
                      'two_gaussian_potential': two_gaussian_potential,
                      'pratyush_voter_2D_potential':
                      pratyush_voter_2D_potential}

    return potential_dict
