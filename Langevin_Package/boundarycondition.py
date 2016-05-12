"""
The boundary condition class defines the properties and functions of 1-D
boundary conditions to be implemented on a potential. For a 2-D potential
two separate 1-D boundary conditions are implemented separately (1 per axis).
"""
import numpy as np
import pdb


class BoundaryCondition:

    def __init__(self):
        """
        Initializes the boundary condition. The default location at which the
        boundary condition is implemented is negative and positive infinity.
        """
        self.type = "BoundaryCondition"
        self.location = np.array([float('inf')*-1, float('inf')])

    def get_potential(self, coords, origin_potential):
        """
        Returns an adjusted potential energy from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_potential :   float
                             Potential energy as if no boundary present

        Returns:
        -----------
        origin_potential :   float
                             Potential energy as if no boundary present.
                             (Default is no effect of boundary on potential)
        """
        return origin_potential

    def get_force(self, coords, origin_force):
        """
        Returns an adjusted force from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        origin_force     :   float
                             Force as if no boundary present.
                             (Default is no effect of boundary on force)
        """
        return origin_force

    def get_new_location(self, coords):
        """
        Returns an adjusted coordinates from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        coords           :   float
                             Coordinates as if no boundary present.
                             (Default is no effect of boundary on coordinates)
        """
        return coords

    def add_depositions(self, coords, history):
        """
        Returns an adjusted hills deposition from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        history           :  float
                             Return unchanges bias deposition locations
                             (Default is no effect of boundary on bias dep)
        """
        return history

    def get_bc_bias(self, coords):
        """
        Returns bias from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        bcbias           :  float
                            Return the difference between potential with and
                            without boundary (for bias calculation)
                             (Default is no effect of boundar)
        """
        bcbias = 0.0
        return bcbias

    def set_new_lower_location(self, lower_location):
        """
        Set location of boundary (low value)

        Parameters:
        -----------
        lower_location      :   float
                                Defines the coordinate of the walker

        """
        self.location[0] = lower_location

    def set_new_upper_location(self, upper_location):
        """
        Set location of boundary (low value)

        Parameters:
        -----------
        upper_location      :   float
                                Defines the coordinate of the walker

        """
        self.location[1] = upper_location


class QuarticBoundary(BoundaryCondition):

    def __init__(self):
        """
        Initializes the Quartic boundary condition. The default location at
        which the boundary condition is implemented is negative and
        positive infinity.
        """
        self.type = 'Quartic'
        self.location = np.array([float('inf')*-1, float('inf')])

    def get_potential(self, coords, origin_potential):
        """
        Returns an adjusted potential energy from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_potential :   float
                             Potential energy as if no boundary present

        Returns:
        -----------
        new_potential    :   float
                             Potential energy accounting for boundary

        """
        bcvalue = 0.0
        if coords > self.location[1]:
            bcvalue = 100.0*(coords - self.location[1])**4.0
        elif coords < self.location[0]:
            bcvalue = 100.0*(coords - self.location[0])**4.0
        new_potential = (bcvalue + origin_potential)
        return new_potential

    def get_force(self, coords,  origin_force):
        """
        Returns an adjusted force from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        origin_force     :   float
                             Force with impact from boundary.

        """
        if coords > self.location[1]:
            origin_force = -400.0*(coords - self.location[1])**3.0
        if coords < self.location[0]:
            origin_force = -400.0*(coords - self.location[0])**3.0
        return origin_force

    def get_bc_bias(self, coords):
        """
        Returns bias from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        bcbias           :  float
                            Return the difference between potential with and
                            without boundary (for bias calculation)
                             (Default is no effect of boundary for quartic)
        """
        bcbias = 0.0
        return bcbias


class PeriodicBoundary(BoundaryCondition):

    def __init__(self):
        """
        Initializes the Periodic boundary condition. The default location at
        which the boundary condition is implemented is negative and
        positive infinity.
        """
        self.type = 'Periodic'
        self.location = np.array([float('inf')*-1, float('inf')])

    def get_new_location(self, coords):
        """
        Returns an adjusted coordinates from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        newcoords         :  float
                             New coordinates from periodic boundary

        """
        movement = self.location[1]-self.location[0]
        if coords > self.location[1]:

            newcoords = coords - movement
        elif coords < self.location[0]:

            newcoords = coords + movement
        else:
            newcoords = coords
        return newcoords

    def add_depositions(self, coords, history):
        """
        Returns an adjusted hills deposition from the boundary

        Parameters:
        -----------
        coords           :   Array of floats
                             Defines the coordinate of the walker
        origin_force     :   float
                             Force as if no boundary present

        Returns:
        -----------
        history           :  float
                             Return additional hills to account for periodic
                             boundary (want hills to wrap around)
                             (Default is no effect of boundary on bias dep)
        """
        movement = self.location[1]-self.location[0]
        history = np.append(history, (coords + movement))
        history = np.append(history, (coords - movement))
        return history


def get_boundary_condition_dict():
    """Return a dictionary of all of the available potential functions."""
    bc_dict = {'Quartic': QuarticBoundary,
               'Periodic': PeriodicBoundary,
               'Empty': BoundaryCondition
               }

    return bc_dict
