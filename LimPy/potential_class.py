"""
This file defines the classes for 1 and 2 dimensional potential energy
surfaces, along with the standard properties and functions to normalize how
they can be used.
"""
import numpy as np
# for importing data to make a PES.
try:
    import cPickle as pickle
except:
    import pickle
import pdb


class PotentialFunction1D:
    """Basic Outline of a 1-D Potential Function"""
    dimension = '1-D Potential'

    def __init__(self):
        """
        Initializes the potential energy surface. Parameters are set to zero
        so there are no adjustments to default barriers and the rare events
        are to set negative and positive infinity so infrequent metadynamics
        is not used.
        """
        self.parameters = np.array([0, 0])
        self.rare_event = np.array([float("inf")*-1, float("inf")])

    def get_potential(self, coords):
        """
        Equation that returns the potential energy at a CV
        To be specified by child class
        """
        pass

    def set_parameters(self, values):
        """
        Overwrites default parameters for potential function.
        This is handy for scaling energy barriers.

        Parameters:
        -----------
        values    :     float
                        Parameters to scale/adjust potential energy values
        """
        self.parameters = values

    def get_force(self, coords):
        """
        Equation that returns the force at a CV
        To be specified by child class
        """
        pass

    def set_lower_rare_event(self, lower):
        """
        Overwrites default lower rare event for potential function
        Parameters:
        -----------
        lower    :      float
                        Set the lower bound for a rare event
        """
        self.rare_event[0] = lower

    def set_upper_rare_event(self, upper):
        """
        Overwrites default upper rare event for potential function
        Parameters:
        -----------
        upper    :      float
                        Set the upper bound for a rare event
        """
        self.rare_event[1] = upper

    def get_triggered(self, coords):
        """
        Determines if rare event has happened, and specifies which one if
        there are multiple
        Parameters:
        -----------
        coords    :      float
                         coordinates of walker.

        Returns:
        ----------
        trigger   :      boolean
                         True if rare event has happended, otherwise False

        event     :      String
                         Distinguishes which path/rare event occurs if multiple
                         paths are present
        """
        trigger = False
        event = "null"
        if coords < self.rare_event[0]:
            trigger = True
            event = "A"
        elif coords > self.rare_event[1]:
            trigger = True
            event = "B"
        return trigger, event


class PotentialFunction2D:
    """
    Basic Outline of a 2-D Potential Function
    """
    dimension = '2-D Potential'

    def __init__(self):
        """
        Initializes the potential energy surface. Parameters are set to zero
        so there are no adjustments to default barriers and the rare events
        are to set negative and positive infinity so infrequent metadynamics
        is not used.
        """
        self.parameters = np.array([0, 0, 0, 0])
        self.rare_event = np.array([[float("inf")*-1, float("inf")],
                                    [float("inf")*-1, float("inf")]])

    def get_potential(self, coords):
        """
        Equation that returns the potential energy at a CV
        To be specified by child class
        """
        pass

    def set_parameters(self, values):
        """
        Overwrites default parameters for potential function.
        This is handy for scaling energy barriers.

        Parameters:
        -----------
        values    :     array of floats
                        Parameters to scale/adjust potential energy values
        """
        self.parameters = values

    def get_force(self, coords):
        """
        Equation that returns the force at a CV
        To be specified by child class
        """
        pass

    def set_lower_rare_event(self, lowerx, lowery):
        """
        Overwrites default lower rare event for potential function
        Parameters:
        -----------
        lowerx    :      float
                        Set the lower bound for a rare event x coordinate

        lowery   :      float
                        Set the lower bound for a rare event y coordinate
        """
        self.rare_event[0][0] = lowerx
        self.rare_event[0][1] = lowery

    def set_upper_rare_event(self, upperx, uppery):
        """
        Overwrites default upper rare event for potential function
        Parameters:
        -----------
        upperx    :      float
                        Set the upper bound for a rare event x coordinate

        uppery   :      float
                        Set the upper bound for a rare event y coordinate
        """
        self.rare_event[1][0] = upperx
        self.rare_event[1][1] = uppery

    def get_triggered(self, coords):
        """
        Determines if rare event has happened, and specifies which one if
        there are multiple
        """
        pass


class CosinePotential(PotentialFunction1D):

    def __init__(self):
        """
        Initializes the 1-D potential energy surface for a cosine function.
        Parameters are set to zero so there are no adjustments to
        default barriers and rare events are initialized as well.
        """
        self.parameters = np.array([1, 0])
        self.rare_event = np.array([-2.0, 8.28])

    def get_potential(self, coords):
        """
        Returns the potential energy at a coordinate.

        Parameters:
        -----------
        coords      :   float
                        coordinate of walker

        Returns:
        ----------
        the_pot     :   float
                        potential energy at the present location.
        """
        the_pot = self.parameters[0]*np.cos(coords)+self.parameters[1]
        return the_pot

    def get_force(self, coords):
        """
        Returns the force acting on the walker at a coordinate.

        Parameters:
        -----------
        coords      :   float
                        coordinate of walker

        Returns:
        ----------
        the_force   :   float
                        forces at the present location.
        """
        the_force = self.parameters[0]*np.sin(coords)
        return the_force


class PieceWiseCosinePotential(PotentialFunction1D):

    def __init__(self):
        """
        Initializes the 1-D potential energy surface for a cosine function.
        Parameters are set to zero so there are no adjustments to
        default barriers and rare events are initialized as well.
        """
        self.parameters = np.array([1, 0, 1, 0])
        self.rare_event = np.array([-2.0, 8.28])

    def get_potential(self, coords):
        """
        Returns the potential energy at a coordinate.

        Parameters:
        -----------
        coords      :   float
                        coordinate of walker

        Returns:
        ----------
        V           :   float
                        potential energy at the present location.
        """
        if hasattr(coords, "__len__") is True:
            V = np.zeros_like(coords)
            for i in range(0, len(coords)):
                if coords[i] <= np.pi:
                    V[i] = (self.parameters[0]*np.cos(coords[i]) +
                            self.parameters[1])

                elif coords[i] > np.pi:
                    V[i] = (self.parameters[2]*np.cos(coords[i]) +
                            self.parameters[3])
            return V
        else:
            if coords <= np.pi:
                V = self.parameters[0]*np.cos(coords)+self.parameters[1]
                return V

            elif coords > np.pi:
                V = self.parameters[2]*np.cos(coords)+self.parameters[3]
                return V

    def get_force(self, coords):
        """
        Returns the force acting on the walker at a coordinate.

        Parameters:
        -----------
        coords      :   float
                        coordinate of walker

        Returns:
        ----------
        the_force   :   float
                        forces at the present location.
        """
        if coords <= np.pi:
            V = self.parameters[0]*np.sin(coords)
            return the_force

        elif coords > np.pi:
            V = self.parameters[2]*np.sin(coords)
            return the_force


class TwoGaussianPotential(PotentialFunction1D):

    def __init__(self):
        """
        Initializes the 1-D potential energy surface for a two gaussian
        function. Parameters are set to zero so there are no adjustments to
        default barriers and rare events are initialized as well.
        """
        self.parameters = np.array([-5, 10])
        self.rare_event = np.array([-2.0, float("inf")])

    def get_potential(self, coords):
        """
        Returns the potential energy at a coordinate.

        Parameters:
        -----------
        coords      :   float
                        coordinate of walker

        Returns:
        ----------
        V           :   float
                        potential energy at the present location.
        """
        V = (self.parameters[0] * np.exp(-(coords - 2/0.75)**2) -
             self.parameters[1] * np.exp(-(coords + 2/0.75)**2))
        return V

    def get_force(self, coords):
        """
        Returns the force acting on the walker at a coordinate.

        Parameters:
        -----------
        coords      :   float
                        coordinate of walker

        Returns:
        ----------
        the_force   :   float
                        forces at the present location.
        """
        the_force = ((self.parameters[0] * 2 * -1 * (coords - 2/0.75) *
                      np.exp(-(coords - 2/0.75)**2) - self.parameters[1] * 2 *
                      (-1) * (coords + 2/0.75) *
                      np.exp(-(coords + 2/0.75)**2)) * (-1))
        return the_force


class MullerBrownPotential(PotentialFunction2D):

    def __init__(self):
        """
        Initializes the 2-D potential energy surface for a Muller Brown
        Potential. Parameters are set to zero so there are no adjustments to
        default barriers and rare events are initialized as well.
        """
        self.parameters = np.array([106])
        self.rare_event = np.array([-0.5, 1, 0.5, 0.1])

    def get_potential(self, coords):
        """
        Returns the potential energy at a coordinate.

        Parameters:
        -----------
        coords      :   array of floats
                        coordinate of walker

        Returns:
        ----------
        V           :   float
                        potential energy at the present location.
        """
        Eb = self.parameters[0]
        A = np.array([-200.0*Eb/106, -100.0*Eb/106, -170.0*Eb/106,
                      15.0*Eb/106])
        a = np.array([-1.0, -1.0, -6.50, 0.7])
        b = np.array([0.0, 0.0, 11.0, 0.6])
        c = np.array([-10.0, -10.0, -6.50, 0.7])
        x0 = np.array([1.0, 0.0, -0.50, -1.0])
        y0 = np.array([0.0, 0.5, 1.50, 1.0])
        x = coords[0]
        y = coords[1]
        V = np.zeros([y.size, x.size])
        if hasattr(x, "__len__") is True:
            for k in range(0, y.size - 1):
                for j in range(0, x.size - 1):
                    V[k, j] = sum(A * np.exp(a*(x[j]-x0)**2 +
                                             b*(x[j]-x0)*(y[k]-y0) +
                                             c*(y[k]-y0)**2))
        else:
            V = sum(A * np.exp(a * (x-x0)**2 + b * (x-x0)*(y-y0) +
                               c * (y-y0)**2))
        return V

    def get_force(self, coords):
        """
        Returns the force acting on the walker at a coordinate.

        Parameters:
        -----------
        coords      :   array of floats
                        coordinate of walker

        Returns:
        ----------
        Fpot        :   float
                        forces at the present location.
        """
        Eb = self.parameters[0]
        A = np.array([-200.0*Eb/106, -100.0*Eb/106, -170.0*Eb/106,
                      15.0*Eb/106])
        a = np.array([-1.0, -1.0, -6.50, 0.7])
        b = np.array([0.0, 0.0, 11.0, 0.6])
        c = np.array([-10.0, -10.0, -6.50, 0.7])
        x0 = np.array([1.0, 0.0, -0.50, -1.0])
        y0 = np.array([0.0, 0.5, 1.50, 1.0])
        x = coords[0]
        y = coords[1]
        Fpotx = (200.0/53.0*np.exp(-(x-1)**2-10*y**2)*Eb*(x-1) +
                 100.0/53.0*np.exp(-x**2-10*(y-0.5)**2)*Eb*x -
                 85.0/53.0*np.exp(-6.5*(0.5+x)**2 +
                                  11*(0.5+x)*(y-1.5) -
                                  6.5*(y-1.5)**2)*Eb*(-13*(x+0.5)+11*(y-1.5)) +
                 15.0/106.0*np.exp(0.7*(1+x)**2 +
                                   0.6*(1+x)*(y-1) +
                                   0.7*(y-1)**2)*Eb*(1.4*(1+x)+0.6*(y-1)))*-1

        Fpoty = ((-85.0/53.0)*np.exp(-6.5*(0.5+x)**2 +
                                     11.0*(0.5+x)*(y-1.5) -
                                     6.5*(-1.5+y)**2)*Eb*(11*(x+0.5) -
                                                          13*(y-1.5)) +
                 (15.0/106.0)*np.exp(0.7*(1+x)**2+0.6*(1+x)*(y-1) +
                                     0.7*(y-1)**2)*Eb*(0.6*(1+x)+1.4*(y-1)) +
                 (1000.0/53.0)*np.exp(-x**2-10*(y-0.5)**2)*Eb*(y-0.5) +
                 (2000.0/53.0)*np.exp(-(x-1)**2-10*y**2)*Eb*y
                 )*-1
        Fpot = np.array([Fpotx, Fpoty])
        return Fpot

    def get_triggered(self, coords):
        """
        Determines if rare event has happened, and specifies which one if
        there are multiple
        Parameters:
        -----------
        coords    :      float
                         coordinates of walker.

        Returns:
        ----------
        trigger   :      boolean
                         True if rare event has happended, otherwise False

        event     :      String
                         Distinguishes which path/rare event occurs if multiple
                         paths are present
        """
        trigger = False
        event = "null"
        if coords[0] < self.rare_event[0] and coords[1] > self.rare_event[1]:
            trigger = True
            event = "A"
        elif coords[0] > self.rare_event[2] and coords[1] < self.rare_event[3]:
            trigger = True
            event = "B"
        return trigger, event


class CHLCPotential(PotentialFunction2D):

    def __init__(self):
        """
        Initializes the 2-D potential energy surface for a CH3l + Cl potential.
        Parameters are set to zero so there are no adjustments to
        default barriers and rare events are initialized as well.
        Imports the potential energy function.
        """
        self.parameters = np.array([0, 0])
        self.rare_event = np.array([0.20E-9, 0.25E-9,
                                    float('inf'), float("inf")*-1])
        self.pot = pickle.load(open("final_c_chl_potential.p", "rb"))
        self.fx = pickle.load(open("final_c_chl_fx.p", "rb"))
        self.fy = pickle.load(open("final_c_chl_fy.p", "rb"))

    def get_potential(self, coords):
        """
        Returns the potential energy at a coordinate.

        Parameters:
        -----------
        coords      :   array of floats
                        coordinate of walker

        Returns:
        ----------
        V           :   float
                        potential energy at the present location.
        """
        x = coords[0]
        y = coords[1]
        if hasattr(x, "__len__") is True:
            V = np.zeros([y.size, x.size])
            for k in range(0, x.size - 1):
                for j in range(0, y.size - 1):
                    V[j, k] = self.pot(y[j], x[k])
        else:

            V = self.pot(y, x)
        return V

    def get_force(self, coords):
        """
        Returns the force acting on the walker at a coordinate.

        Parameters:
        -----------
        coords      :   array of floats
                        coordinate of walker

        Returns:
        ----------
        Fpot        :   array of floats
                        forces at the present location.
        """
        x = coords[0]
        y = coords[1]
        # pdb.set_trace()
        Fx = self.fx(y, x)
        Fy = self.fy(y, x)
        Fpotx = Fx
        Fpoty = Fy
        Fpot = np.array([Fpotx, Fpoty])
        # print Fpot
        return Fpot

    def get_triggered(self, coords):
        """
        Determines if rare event has happened, and specifies which one if
        there are multiple
        Parameters:
        -----------
        coords    :      array of floats
                         coordinates of walker.

        Returns:
        ----------
        trigger   :      boolean
                         True if rare event has happended, otherwise False

        event     :      String
                         Distinguishes which path/rare event occurs if multiple
                         paths are present
        """
        trigger = False
        event = "null"
        if coords[0] > self.rare_event[0] and coords[1] < self.rare_event[1]:
            trigger = True
            event = "A"
        elif coords[0] > self.rare_event[2] and coords[1] < self.rare_event[3]:
            trigger = True
            event = "B"
        return trigger, event


class CyclePotential(PotentialFunction2D):

    def __init__(self):
        """
        Initializes the 2-D potential energy surface for a cyclization reaction
        potential. Parameters are set to zero so there are no adjustments to
        default barriers and rare events are initialized as well.
        Imports the potential energy function.
        """
        self.parameters = np.array([0, 0])
        self.rare_event = np.array([0.10E-9, 0.2E-9,
                                   float('inf'), float("inf")*-1])
        self.pot = pickle.load(open("cycle_potential.p", "rb"))
        self.fx = pickle.load(open("cycle_fx.p", "rb"))
        self.fy = pickle.load(open("cycle_fy.p", "rb"))

    def get_potential(self, coords):
        """
        Returns the potential energy at a coordinate.

        Parameters:
        -----------
        coords      :   array of floats
                        coordinate of walker

        Returns:
        ----------
        V           :   float
                        potential energy at the present location.
        """
        x = coords[0]
        y = coords[1]
        if hasattr(x, "__len__") is True:
            V = np.zeros([y.size, x.size])
            for k in range(0, x.size - 1):
                for j in range(0, y.size - 1):
                    V[j, k] = self.pot(y[j], x[k])
        else:

            V = self.pot(y, x)
        return V

    def get_force(self, coords):
        """
        Returns the force acting on the walker at a coordinate.

        Parameters:
        -----------
        coords      :   array of floats
                        coordinate of walker

        Returns:
        ----------
        Fpot        :   array of floats
                        forces at the present location.
        """
        x = coords[0]
        y = coords[1]
        # pdb.set_trace()
        Fx = self.fx(y, x)
        Fy = self.fy(y, x)
        Fpotx = Fx
        Fpoty = Fy
        Fpot = np.array([Fpotx, Fpoty])
        # print Fpot
        return Fpot

    def get_triggered(self, coords):
        """
        Determines if rare event has happened, and specifies which one if
        there are multiple
        Parameters:
        -----------
        coords    :      array of floats
                         coordinates of walker.

        Returns:
        ----------
        trigger   :      boolean
                         True if rare event has happended, otherwise False

        event     :      String
                         Distinguishes which path/rare event occurs if multiple
                         paths are present
        """
        trigger = False
        event = "null"
        if coords[0] < self.rare_event[0] and coords[1] > self.rare_event[1]:
            trigger = True
            event = "A"
        elif coords[0] > self.rare_event[2] and coords[1] < self.rare_event[3]:
            trigger = True
            event = "B"
        return trigger, event


def get_potential_dict():
    """Return a dictionary of all of the available potential functions."""
    potential_dict = {'cosine_potential': CosinePotential,
                      'two_gaussian_potential': TwoGaussianPotential,
                      'muller_brown_potential': MullerBrownPotential,
                      'piece_wise_cosine': PieceWiseCosinePotential,
                      'c_chl_potential': CHLCPotential,
                      'cycle_potential': CyclePotential
                      }

    return potential_dict


def get_GUI_presets_dict():
    """Return a dictionary of all of the available potential functions."""
    preset_dict = {'cosine_potential': np.array([3.14, -6.28, 12.57, 0.01, 0,
                                                 0, 0, 0]).astype(str),
                   'two_gaussian_potential': np.array([2.67, -4, 4, 0.01,
                                                       0, 0, 0,
                                                       0]).astype(str),
                   'pv_2D_potential': np.array([1.5, 0, 3.0, 0.01, 0.6, -2.0,
                                                2.0, 0.01]).astype(str),
                   'muller_brown_potential': np.array([0, 0, 0, 0, 0, 0, 0,
                                                       0]).astype(str),
                   'C_Cl_potential': np.array([0, 0, 0, 0, 0, 0, 0,
                                               0]).astype(str)
                   }
    return preset_dict
