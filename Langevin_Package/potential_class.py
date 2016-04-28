"""Defines the Potential Functions to be used as a class"""
import numpy as np
try:
    import cPickle as pickle
except:
    import pickle
import pdb


class PotentialFunction1D:
    """Basic Outline of a 1-D Potential Function"""
    dimension='1-D Potential'
    def __init__(self):
        self.parameters=np.array([0,0])
        self.rare_event = np.array([float("inf")*-1,float("inf")])

    def get_potential(self, coords):
        """
        Equation that returns the potential energy at a CV
        To be specified by child class
        """
        pass

    def set_parameters(self, values):
        """Overwrites default parameters for potential function"""
        self.parameters = values

    def get_force(self, coords):
        """
        Equation that returns the force at a CV
        To be specified by child class
        """
        pass

    def set_lower_rare_event(self, lower):
        """Overwrites default lower rare event for potential function"""
        self.rare_event[0] = lower

    def set_upper_rare_event(self, upper):
        """Overwrites default upper rare event for potential function"""
        self.rare_event[1] = upper

    def get_triggered(self, coords):
        """
        Determines if rare event has happened, and specifies which one if
        there are multiple
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

    dimension='2-D Potential'
    def __init__(self):
        self.parameters=np.array([0,0,0,0])
        self.rare_event = np.array([[float("inf")*-1,float("inf")],
                                    [float("inf")*-1,float("inf")]])
        #self.dimension='1-D Potential'
    def get_potential(self, coords):
        pass

    def set_parameters(self, values):
        self.parameters = values

    def get_force(self, coords):
        pass

    def set_lower_rare_event(self, lowerx, lowery):
        self.rare_event[0][0] = lowerx
        self.rare_event[0][1] = lowery

    def set_upper_rare_event(self, upperx, uppery):
        self.rare_event[1][0] = upperx
        self.rare_event[1][1] = uppery

    def get_triggered(self, coords):
        pass


class CosinePotential(PotentialFunction1D):

    def __init__(self):
        self.parameters=np.array([1,0])
        self.rare_event = np.array([-2.0, 8.28])
        #self.dimension='1-D Potential'
    def get_potential(self, coords):
        return self.parameters[0]*np.cos(coords)+self.parameters[1]

    def get_force(self, coords):
        return self.parameters[0]*np.sin(coords)

class PieceWiseCosinePotential(PotentialFunction1D):

    def __init__(self):
        self.parameters=np.array([1,0, 1,0])
        self.rare_event = np.array([-2.0, 8.28])

    def get_potential(self, coords):

        if hasattr(coords, "__len__") is True:
            V = np.zeros_like(coords)
            for i in range(0,len(coords)):
                if coords[i] < np.pi:
                    V[i] = (self.parameters[0]*np.cos(coords[i])+
                            self.parameters[1])

                elif coords[i] > np.pi:
                    V[i] = (self.parameters[2]*np.cos(coords[i])+
                            self.parameters[3])
            return V
        else:
            if coords < np.pi:
                return self.parameters[0]*np.cos(coords)+self.parameters[1]

            elif coords > np.pi:
                return self.parameters[2]*np.cos(coords)+self.parameters[3]

    def get_force(self, coords):
        return self.parameters[0]*np.sin(coords)

class TwoGaussianPotential(PotentialFunction1D):

    def __init__(self):
        self.parameters=np.array([-5,10])
        self.rare_event = np.array([-2.0, float("inf")])
    def get_potential(self, coords):
        return (self.parameters[0]* np.exp(-(coords - 2/0.75)**2) -
                self.parameters[1]*np.exp(-(coords + 2/0.75)**2))

    def get_force(self, coords):
        return ((self.parameters[0] * 2 * -1 * (coords - 2/0.75) *
                 np.exp(-(coords - 2/0.75)**2) -  self.parameters[1] * 2 *
                 -1 * (coords + 2/0.75) * np.exp(-(coords + 2/0.75)**2)) * (-1))

class MullerBrownPotential(PotentialFunction2D):

    def __init__(self):
        self.parameters=np.array([106])
        self.rare_event = np.array([-0.05,1,0.5,0.0])
        #self.dimension='1-D Potential'
    def get_potential(self, coords):
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
        self.parameters=np.array([0,0])
        self.rare_event = np.array([0.20E-9,0.25E-9,float('inf'),float("inf")*-1])
        self.pot = pickle.load(open("c_chl_potential.p", "rb"))
        self.fx = pickle.load(open("c_chl_fdx.p", "rb"))
        self.fy = pickle.load(open("c_chl_fdy.p", "rb"))
        #self.dimension='1-D Potential'
    def get_potential(self, coords):
        x = coords[0]*1E9
        y = coords[1]*1E9
        if hasattr(x, "__len__") is True:
            V = np.zeros([x.size, y.size])
            for k in range(0, x.size - 1):
                for j in range(0, y.size - 1):
                    V[k, j] = self.pot(x[k], y[j])
        else:

            V= self.pot(x, y)
        return V

    def get_force(self, coords):

        x = coords[0]*1E9
        y = coords[1]*1E9
        # pdb.set_trace()
        Fx = self.fx(x, y)
        Fy = self.fy(x, y)
        Fpotx = Fx * -1
        Fpoty = Fy * -1
        Fpot = np.array([Fpotx, Fpoty])
        # print Fpot
        return Fpot

    def get_triggered(self, coords):
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
                      'c_chl_potential': CHLCPotential
                     }


    return potential_dict

def get_GUI_presets_dict():
    """Return a dictionary of all of the available potential functions."""
    preset_dict = {'cosine_potential': np.array([3.14,-6.28,12.57,0.01,0,
                                                    0,0,0]).astype(str),
                      'two_gaussian_potential': np.array([2.67,-4,4,0.01,
                                                          0,0,0,0]).astype(str),
                      'pv_2D_potential': np.array([1.5,0,3.0,0.01,0.6,-2.0,
                                                   2.0,0.01]).astype(str),
                      'muller_brown_potential': np.array([0,0,0,0,0,0,0,
                                                          0]).astype(str),
                      'C_Cl_potential': np.array([0,0,0,0,0,0,0,0]).astype(str)
                    }
    return preset_dict
