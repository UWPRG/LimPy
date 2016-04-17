import numpy as np

class PotentialFunction1D:

    dimension='1-D Potential'
    def __init__(self):
        self.parameters=np.array([0,0])
        self.rare_event = np.array([float("inf")*-1,float("inf")])
        #self.dimension='1-D Potential'
    def get_potential(self, coords):
        pass

    def set_parameters(self, values):
        self.parameters = values

    def get_force(self, coords):
        pass

    def set_lower_rare_event(self, lower):
        self.rare_event[0] = lower

    def set_upper_rare_event(self, upper):
        self.rare_event[1] = upper

    def get_triggered(self, coords):
        trigger = False
        event = "null"
        if coords < rare_event[0]:
            trigger = True
            event = "A"
        elif coords > rare_event[1]:
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
        self.parameters=np.array([1,0],[1,0])
        #self.dimension='1-D Potential'

    def set_switchpoint(self, divider):
        self.divider=divider

    def get_potential(self, coords):
        if coords < self.divider:
            return self.parameters[0,0]*np.cos(coords)+self.parameters[0,1]

        elif coords < self.divider:
            return self.parameters[1,0]*np.cos(coords)+self.parameters[1,1]

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
        self.parameters=np.array([0,0,0,0])
        self.rare_event = np.array([[float("inf")*-1,float("inf")],
                                    [float("inf")*-1,float("inf")]])
        #self.dimension='1-D Potential'
    def get_potential(self, coords):
        A = np.array([-200.0, -100.0, -170.0, 15.0])
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
                                             c*(y[k]-y0)**2))/4.184
        else:
            V = sum(A * np.exp(a * (x-x0)**2 + b * (x-x0)*(y-y0) +
                               c * (y-y0)**2))/4.184
        return V

    def get_force(self, coords):
        A = np.array([-200.0, -100.0, -170.0, 15.0])
        a = np.array([-1.0, -1.0, -6.50, 0.7])
        b = np.array([0.0, 0.0, 11.0, 0.6])
        c = np.array([-10.0, -10.0, -6.50, 0.7])
        x0 = np.array([1.0, 0.0, -0.50, -1.0])
        y0 = np.array([0.0, 0.5, 1.50, 1.0])
        x = coords[0]
        y = coords[1]
        Fpotx = (-400*np.exp(-1*(-1+x)**2 - 10*y**2)*(-1+x) -
                 200*np.exp(-(x**2)-10*(-0.5+y)**2)*x +
                 170*np.exp(-6.5*(0.5+x)**2+11*(0.5+x)*(-1.5+y) -
                 6.5*(-1.5+y)**2)*(-13*(0.5+x)+11*(-1.5+y)) -
                 15*np.exp(0.7*(1+x)**2+0.6*(1+x)*(y-1) +
                 0.7*(y-1)**2)*(1.4*(1+x)+0.6*(y-1)))/4.184
        Fpoty = (170*np.exp(-6.5*(0.5+x)**2 +
                 11*(0.5+x)*(-1.5+y) -
                 6.5*(y-1.5)**2)*(11*(0.5+x)-13*(y-1.5)) -
                 15*np.exp(0.7*(1+x)**2+0.6*(1+x)*(y-1) +
                 0.7*(y-1)**2)*(0.6*(x+1)+1.4*(y-1)) -
                 2000*np.exp(-x**2-10*(y-0.5)**2)*(y-0.5) -
                 4000*np.exp(-1*(x-1)**2-10*y**2)*y)/4.184
        Fpot = np.array([Fpotx, Fpoty])
        return Fpot

    def get_triggered(self, coords):
        return (False, 'Null')

def get_potential_dict():
    """Return a dictionary of all of the available potential functions."""
    potential_dict = {'cosine_potential': CosinePotential,
                      'two_gaussian_potential': TwoGaussianPotential,
                      'muller_brown_potential': MullerBrownPotential
                     }


    return potential_dict
