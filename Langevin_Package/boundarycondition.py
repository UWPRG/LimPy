import numpy as np
import pdb
class BoundaryCondition:

    def __init__(self):
        self.type = "BoundaryCondition"
        self.location = np.array([float('inf')*-1,float('inf')])
    def get_potential(self, coords, origin_potential):
        return origin_potential

    def get_force(self, coords, origin_force):
        return origin_force

    def get_new_location(self, coords):
        return coords

    def add_depositions(self, coords, history):
        return history

    def get_bc_bias(self, coords):
        return 0.0

    def set_new_lower_location(self, lower_location):
        self.location[0] = lower_location

    def set_new_upper_location(self, upper_location):
        self.location[1] = upper_location

class QuarticBoundary(BoundaryCondition):

    def __init__(self):
        self.type='Quartic'
        self.location = np.array([float('inf')*-1,float('inf')])
    def get_potential(self, coords, origin_potential):
        bcvalue = 0.0
        if coords > self.location[1]:
             bcvalue = 100.0*(coords - self.location[1])**4.0
        elif coords < self.location[0]:
             bcvalue = 100.0*(coords - self.location[0])**4.0

        return (bcvalue + origin_potential)

    def get_force(self, coords,  origin_force):
        if coords > self.location[1]:
            origin_force = -400.0*(coords - self.location[1])**3.0
        if coords < self.location[0]:
            origin_force = -400.0*(coords - self.location[0])**3.0
        return origin_force

    def get_bc_bias(self, coords):
        return 100.0*(coords - self.location)**4.0

class PeriodicBoundary(BoundaryCondition):

    def __init__(self):
        self.type='Periodic'
        self.location = np.array([float('inf')*-1,float('inf')])
    def get_new_location(self, coords):
        movement = self.location[1]-self.location[0]
        if coords > self.location[1]:

            newcoords = coords - movement
        elif coords < self.location[0]:

            newcoords = coords + movement
        else:
            newcoords = coords
        return newcoords
    def add_depositions(self, coords, history):
        movement = self.location[1]-self.location[0]
        history = np.append(history,(coords+movement))
        history = np.append(history,(coords-movement))
        return history
def get_boundary_condition_dict():
    """Return a dictionary of all of the available potential functions."""
    bc_dict = {'Quartic': QuarticBoundary,
               'Periodic': PeriodicBoundary,
               'Empty': BoundaryCondition
               }

    return bc_dict
