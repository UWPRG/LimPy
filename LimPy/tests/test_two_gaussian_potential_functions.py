"""Unit tests for potential functions"""

import unittest
import numpy as np
from ..potential_functions import (two_gaussian_potential,
                                   two_gaussian_potential_bc)


coords = np.array([0.5, 0.5])
vnew = 0.5
f2 = 1.0


class TestTwoGaussianPotentialFunctions(unittest.TestCase):
    """Tests for the various potential functions"""

    def test_two_gaussian_potential_correct_force(self):
        """Check two_gaussian_potential returns correct force value"""

        f = round(two_gaussian_potential(coords[0])[1], 5)
        self.assertEqual(f, 0.19537)

    def test_two_gaussian_potential_no_trigger(self):
        """Check two_gaussian_potential doesn't trigger at non-rare event"""

        trigger = two_gaussian_potential(coords[0])[2]
        self.assertFalse(trigger)

    def test_two_gaussian_potential_trigger(self):
        """Check two_gaussian_potential trigger activates at rare event"""

        trigger2 = two_gaussian_potential(coords[0]-5)[2]
        self.assertTrue(trigger2)

    def test_two_gaussian_potential_correct_potential(self):
        """Check two_gaussian_potential returns correct potential value"""

        vpot = round(two_gaussian_potential(coords[0])[0], 5)
        self.assertEqual(vpot, -0.04617)

    def test_two_gaussian_potential_correct_bc_coords_nochange(self):
        """Check two_gaussian_potential_bc is only applied at boundary"""

        newcoords = two_gaussian_potential_bc(vnew, f2, coords[0]+5)[2]
        self.assertEqual(newcoords, coords[0]+5)

        newcoords = two_gaussian_potential_bc(vnew, f2, coords[0]-5)[2]
        self.assertEqual(newcoords, coords[0]-5)

        newcoords = two_gaussian_potential_bc(vnew, f2, coords[0])[2]
        self.assertEqual(newcoords, coords[0])

    def test_two_gaussian_potential_correct_bc_potential(self):
        """Check two_gaussian_potential_bc doesn't affect Energy"""

        newcoords1 = round(two_gaussian_potential_bc(vnew, f2, 4.5)[0], 5)
        self.assertEqual(newcoords1, 5.40493)
        newcoords2 = round(two_gaussian_potential_bc(vnew, f2, -4.5)[0], 5)
        self.assertEqual(newcoords2, 4.55987)
        newcoords3 = round(two_gaussian_potential_bc(vnew, f2, 0.5)[0], 5)
        self.assertEqual(newcoords3, vnew)

    def test_two_gaussian_potential_correct_bc_force(self):
        """Check two_gaussian_potential_bc doesn't affect force"""

        newcoords1 = round(two_gaussian_potential_bc(vnew, f2, 4.5)[1], 5)
        self.assertEqual(newcoords1, -50.0)
        newcoords2 = round(two_gaussian_potential_bc(vnew, f2, -4.5)[1], 5)
        self.assertEqual(newcoords2, 50.0)
        newcoords3 = round(two_gaussian_potential_bc(vnew, f2, 0.5)[1], 5)
        self.assertEqual(newcoords3, f2)

    def test_two_gaussian_potential_correct_bc_bcbias(self):
        """Check two_gaussian_potential_bc applies no additional bias"""

        newcoords1 = round(two_gaussian_potential_bc(vnew, f2, 4.5)[3], 5)
        self.assertEqual(newcoords1, 5.40493 - vnew)
        newcoords2 = round(two_gaussian_potential_bc(vnew, f2, -4.5)[3], 5)
        self.assertEqual(newcoords2, 4.55987 - vnew)
        newcoords3 = round(two_gaussian_potential_bc(vnew, f2, 0.5)[3], 5)
        self.assertEqual(newcoords3, 0.0)
