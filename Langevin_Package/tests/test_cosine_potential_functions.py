"""Unit tests for potential functions"""

import unittest
import numpy as np
from ..potential_functions import cosine_potential, cosine_potential_bc


coords = np.array([0.5, 0.5])
vnew = 0.5
f2 = 1.0


class TestCosinePotentialFunctions(unittest.TestCase):
    """Tests for the various potential functions"""

    def test_cosine_correct_force(self):
        """Check cosine potential returns correct force value"""

        f = cosine_potential(coords[0])[1]
        self.assertEqual(f, np.sin(coords[0])*2.5)

    def test_cosine_no_trigger(self):
        """Check cosine potential doesn't trigger at non-rare event"""

        trigger = cosine_potential(coords[0])[2]
        self.assertFalse(trigger)

    def test_cosine_trigger(self):
        """Check cosine potential trigger activates at rare event"""

        trigger1 = cosine_potential(coords[0]-np.pi)[2]
        self.assertTrue(trigger1)
        trigger2 = cosine_potential(coords[0]+2*np.pi)[2]
        self.assertTrue(trigger2)

    def test_cosine_correct_potential(self):
        """Check cosine potential returns correct potential value"""

        vpot = cosine_potential(coords[0])[0]
        self.assertEqual(vpot, np.cos(coords[0])*2.5)

    def test_cosine_correct_bc_coords_nochange(self):
        """Check cosine boundary condition is only applied at boundary"""

        newcoords = cosine_potential_bc(vnew, f2, coords[0])[2]
        self.assertEqual(newcoords, coords[0])

    def test_cosine_correct_bc_coords_change(self):
        """Check cosine boundary condition is periodic"""

        newcoords = cosine_potential_bc(vnew, f2, coords[0]-2*np.pi)[2]
        self.assertEqual(newcoords, coords[0])

        newcoords = cosine_potential_bc(vnew, f2, coords[0]+2*np.pi)[2]
        self.assertEqual(newcoords, coords[0])

    def test_cosine_correct_bc_potential_nochange(self):
        """Check cosine boundary condition doesn't affect potential"""

        newcoords = cosine_potential_bc(vnew, f2, coords[0]+2*np.pi)[0]
        self.assertEqual(newcoords, vnew)
        newcoords = cosine_potential_bc(vnew, f2, coords[0]-2*np.pi)[0]
        self.assertEqual(newcoords, vnew)
        newcoords = cosine_potential_bc(vnew, f2, coords[0])[0]
        self.assertEqual(newcoords, vnew)

    def test_cosine_correct_bc_force_nochange(self):
        """Check cosine boundary condition doesn't affect force"""

        newcoords = cosine_potential_bc(vnew, f2, coords[0]+2*np.pi)[1]
        self.assertEqual(newcoords, f2)
        newcoords = cosine_potential_bc(vnew, f2, coords[0]-2*np.pi)[1]
        self.assertEqual(newcoords, f2)
        newcoords = cosine_potential_bc(vnew, f2, coords[0])[1]
        self.assertEqual(newcoords, f2)

    def test_cosine_correct_bc_no_bcbias(self):
        """Check cosine boundary condition applies no additional bias"""

        newcoords = cosine_potential_bc(vnew, f2, coords[0]+2*np.pi)[3]
        self.assertEqual(newcoords, 0)
        newcoords = cosine_potential_bc(vnew, f2, coords[0]-2*np.pi)[3]
        self.assertEqual(newcoords, 0)
        newcoords = cosine_potential_bc(vnew, f2, coords[0])[3]
        self.assertEqual(newcoords, 0)
