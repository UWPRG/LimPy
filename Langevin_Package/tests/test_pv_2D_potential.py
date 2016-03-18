"""Unit tests for potential functions"""

import unittest
import numpy as np
import pdb
from ..potential_functions import (pv_2D_potential,
                                   pv_2D_potential_bc)


coords = np.array([0.5, 0.5])
vnew = 0.5
f2 = np.array([1.0, 1.0])


class TestPotentialFunctions(unittest.TestCase):
    """Tests for the various potential functions"""

    def test_pv_2D_potential_correct_force(self):
        """Check pv_2D_potential returns correct force value"""

        f1 = round(pv_2D_potential(coords[0], coords[1])[1][0], 5)
        self.assertEqual(f1, -1.36035)
        f2 = round(pv_2D_potential(coords[0], coords[1])[1][1], 5)
        self.assertEqual(f2, 0.85841)

    def test_pv_2D_potential_no_trigger(self):
        """Check pv_2D_potential doesn't trigger at non-rare event"""

        trigger = pv_2D_potential(1.50, 0.60)[2]
        self.assertFalse(trigger)

    def test_pv_2D_potential_trigger(self):
        """Check pv_2D_potential trigger activates at rare event"""
        trigger2 = pv_2D_potential(float(0.50), float(0.50))[2]
        self.assertTrue(trigger2)
        trigger2 = pv_2D_potential(float(2.50), float(0.50))[2]
        self.assertTrue(trigger2)

    def test_pv_2D_potential_correct_potential(self):
        """Check pv_2D_potential returns correct potential value"""

        vpot = round(pv_2D_potential(coords[0], coords[1])[0], 4)
        self.assertEqual(vpot, -2.5896)

    def test_pv_2D_potential_correct_bc_coords_change(self):
        """Check pv_2D_potential_bc is only applied at boundary"""

        newcoords1 = pv_2D_potential_bc(vnew, f2, np.array([3.7, 0.5]))[2]
        self.assertEqual(round(newcoords1[0], 5), 0.7)
        self.assertEqual(newcoords1[1], coords[1])

        newcoords2 = pv_2D_potential_bc(vnew, f2, np.array([-2.7, 0.5]))[2]
        self.assertEqual(round(newcoords2[0], 5), 0.3)
        self.assertEqual(newcoords2[1], coords[1])

        newcoords3 = pv_2D_potential_bc(vnew, f2, coords)[2]
        self.assertEqual(newcoords3[0], coords[0])
        self.assertEqual(newcoords3[1], coords[1])

        newcoords4 = pv_2D_potential_bc(vnew, f2, np.array([0.5, 3.5]))[2]
        self.assertEqual(newcoords4[0], coords[0])
        self.assertEqual(newcoords4[1], 3.5)

        newcoords5 = pv_2D_potential_bc(vnew, f2, np.array([0.5, -2.5]))[2]
        self.assertEqual(newcoords5[0], coords[0])
        self.assertEqual(newcoords5[1], -2.5)

    def test_pv_2D_potential_correct_bc_potential(self):
        """Check pv_2D_potential_bc doesn't affect Energy"""

        newv1 = pv_2D_potential_bc(vnew, f2, np.array([3.5, 0.5]))[0]
        self.assertEqual(newv1, vnew)
        newv2 = pv_2D_potential_bc(vnew, f2, np.array([-2.5, 0.5]))[0]
        self.assertEqual(newv2, vnew)
        newv3 = pv_2D_potential_bc(vnew, f2, np.array([0.5, 0.5]))[0]
        self.assertEqual(newv3, vnew)

    def test_pv_2D_potential_correct_bc_force(self):
        """Check pv_2D_potential_bc doesn't affect force"""
        newf1x = pv_2D_potential_bc(vnew, f2, np.array([3.5, 0.5]))[1][0]
        self.assertEqual(newf1x, f2[0])
        newf2x = pv_2D_potential_bc(vnew, f2, np.array([-2.5, 0.5]))[1][0]
        self.assertEqual(newf2x, f2[0])
        newf3x = pv_2D_potential_bc(vnew, f2, np.array([0.5, 0.5]))[1][0]
        self.assertEqual(newf3x, f2[0])

        newf1y = pv_2D_potential_bc(vnew, f2, np.array([3.5, 0.5]))[1][1]
        self.assertEqual(newf1y, f2[1])
        newf2y = pv_2D_potential_bc(vnew, f2, np.array([-2.5, 0.5]))[1][1]
        self.assertEqual(newf2y, f2[1])
        newf3y = pv_2D_potential_bc(vnew, f2, np.array([0.5, 0.5]))[1][1]
        self.assertEqual(newf3y, f2[1])

    def test_pv_2D_potential_correct_bc_bcbias(self):
        """Check pv_2D_potential_bc applies no additional bias"""

        bcbias1 = pv_2D_potential_bc(vnew, f2, np.array([3.5, 0.5]))[3]
        self.assertEqual(bcbias1, 0.0)
        bcbias2 = pv_2D_potential_bc(vnew, f2, np.array([-2.5, 0.5]))[3]
        self.assertEqual(bcbias2, 0.0)
        bcbias3 = pv_2D_potential_bc(vnew, f2, np.array([0.5, 0.5]))[3]
        self.assertEqual(bcbias3, 0.0)
