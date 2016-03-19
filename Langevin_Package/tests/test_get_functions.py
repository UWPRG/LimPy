"""Unit tests for potential functions"""

import unittest
from ..potential_functions import (get_potential_dict,
                                   get_boundary_condition_dict)


class TestGetFunctions(unittest.TestCase):

    def test_get_potential_dict(self):
        """Check that dictionary of potentials can be retrieved"""
        pot_dict = get_potential_dict()
        self.assertTrue(type(pot_dict) == dict)
        self.assertTrue(len(pot_dict) > 1)

    def test_get_boundary_condition_dict(self):
        """Check that dictionary of bc's can be retrieved"""
        bc_dict = get_boundary_condition_dict()
        self.assertTrue(type(bc_dict) == dict)
        self.assertTrue(len(bc_dict) > 1)
