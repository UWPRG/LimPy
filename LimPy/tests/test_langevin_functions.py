import unittest

import pandas as pd
import numpy as np
import os
import os.path as op

from ..potential_functions import (get_potential_dict,
                                   get_boundary_condition_dict)
from ..langevin_functions import get_parameters

oneDinps = np.array([1000000, 0.01, 2.666, 300, 1, -4, 4, 0.01, 0.5, -2, 2,
                     0.01, 5])
oneDmdps = np.array([0.25, 0.1, 250, 10000, 100])
oneDdimension = '1-D Potential'
oneDmethod = 'Infrequent WT MetaD'
oneDpotfunc = 'two_gaussian_potential'
oneDfiletitle = '1dcheck'
oneDmakeplot = 'False'
path = op.split(os.getcwd())[0]
path = op.join(path, 'sample_inputs/inputs2gaussian.csv')


class TestLangevinFunctions(unittest.TestCase):
    """Tests for the various langevin functions"""

    def test_get_parameters(self):
        """Test that datafile is read in correctly"""
        received = get_parameters(path)
        np.testing.assert_almost_equal(received[0], oneDinps, decimal=5)
        np.testing.assert_almost_equal(received[1], oneDmdps, decimal=5)
        self.assertEqual(received[2], oneDdimension)
        self.assertEqual(received[3], oneDmethod)
        self.assertEqual(received[4], oneDpotfunc)
        self.assertEqual(received[5], oneDfiletitle)
        self.assertEqual(received[6], oneDmakeplot)
