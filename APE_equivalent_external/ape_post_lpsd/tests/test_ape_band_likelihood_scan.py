"""Regression tests for band-wide mass selection and resume identity."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ape_band_likelihood_scan import selected_masses, config_hash

class BandScanTests(unittest.TestCase):
    def test_selection(self):
        f=selected_masses(64,65.1,8000.)
        self.assertEqual(len(f),64)
        self.assertAlmostEqual(float(f[0]),65.1)
        self.assertAlmostEqual(float(f[-1]),8000.)
        self.assertTrue(np.all(np.diff(f)>0))
    def test_configuration_hash(self):
        f=selected_masses(8,65.1,8000.)
        self.assertEqual(config_hash(f,{'seeds':25}),config_hash(f,{'seeds':25}))
        self.assertNotEqual(config_hash(f,{'seeds':25}),config_hash(f,{'seeds':30}))
    def test_invalid_bounds(self):
        with self.assertRaises(ValueError):selected_masses(1,65.1,8000.)
        with self.assertRaises(ValueError):selected_masses(8,65.,8000.)
