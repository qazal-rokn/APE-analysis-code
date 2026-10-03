"""Synthetic-only plotting tests; deliberately not executed in this stage."""
import unittest
import numpy as np
from ape_c_lpsd.ingestion.plotting import _envelope
class PlottingTests(unittest.TestCase):
 def test_envelope_preserves_bounds(self):
  i,lo,hi,width=_envelope(np.array([1.,4.,2.,3.]),2)
  self.assertEqual(width,2); self.assertTrue(np.array_equal(lo,[1.,2.])); self.assertTrue(np.array_equal(hi,[4.,3.]))
