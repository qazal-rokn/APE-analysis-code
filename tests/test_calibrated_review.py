"""Synthetic reference-equivalence test; intentionally opens no Dewesoft data."""
from pathlib import Path
import sys
import unittest
import numpy as np

REFERENCE = Path("/Users/user/Desktop/APE Research/APE-data-analysis")
sys.path.insert(0, str(REFERENCE / "src"))

from ape_analysis.calibration import (
    PolarizationCalibrationParameters, RawPemQuadratureChannel, build_polarization_calibration,
)
from ape_c_lpsd.calibrated_review import calibrated_quadratures


class _ReducedRecord:
    def __init__(self, i0_w, eta0): self.i0_w=i0_w; self.eta0=eta0


class _Quality:
    records=(_ReducedRecord(1.3e-3, .21), _ReducedRecord(1.1e-3, .18))
    combined_good=(True, True)


class CalibratedReviewTests(unittest.TestCase):
 def test_component_norm_matches_frozen_reference_rho_for_valid_samples(self):
  parameters=PolarizationCalibrationParameters()
  x=np.array([.1,-.2,.3,-.4,.5,-.6,.7,-.8],dtype=float)
  y=np.array([.9,.8,-.7,-.6,.5,.4,-.3,-.2],dtype=float)
  quality=_Quality()
  # Repeat the two reduced records at the frozen 1000:1 mapping without a .dxd source.
  x=np.tile(x,250); y=np.tile(y,250)
  valid=np.ones_like(x,dtype=bool)
  d,xe,ye=calibrated_quadratures(x,y,quality,0,parameters,valid)
  reference=build_polarization_calibration(
      RawPemQuadratureChannel("pem1_x",x,20000.,"synthetic",0.),
      RawPemQuadratureChannel("pem1_Y",y,20000.,"synthetic",0.),
      (1.3e-3,1.1e-3),(.21,.18),(True,True),parameters).records
  rho=np.array([record.rho for record in reference],dtype=float)
  # Both paths use IEEE float64 operations; allow only round-off level difference.
  np.testing.assert_allclose(np.hypot(xe,ye),rho,rtol=1e-14,atol=1e-15)
  assert np.all(np.isfinite(d))


if __name__ == "__main__":
 unittest.main()
