"""Synthetic validation for the approved APE calibration equations."""
from pathlib import Path
import unittest

import numpy as np

from ape_c_lpsd.calibration import calibrate_chunk
from ape_c_lpsd.ingestion.config import CalibrationConfig
from ape_c_lpsd.ingestion.models import APEDXDChunk


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.cfg = CalibrationConfig(
            cavity_transmission_w_per_v=0.020,
            responsivity_a_per_w=0.7,
            transimpedance_gain_v_per_a=1.0e6,
            second_harmonic_pd_response=0.7,
        )

    def test_recovers_known_injected_physics(self):
        n = 64
        i0 = np.linspace(0.008, 0.012, n)
        eta = np.linspace(0.020, 0.024, n)
        rho_x = np.linspace(-2e-6, 3e-6, n)
        rho_y = np.linspace(4e-6, -1e-6, n)

        q = self.cfg.responsivity_a_per_w
        g = self.cfg.transimpedance_gain_v_per_a
        h2 = self.cfg.second_harmonic_pd_response
        ctrans = self.cfg.cavity_transmission_w_per_v

        vtrans = i0 / ctrans
        r2f = h2 * q * g * i0 * eta**2 / (2.0 * np.sqrt(2.0))
        x2 = r2f
        y2 = np.zeros_like(r2f)
        x1 = np.sqrt(2.0) * q * g * i0 * eta * rho_x
        y1 = np.sqrt(2.0) * q * g * i0 * eta * rho_y

        chunk = APEDXDChunk(
            source_file=Path("synthetic.dxd"),
            start_index=100,
            stop_index=100 + n,
            sample_rate_hz=20000.0,
            x_1f=x1,
            y_1f=y1,
            x_2f=x2,
            y_2f=y2,
            cavity_transmission=vtrans,
            first_timestamp_s=100 / 20000.0,
            last_timestamp_s=(100 + n - 1) / 20000.0,
        )

        out = calibrate_chunk(chunk, self.cfg)
        self.assertTrue(np.all(out.valid_mask))
        np.testing.assert_allclose(out.i0_w, i0, rtol=1e-14, atol=0.0)
        np.testing.assert_allclose(out.eta0, eta, rtol=1e-14, atol=0.0)
        np.testing.assert_allclose(out.rho_x, rho_x, rtol=1e-14, atol=1e-18)
        np.testing.assert_allclose(out.rho_y, rho_y, rtol=1e-14, atol=1e-18)
        self.assertEqual(out.start_index, 100)
        self.assertEqual(out.stop_index, 164)

    def test_invalid_physical_samples_remain_indexed_and_nan(self):
        n = 4
        chunk = APEDXDChunk(
            source_file=Path("synthetic.dxd"),
            start_index=0,
            stop_index=n,
            sample_rate_hz=20000.0,
            x_1f=np.ones(n),
            y_1f=np.ones(n),
            x_2f=np.array([1.0, 0.0, 1.0, 1.0]),
            y_2f=np.zeros(n),
            cavity_transmission=np.array([0.5, 0.5, 0.0, np.nan]),
            first_timestamp_s=0.0,
            last_timestamp_s=3 / 20000.0,
        )
        out = calibrate_chunk(chunk, self.cfg)
        self.assertEqual(out.valid_mask.tolist(), [True, False, False, False])
        self.assertTrue(np.isfinite(out.rho_x[0]))
        self.assertTrue(np.isnan(out.rho_x[1:]).all())
        self.assertTrue(np.isnan(out.rho_y[1:]).all())
        self.assertTrue(np.isnan(out.eta0[1:]).all())


if __name__ == "__main__":
    unittest.main()
