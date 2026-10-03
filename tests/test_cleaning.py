"""Synthetic validation for staged local-median APE cleaning."""
from pathlib import Path
import unittest
import numpy as np

from ape_c_lpsd.cleaning import clean_chunk
from ape_c_lpsd.calibration import calibrate_chunk
from ape_c_lpsd.ingestion.config import CalibrationConfig, CleaningConfig
from ape_c_lpsd.ingestion.models import APEDXDChunk


class CleaningTests(unittest.TestCase):
    def setUp(self):
        self.fs = 100.0
        self.cleaning = CleaningConfig(
            block_size_seconds=1.0,
            reference_window_seconds=20.0,
            cavity_lower_fraction=0.90,
            eta0_lower_fraction=0.90,
        )
        self.calibration = CalibrationConfig(
            cavity_transmission_w_per_v=0.020,
            responsivity_a_per_w=0.7,
            transimpedance_gain_v_per_a=1.0e6,
            second_harmonic_pd_response=0.7,
        )

    def _chunk(self):
        n = 12000
        t = np.arange(n) / self.fs
        vtrans = 0.50 - 0.04 * (t / t[-1])
        eta = 0.022 - 0.001 * (t / t[-1])

        cavity_glitch = slice(3500, 3600)
        eta_glitch = slice(8000, 8100)
        vtrans[cavity_glitch] *= 0.50
        eta[eta_glitch] *= 0.50

        i0 = self.calibration.cavity_transmission_w_per_v * vtrans
        q = self.calibration.responsivity_a_per_w
        g = self.calibration.transimpedance_gain_v_per_a
        h2 = self.calibration.second_harmonic_pd_response
        r2f = h2 * q * g * i0 * eta**2 / (2.0 * np.sqrt(2.0))

        rho_x = np.full(n, 2.0e-6)
        rho_y = np.full(n, -3.0e-6)
        x1 = np.sqrt(2.0) * q * g * i0 * eta * rho_x
        y1 = np.sqrt(2.0) * q * g * i0 * eta * rho_y

        return APEDXDChunk(
            source_file=Path("synthetic.dxd"),
            start_index=0,
            stop_index=n,
            sample_rate_hz=self.fs,
            x_1f=x1,
            y_1f=y1,
            x_2f=r2f,
            y_2f=np.zeros(n),
            cavity_transmission=vtrans,
            first_timestamp_s=0.0,
            last_timestamp_s=(n - 1) / self.fs,
        ), cavity_glitch, eta_glitch

    def test_rejects_downward_glitches_but_preserves_slow_drift(self):
        chunk, cavity_glitch, eta_glitch = self._chunk()
        result = clean_chunk(chunk, self.cleaning, self.calibration)
        self.assertFalse(np.any(result.cavity_good[cavity_glitch]))
        self.assertTrue(np.all(result.cavity_good[eta_glitch]))
        self.assertFalse(np.any(result.eta_good[eta_glitch]))

        safe = np.ones(chunk.n_samples, dtype=bool)
        safe[cavity_glitch] = False
        safe[eta_glitch] = False
        self.assertGreater(np.mean(result.final_good[safe]), 0.999)

    def test_full_calibration_is_gated_by_final_cleaning_mask(self):
        chunk, cavity_glitch, eta_glitch = self._chunk()
        result = clean_chunk(chunk, self.cleaning, self.calibration)
        calibrated = calibrate_chunk(
            chunk,
            self.calibration,
            input_valid_mask=result.final_good,
        )
        self.assertTrue(np.isnan(calibrated.rho_x[cavity_glitch]).all())
        self.assertTrue(np.isnan(calibrated.rho_x[eta_glitch]).all())
        self.assertTrue(np.isfinite(calibrated.rho_x[result.final_good]).all())


if __name__ == "__main__":
    unittest.main()
