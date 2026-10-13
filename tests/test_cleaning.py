"""Synthetic validation for approved cavity + protected-eta0 APE cleaning."""
from pathlib import Path
import unittest

import numpy as np

from ape_c_lpsd.calibration import calibrate_chunk
from ape_c_lpsd.cleaning import (
    build_cavity_reference_plan,
    build_cleaning_reference_plan,
    clean_chunk,
    clean_chunk_with_plan,
)
from ape_c_lpsd.ingestion.config import CalibrationConfig, CleaningConfig
from ape_c_lpsd.ingestion.models import APEDXDChunk


class CleaningTests(unittest.TestCase):
    def setUp(self):
        self.fs = 100.0
        self.cleaning = CleaningConfig(
            block_size_seconds=1.0,
            reference_window_seconds=20.0,
            reference_quantile=0.90,
            refinement_iterations=3,
            cavity_lower_fraction=0.90,
            eta0_lower_fraction=0.80,
        )
        self.calibration = CalibrationConfig(
            cavity_transmission_w_per_v=0.020,
            responsivity_a_per_w=0.7,
            transimpedance_gain_v_per_a=1.0e6,
            second_harmonic_pd_response=0.7,
        )

    def _chunk(self, *, nominal_eta=0.022):
        n = 12000
        t = np.arange(n) / self.fs

        # Slow drifts in both observables should track their local q90 references.
        vtrans = 0.50 - 0.04 * (t / t[-1])
        eta = nominal_eta - 0.001 * (t / t[-1])

        cavity_glitch = slice(3500, 3600)  # exactly one 1 s reduced block
        eta_glitch = slice(8000, 8100)     # exactly one 1 s reduced block
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

    def test_rejects_block_glitches_but_preserves_slow_drift(self):
        chunk, cavity_glitch, eta_glitch = self._chunk()
        result = clean_chunk(chunk, self.cleaning, self.calibration)

        self.assertFalse(np.any(result.cavity_good[cavity_glitch]))
        self.assertTrue(np.all(result.cavity_good[eta_glitch]))
        self.assertFalse(np.any(result.eta_good[eta_glitch]))

        safe = np.ones(chunk.n_samples, dtype=bool)
        safe[cavity_glitch] = False
        safe[eta_glitch] = False
        self.assertGreater(np.mean(result.final_good[safe]), 0.99)

    def test_cavity_decision_uses_block_minimum_not_block_average(self):
        chunk, _, _ = self._chunk()
        v = chunk.cavity_transmission.copy()

        # A single deep sample is far too short to lower the 1 s block average
        # appreciably, but the old validated cavity rule deliberately uses the
        # block minimum for the lock decision.
        target = 5000
        v[target] *= 0.20
        changed = APEDXDChunk(
            source_file=chunk.source_file,
            start_index=chunk.start_index,
            stop_index=chunk.stop_index,
            sample_rate_hz=chunk.sample_rate_hz,
            x_1f=chunk.x_1f,
            y_1f=chunk.y_1f,
            x_2f=chunk.x_2f,
            y_2f=chunk.y_2f,
            cavity_transmission=v,
            first_timestamp_s=chunk.first_timestamp_s,
            last_timestamp_s=chunk.last_timestamp_s,
        )
        result = clean_chunk(changed, self.cleaning, self.calibration)
        block = slice(5000, 5100)
        self.assertFalse(np.any(result.cavity_good[block]))

    def test_eta0_is_relative_not_fixed_0p019(self):
        # A stable eta0 near 0.010 is below the old fixed 0.019 threshold but
        # should pass the new researcher-approved relative q90 rule.
        chunk, cavity_glitch, eta_glitch = self._chunk(nominal_eta=0.010)
        result = clean_chunk(chunk, self.cleaning, self.calibration)
        safe = np.ones(chunk.n_samples, dtype=bool)
        safe[cavity_glitch] = False
        safe[eta_glitch] = False
        self.assertGreater(np.mean(result.eta_good[safe]), 0.99)

    def test_eta0_decision_uses_sample_level_block_minimum(self):
        chunk, _, _ = self._chunk()
        x2 = chunk.x_2f.copy()

        # One deep eta0 sample should reject its whole reduced block even though
        # the block-mean eta0 remains close to nominal.
        target = 5000
        x2[target] *= 0.04
        changed = APEDXDChunk(
            source_file=chunk.source_file,
            start_index=chunk.start_index,
            stop_index=chunk.stop_index,
            sample_rate_hz=chunk.sample_rate_hz,
            x_1f=chunk.x_1f,
            y_1f=chunk.y_1f,
            x_2f=x2,
            y_2f=chunk.y_2f,
            cavity_transmission=chunk.cavity_transmission,
            first_timestamp_s=chunk.first_timestamp_s,
            last_timestamp_s=chunk.last_timestamp_s,
        )
        result = clean_chunk(changed, self.cleaning, self.calibration)
        block = slice(5000, 5100)
        self.assertFalse(np.any(result.eta_good[block]))

    def test_sustained_eta0_collapse_cannot_redefine_reference(self):
        chunk, _, _ = self._chunk()
        x2 = chunk.x_2f.copy()

        # 25 s is longer than the 20 s protected-history window used in this
        # synthetic test. The failed state must remain rejected throughout,
        # rather than becoming its own local reference.
        collapse = slice(5000, 7500)
        x2[collapse] *= 0.01
        changed = APEDXDChunk(
            source_file=chunk.source_file,
            start_index=chunk.start_index,
            stop_index=chunk.stop_index,
            sample_rate_hz=chunk.sample_rate_hz,
            x_1f=chunk.x_1f,
            y_1f=chunk.y_1f,
            x_2f=x2,
            y_2f=chunk.y_2f,
            cavity_transmission=chunk.cavity_transmission,
            first_timestamp_s=chunk.first_timestamp_s,
            last_timestamp_s=chunk.last_timestamp_s,
        )
        result = clean_chunk(changed, self.cleaning, self.calibration)

        self.assertTrue(np.all(result.cavity_good[collapse]))
        self.assertFalse(np.any(result.eta_good[collapse]))
        self.assertTrue(np.all(result.eta_good[7500:7900]))
        self.assertGreater(np.nanmedian(result.eta0_reference[collapse]), 0.015)

    def test_cavity_only_plan_matches_full_plan_cavity_stage(self):
        whole, _, _ = self._chunk()

        def factory():
            yield whole

        cavity_plan = build_cavity_reference_plan(factory, self.cleaning)
        full_plan = build_cleaning_reference_plan(
            factory, self.cleaning, self.calibration
        )

        np.testing.assert_allclose(
            cavity_plan.cavity_reference_blocks,
            full_plan.cavity_reference_blocks,
            equal_nan=True,
        )
        np.testing.assert_array_equal(
            cavity_plan.cavity_good_blocks,
            full_plan.cavity_good_blocks,
        )

    def test_full_record_reference_plan_crosses_chunk_boundaries(self):
        whole, cavity_glitch, eta_glitch = self._chunk()
        chunk_samples = 2000

        def factory():
            for start in range(0, whole.n_samples, chunk_samples):
                stop = min(whole.n_samples, start + chunk_samples)
                yield APEDXDChunk(
                    source_file=whole.source_file,
                    start_index=start,
                    stop_index=stop,
                    sample_rate_hz=whole.sample_rate_hz,
                    x_1f=whole.x_1f[start:stop],
                    y_1f=whole.y_1f[start:stop],
                    x_2f=whole.x_2f[start:stop],
                    y_2f=whole.y_2f[start:stop],
                    cavity_transmission=whole.cavity_transmission[start:stop],
                    first_timestamp_s=start / whole.sample_rate_hz,
                    last_timestamp_s=(stop - 1) / whole.sample_rate_hz,
                )

        plan = build_cleaning_reference_plan(
            factory, self.cleaning, self.calibration
        )
        masks = []
        for chunk in factory():
            result = clean_chunk_with_plan(
                chunk, plan, self.cleaning, self.calibration
            )
            masks.append(result.final_good)
        final_good = np.concatenate(masks)

        self.assertFalse(np.any(final_good[cavity_glitch]))
        self.assertFalse(np.any(final_good[eta_glitch]))

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
