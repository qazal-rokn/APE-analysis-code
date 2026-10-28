"""Focused analytic checks for the paper-faithful injection; no native C run."""
import unittest

import numpy as np

from tools.run_paper_axion_like_injection import (
    F_REST, Q, SEED, sample_modes, window_power_kernel,
    baseband_envelope,
)


class PaperAxionInjectionTests(unittest.TestCase):
    def test_halo_positive_broadened_and_reproducible(self):
        a, phi = sample_modes(np.random.default_rng(SEED), 512)
        b, phi2 = sample_modes(np.random.default_rng(SEED), 512)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(phi, phi2)
        self.assertTrue(np.all(a >= F_REST))
        self.assertGreater(np.ptp(a), F_REST/Q)
        self.assertTrue(np.all(np.diff(a) > 0))

    def test_c_kernel_normalization_without_unit_sum(self):
        f = np.array([F_REST, F_REST*(1+1/Q)])
        W = window_power_kernel(f, np.array([F_REST]), q=Q, psll=238.13)
        self.assertAlmostEqual(float(W[0,0]), 1., places=12)
        self.assertGreater(float(W.sum()), 1.)  # overlapping window; DO NOT renormalize

    def test_single_mode_has_continuous_phase(self):
        f = np.array([F_REST+0.004])
        p = np.array([0.])
        t, envelope = baseband_envelope(f, p, np.array([2.]), 3.,
                                        f_rest=F_REST, duration=3.,
                                        envelope_rate=8.)
        np.testing.assert_allclose(abs(envelope), 6., atol=1e-10)
        np.testing.assert_allclose(np.angle(envelope[0]), 0., atol=1e-12)


if __name__ == "__main__":
    unittest.main()
