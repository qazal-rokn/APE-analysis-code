"""Non-destructive dimensional and guard tests for explicit unit boundary."""
import unittest
import numpy as np
from ape_c_lpsd.analysis.spectral_units import (
    archived_psd_to_ps, amplitude_from_power, combined_xy_amplitude,
    coupling_from_signal_amplitude, coupling_from_fitted_power,
)


class SpectralUnitsTests(unittest.TestCase):
    def test_archived_psd_requires_external_bandwidth(self):
        np.testing.assert_allclose(
            archived_psd_to_ps([4., 9.], effective_bandwidth_hz=[.25, .5]),
            [1., 4.5])
        with self.assertRaises(TypeError):
            archived_psd_to_ps([4.])
        with self.assertRaises(ValueError):
            archived_psd_to_ps([4.], effective_bandwidth_hz=0)

    def test_xy_diagnostic_and_amplitude(self):
        np.testing.assert_allclose(combined_xy_amplitude([3.], [6.]), [3.])
        np.testing.assert_allclose(amplitude_from_power([4., 9.]), [2., 3.])
        with self.assertRaises(ValueError):
            amplitude_from_power([-1.])

    def test_amplitude_and_mu_coupling_agree(self):
        g_ref=1e-10
        mu=np.array([0., 1., 4.])
        np.testing.assert_allclose(
            coupling_from_fitted_power(mu, g_ref_gev_inv=g_ref),
            coupling_from_signal_amplitude(np.sqrt(mu),
                response_rad_per_gev_inv=1/g_ref))
        np.testing.assert_allclose(
            coupling_from_fitted_power(mu, g_ref_gev_inv=g_ref),
            [0, 1e-10, 2e-10])

    def test_invalid_quadratures(self):
        with self.assertRaises(ValueError):
            combined_xy_amplitude([1, 2], [1])
        with self.assertRaises(ValueError):
            coupling_from_signal_amplitude([1], response_rad_per_gev_inv=0)


if __name__ == "__main__":
    unittest.main()
