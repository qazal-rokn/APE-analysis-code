"""Step 5 deterministic APE-QWP template; no real data, no vendor edits."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ape_post_lpsd.ape_axion_signal_psd import (
    APEModelError, build_response_psd, inject_expected_psd,
    save_provisional_response,
)
from ape_post_lpsd.physics.qwp_signal import (
    qwp_amplitude_per_coupling, qwp_signal_power_template,
)
from ape_post_lpsd.physics.axion_halo import maxwell_halo_bin_weights


class APEStep5Tests(unittest.TestCase):
    def setUp(self):
        self.native = np.array([4999.995, 5000.0, 5000.005, 5000.010])
        self.halo = np.array([5000.0, 5000.005])
        self.weights = np.array([0.25, 0.75])
        self.ps_window = np.array([[[1, 0], [.25, .75]],
                                   [[1, 0], [.5, .5]]], dtype=float)
        self.conversion = np.array([7.0, 11.0])
        self.kernel = self.ps_window * self.conversion[None, :, None]
        self.physics = dict(finesse_s=1510.0, tau_s=2.5e-9, tau1_s=1e-8,
                            rho_dm_gev4=3.1e-42, unit_conversion=1.519e24)
        self.args = dict(native_frequency_hz=self.native, output_bin_indices=[1, 2],
                         halo_frequency_hz=self.halo,
                         normalized_halo_weights=self.weights,
                         native_c_psd_kernel_per_variance_hz=self.kernel,
                         quadrature_power_fraction=[0.8, 0.2],
                         reference_coupling_gev_inv=1e-10, **self.physics)

    def test_matches_existing_ape_peak_power_model_with_explicit_psd_normalization(self):
        computed = build_response_psd(**self.args)
        original_ps = qwp_signal_power_template(
            signal_frequency_hz=self.halo, normalized_halo_weights=self.weights,
            window_power_response=self.ps_window,
            quadrature_power_fraction=[0.8, 0.2],
            reference_coupling_gev_inv=1e-10, **self.physics)
        np.testing.assert_allclose(computed[:, [1, 2]],
                                   original_ps * self.conversion[None, :], rtol=1e-13)
        np.testing.assert_array_equal(computed[:, [0, 3]], np.zeros((2, 2)))

    def test_coupling_squared_deterministic_injection(self):
        signal = build_response_psd(**self.args)
        background = np.ones_like(signal) * 1e-25
        np.testing.assert_array_equal(inject_expected_psd(background, signal, mu=0), background)
        tested = inject_expected_psd(background, signal, mu=4)
        np.testing.assert_allclose(tested - background, 4 * signal, rtol=1e-14, atol=0)

    def test_zero_and_negative_kernel_refused(self):
        args = dict(self.args)
        args['native_c_psd_kernel_per_variance_hz'] = np.zeros_like(self.kernel)
        with self.assertRaisesRegex(APEModelError, 'zero'):
            build_response_psd(**args)
        args['native_c_psd_kernel_per_variance_hz'] = -np.ones_like(self.kernel)
        with self.assertRaisesRegex(APEModelError, 'kernel'):
            build_response_psd(**args)

    def test_wrong_grid_or_halo_weights_refused(self):
        args = dict(self.args)
        args['native_frequency_hz'] = self.native[::-1]
        with self.assertRaisesRegex(APEModelError, 'Native C frequency'):
            build_response_psd(**args)
        args = dict(self.args)
        args['normalized_halo_weights'] = np.array([0.4, 0.4])
        with self.assertRaisesRegex(APEModelError, 'sum to one'):
            build_response_psd(**args)

    def test_separate_quadratures_no_statistical_independence_assumed(self):
        x, y = build_response_psd(**self.args)
        self.assertGreater(np.max(x), 0)
        self.assertGreater(np.max(y), 0)
        self.assertFalse(np.array_equal(x, y))

    def test_npz_is_step_three_compatible_and_explicitly_provisional(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'provisional.npz'
            meta = save_provisional_response(
                out, **self.args, halo_source='synthetic halo fixture',
                native_c_kernel_source='synthetic explicitly scaled PS->PSD fixture',
                calibration_version='synthetic-only')
            with np.load(out, allow_pickle=False) as z:
                self.assertEqual(set(z.files), {'frequency_hz', 'peak_norm_x',
                                                'peak_norm_y', 'metadata_json'})
                np.testing.assert_array_equal(z['frequency_hz'], self.native)
                self.assertEqual(json.loads(str(z['metadata_json'].item())), meta)
            self.assertEqual(meta['response_units'], 'rad^2/Hz per unit mu')
            self.assertEqual(meta['optical_response_status'], 'PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS')
            with self.assertRaisesRegex(APEModelError, 'overwrites'):
                save_provisional_response(out, **self.args, halo_source='synthetic',
                    native_c_kernel_source='synthetic', calibration_version='synthetic')

    def test_original_unboosted_halo_helper_not_authorized_as_full_shm(self):
        # Existing module provides only a *benchmark* no lab/escape effects.
        with self.assertRaisesRegex(ValueError, 'misses too much'):
            maxwell_halo_bin_weights(frequency_edges_hz=[5000, 5000.00001],
                rest_frequency_hz=5000, velocity_scale_m_s=229000.,
                light_speed_m_s=299792458., min_captured_probability=0.99)

    def test_wrong_covariance_or_stochastic_injection_not_added(self):
        with self.assertRaisesRegex(APEModelError, 'nonnegative'):
            inject_expected_psd(np.ones((2, 3)), np.ones((2, 3)), mu=-1)


if __name__ == '__main__':
    unittest.main()
