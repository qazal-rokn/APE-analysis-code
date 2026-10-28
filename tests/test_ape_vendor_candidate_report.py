"""Local q0 candidate bridge tests; no native C or APE data."""
import unittest
from unittest.mock import patch
import numpy as np
from ape_c_lpsd.analysis import ape_vendor_candidate_report as bridge


class Tests(unittest.TestCase):
    def test_terminal_cluster_and_local_q0_exact(self):
        # This fixture tests q0 and grouping, not any physical axion distribution.
        fitted = [
            dict(log_likelihood_null=0., log_likelihood_max=0., mu_hat=0.,
                 expected_upper_limit={'sqrt_mu_ul': 1.0}),
            dict(log_likelihood_null=0., log_likelihood_max=3., mu_hat=1.,
                 expected_upper_limit={'sqrt_mu_ul': 1.0}),
            dict(log_likelihood_null=0., log_likelihood_max=14., mu_hat=2.,
                 expected_upper_limit={'sqrt_mu_ul': 1.0}),
        ]
        with patch.object(bridge, 'fit_vendor_mle_and_expected_ul', side_effect=fitted) as call:
            result = bridge.fit_and_report_single_bin_candidates(
                candidate_frequencies_hz=[4999., 5000., 5001.],
                prepared_windows=[{}, {}, {}], external_root='/source',
                min_log10mu=-3, max_log10mu=2, n_seed_lkl=120,
                local_sigma_threshold=3.)
        self.assertEqual(call.call_count, 3)
        self.assertEqual(result['selection']['flagged_bins'], 1)
        self.assertEqual(len(result['selection']['candidates']), 1)
        self.assertEqual(result['selection']['candidates'][0]['frequency_hz'], 5001.)
        self.assertAlmostEqual(result['selection']['candidates'][0]['q0'], 28.)
        self.assertEqual(result['selection']['candidates'][0]['cluster_end_index_exclusive'], 3)
        self.assertFalse(result['global_significance_calibrated'])

    def test_real_original_vendor_fit_to_candidate_report(self):
        # Identical deterministic vendor-model fixture previously checked by user.
        rng = np.random.default_rng(20261028)
        n = 64
        prepared = dict(
            Y=np.log(1.5e-16) + rng.normal(0, .25, (n, 1)),
            bkg=np.log(np.full((n, 1), 1e-16)),
            peak_norm=np.full((n, 1), .5e-16),
            peak_shape=np.ones(1),
            model_args=np.tile([0., .25, 0.], (n, 1)),
            reference_coupling_gev_inv=1e-10,
        )
        result = bridge.fit_and_report_single_bin_candidates(
            candidate_frequencies_hz=[5000.], prepared_windows=[prepared],
            external_root='/source', min_log10mu=-3,
            max_log10mu=2, n_seed_lkl=120, local_sigma_threshold=5.)
        self.assertAlmostEqual(result['fits'][0]['mu_hat'], 1.0359235386660495, places=6)
        self.assertAlmostEqual(result['fits'][0]['q0'], 178.37717433895628, places=5)
        self.assertEqual(len(result['selection']['candidates']), 1)
        self.assertAlmostEqual(result['fits'][0]['expected_upper_limit']['sqrt_mu_ul'],
                               .3950310474302729, places=6)

    def test_bad_frequency_alignment_rejected(self):
        for f in ([5000., 5000.], [5001., 5000.], [-1., 5000.]):
            with self.subTest(f=f), self.assertRaises(ValueError):
                bridge.fit_and_report_single_bin_candidates(
                    candidate_frequencies_hz=f, prepared_windows=[{}, {}],
                    external_root='/source', min_log10mu=-3,
                    max_log10mu=2, n_seed_lkl=120,
                    local_sigma_threshold=3.)


if __name__ == '__main__':
    unittest.main()
