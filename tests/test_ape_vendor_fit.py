"""Direct vendor optimizer parity with APE bridge. No synthetic C or real data."""
import unittest
import numpy as np
from ape_c_lpsd.analysis.ape_vendor_fit import fit_vendor_mle_and_expected_ul
from scalar_dark_matter_lpsd.lvk.max_lkl_executable import calc_max_lkl


class Tests(unittest.TestCase):
    def test_optimizer_and_expected_ul_parity(self):
        rng = np.random.default_rng(20261028)
        n = 64
        bkg = np.log(np.full((n, 1), 1e-16))
        Y = np.log(1.5e-16) + rng.normal(0, 0.25, (n, 1))
        norm = np.full((n, 1), 0.5e-16)
        pars = np.tile([0.0, 0.25, 0.0], (n, 1))
        prepared = dict(Y=Y, bkg=bkg, peak_norm=norm, peak_shape=np.ones(1),
                        model_args=pars, reference_coupling_gev_inv=1e-10)
        direct = calc_max_lkl(Y, bkg, pars, norm, np.ones(1),
                              min_log10mu=-3, max_log10mu=2, n_seed_lkl=120)
        result = fit_vendor_mle_and_expected_ul(
            prepared=prepared, external_root='/source', min_log10mu=-3,
            max_log10mu=2, n_seed_lkl=120)
        for key, number in zip(('log_likelihood_null', 'log_likelihood_max',
                                 'mu_hat', 'sigma_mu'), direct):
            self.assertAlmostEqual(result[key], number, places=8)
        self.assertAlmostEqual(result['mu_hat'], 1.0359235386660495, places=6)
        self.assertAlmostEqual(result['sigma_mu'], 0.09487137692882767, places=6)
        self.assertAlmostEqual(result['q0'], 178.37717433895628, places=5)
        self.assertAlmostEqual(result['expected_upper_limit']['sqrt_mu_ul'],
                               0.3950310474302729, places=6)
        self.assertAlmostEqual(result['expected_upper_limit']['g_ul_gev_inv'],
                               3.950310474302729e-11, places=16)
        self.assertIn('not_ape_coverage_validated', result['status'])

    def test_multibin_pilot_refused(self):
        with self.assertRaises(ValueError):
            fit_vendor_mle_and_expected_ul(
                prepared=dict(Y=np.ones((1, 2)), bkg=np.ones((1, 2)),
                              peak_norm=np.ones((1, 2)), peak_shape=[0.5, 0.5],
                              model_args=[[0, 1, 0]], reference_coupling_gev_inv=1e-10),
                external_root='/source', min_log10mu=-3,
                max_log10mu=2, n_seed_lkl=120)


if __name__ == '__main__':
    unittest.main()
