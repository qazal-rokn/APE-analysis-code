"""Direct vendor expected-UL parity checks; no C or measured data."""
import unittest
import numpy as np
from ape_c_lpsd.analysis.ape_expected_upper_limit import expected_upper_limit_ape


class Tests(unittest.TestCase):
    def test_vendor_known_result_and_coupling_scaling(self):
        ref = 1e-10
        result = expected_upper_limit_ape(
            sigma_mu=1.0, g_ref_gev_inv=ref, external_root="/source")
        self.assertAlmostEqual(result["sqrt_mu_ul"], 1.28251847040761, places=9)
        self.assertAlmostEqual(result["sigma_sqrt_mu"], 0.38985793307217637, places=9)
        self.assertAlmostEqual(result["g_ul_gev_inv"] / ref, result["sqrt_mu_ul"])
        self.assertTrue("not_ape_coverage_validated" in result["status"])

    def test_requires_valid_inputs(self):
        for sigma, ref in ((0, 1e-10), (-1, 1e-10), (1, 0), (np.nan, 1e-10)):
            with self.subTest(sigma=sigma, ref=ref), self.assertRaises(ValueError):
                expected_upper_limit_ape(
                    sigma_mu=sigma, g_ref_gev_inv=ref, external_root="/source")


if __name__ == "__main__":
    unittest.main()
