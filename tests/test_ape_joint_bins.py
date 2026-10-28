"""Bounded mathematical checks; no native C, real data or generated caches."""
import unittest
import numpy as np
from scipy.stats import skewnorm
from tools.validate_ape_joint_bins import (
    draw_noise, kaiser_correlation, loglikelihood_grid, make_summary,
)


class JointBinBenchmarkTests(unittest.TestCase):
    def test_window_covariance_and_shared_noise(self):
        f = 5000.0 + np.arange(7)*.005
        gamma, duration, alpha = kaiser_correlation(f, 238.13)
        self.assertAlmostEqual(duration, 200., places=5)
        self.assertAlmostEqual(alpha, 9.549175266675151, places=10)
        self.assertAlmostEqual(abs(gamma[0, 1])**2, .849375, delta=3e-4)
        z = draw_noise(np.random.default_rng(919), 2000, gamma)
        self.assertEqual(z.shape, (2000, 2, 7))
        self.assertGreater(np.corrcoef(abs(z[:, 0, 0])**2, abs(z[:, 0, 1])**2)[0, 1], .7)

    def test_vectorization_exact_original_product_formula(self):
        power = np.array([[[2., 3., 4.], [5., 6., 7.]]])
        bg = np.ones((2, 3))
        template = np.array([[.1, .2, .3], [.3, .2, .1]])
        mus = np.array([0., 1., 2.])
        pars = (1.2, -.5, 1.1)
        expected = []
        for mu in mus:
            residual = np.log(power[0]) - np.log(bg + mu*template)
            expected.append(float(np.sum(skewnorm.logpdf(
                residual, a=pars[0], loc=pars[1], scale=pars[2]))))
        actual = loglikelihood_grid(power, bg, template, pars, mus)
        np.testing.assert_allclose(actual[0], expected, atol=1e-13)

    def test_no_asymptotic_detection_or_coverage_claim(self):
        grid = np.linspace(0, 4, 9)
        ll = np.stack([-np.square(grid - a) for a in [0., 1., 2., 1.]])
        got = make_summary(ll, ll, grid, 1.)
        self.assertIn("nominal_95pct_profile_coverage", got["injection"])
        self.assertEqual(got["null"]["n"], 4)


if __name__ == "__main__":
    unittest.main()
