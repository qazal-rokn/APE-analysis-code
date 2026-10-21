"""No-vendor modification tests for thin external statistics integration."""
import unittest
from unittest.mock import patch
import numpy as np
from ape_c_lpsd.analysis import vendor_statistics as adapter


class DummyFit:
    @staticmethod
    def process_iteration(args):
        _, x, y, config = args
        assert len(x) == len(y) and config["buffer"] == 50
        return 0, np.array([x[0], x[-1], y[0], y[-1]]), [0.0, 0.5, 2.0], 1.5

    @staticmethod
    def bkg_model(x, knots):
        return np.full_like(x, -4.0)


class DummyStats:
    called = False

    @classmethod
    def log_likelihood(cls, mu, y, b, beta, weights, params):
        cls.called = True
        assert y.shape == b.shape == beta.shape == (1, 5)
        assert params.shape == (1, 3)
        assert weights.shape == (1,)
        return float(-np.sum((y-b-mu*beta)**2))


class IntegrationTests(unittest.TestCase):
    def test_background_delegates_without_altering_input(self):
        f = np.linspace(100, 101, 1100)
        p = np.exp(np.full(1100, -4.0))
        frozen = p.copy()
        with patch.object(adapter, "load_external", return_value=(DummyFit, DummyStats, {})):
            out = adapter.fit_window(f, p, "/unused")
        np.testing.assert_array_equal(p, frozen)
        self.assertEqual(out["vendor_reduced_chi2"], 1.5)
        self.assertEqual(out["skewnormal_loc_scale_alpha"].tolist(), [0.0, 0.5, 2.0])
        self.assertTrue(out["vendor_quality_pass"])
        self.assertAlmostEqual(float(np.median(out["background_psd"])), np.exp(-4.0))

    def test_likelihood_delegates_without_default_signal_model(self):
        data = np.full(5, 2.0)
        bg = np.full(5, 1.0)
        beta = np.full(5, 0.5)
        pars = np.array([0.0, 1.0, 0.0])
        DummyStats.called = False
        with patch.object(adapter, "load_external", return_value=(DummyFit, DummyStats, {})):
            self.assertTrue(np.isfinite(adapter.log_likelihood(
                0.3, data, bg, beta, pars, "/unused")))
        self.assertTrue(DummyStats.called)
        with self.assertRaises(ValueError):
            adapter.log_likelihood(0, data, bg, np.array([-1]*5), pars, "/unused")


if __name__ == "__main__":
    unittest.main()
