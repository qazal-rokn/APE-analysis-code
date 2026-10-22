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
    def log_likelihood(cls, mu, y, b, peak_norm, peak_shape, params):
        cls.called = True
        assert y.shape == b.shape == (1, 5)
        assert peak_shape.shape == (5,)
        assert np.isscalar(peak_norm) and peak_norm > 0
        assert params.shape == (1, 3)
        return float(-np.sum(
            (y - np.log(np.exp(b) + peak_norm*peak_shape[None, :]*mu))**2))


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
            actual = adapter.log_likelihood(0.3, data, bg, beta, pars, "/unused")
            direct = DummyStats.log_likelihood(
                0.3, np.log(data)[None, :], np.log(bg)[None, :],
                1.0, beta, pars[None, :])
            self.assertAlmostEqual(actual, direct)
            actual_scaled = adapter.log_likelihood(
                0.3, data, bg, beta, pars, "/unused", weights=np.array([0.7]))
            direct_scaled = DummyStats.log_likelihood(
                0.3, np.log(data)[None, :], np.log(bg)[None, :],
                0.7, beta, pars[None, :])
            self.assertAlmostEqual(actual_scaled, direct_scaled)
        self.assertTrue(DummyStats.called)
        with self.assertRaises(ValueError):
            adapter.log_likelihood(0, data, bg, np.array([-1]*5), pars, "/unused")


if __name__ == "__main__":
    unittest.main()
