"""Synthetic-only tests: APE-specific normalization boundary, vendor source untouched."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ape_signal_response import APEResponseError, load_peak_norm, vendor_likelihood_inputs


class APESignalResponseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "response.npz"
        self.freq = np.array([100.0, 100.1, 100.2])
        self.x = np.array([1e-12, 2e-12, 1e-12])
        self.y = np.array([3e-12, 2e-12, 3e-12])
        self.meta = {
            "mu_definition": "synthetic test scalar amplitude (NOT an axion coupling)",
            "response_source": "synthetic, no real detector calibration",
            "calibration_version": "synthetic-v0",
            "spectral_shape_convention": "unit sum, applied downstream by vendor likelihood",
            "response_units": "rad^2/Hz per unit mu",
        }
        self.save()

    def save(self, **replacements):
        fields = dict(frequency_hz=self.freq, peak_norm_x=self.x, peak_norm_y=self.y,
                      metadata_json=json.dumps(self.meta))
        fields.update(replacements)
        np.savez(self.path, **fields)

    def test_loads_exact_grid_and_separate_power_responses(self):
        norms, meta = load_peak_norm(self.path, self.freq)
        np.testing.assert_array_equal(norms, np.stack((self.x, self.y)))
        self.assertEqual(meta, self.meta)

    def test_rejects_frequency_interpolation(self):
        with self.assertRaisesRegex(APEResponseError, "frequency grid"):
            load_peak_norm(self.path, self.freq + np.array([0, 0.001, 0]))

    def test_rejects_negative_response(self):
        self.save(peak_norm_x=np.array([1e-12, -1e-12, 1e-12]))
        with self.assertRaisesRegex(APEResponseError, "nonnegative"):
            load_peak_norm(self.path, self.freq)

    def test_rejects_missing_physical_definition(self):
        self.meta.pop("mu_definition")
        self.save()
        with self.assertRaisesRegex(APEResponseError, "mu_definition"):
            load_peak_norm(self.path, self.freq)

    def test_passes_same_arrays_to_original_likelihood_interface(self):
        response, _ = load_peak_norm(self.path, self.freq)
        y = np.log(np.array([[1e-10, 1.1e-10, .9e-10], [2e-10, 1.9e-10, 2.1e-10]]))
        b = y.copy()
        pars = np.array([[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
        shape = np.array([0.2, 0.6, 0.2])
        args = vendor_likelihood_inputs(y, b, pars, response, shape)
        for actual, expected in zip(args, (y, b, response, shape, pars)):
            np.testing.assert_array_equal(actual, expected)

    def test_rejects_double_applied_or_unnormalized_shape(self):
        with self.assertRaisesRegex(APEResponseError, "normalized peak shape"):
            vendor_likelihood_inputs(np.zeros((2, 3)), np.zeros((2, 3)),
                                     np.ones((2, 3)), np.ones((2, 3)), np.ones(3))


if __name__ == "__main__":
    unittest.main()
