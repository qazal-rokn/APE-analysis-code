"""Focused opt-in multi-bin template integration tests. No C rerun."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from ape_c_lpsd.analysis.multibin_diagnostic import evaluate_multibin_windows


class MultiBinTemplateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "physical.npz"
        self.frequency = np.array([4999.990, 4999.995, 5000., 5000.005, 5000.010])
        self.W = np.array([[[[.25, .15], [.75, .6], [1., .92]]]*2])
        self.W = np.broadcast_to(self.W, (1, 2, 3, 2)).copy()
        self.make_input()
        self.config = {
            "multibin_template_file": str(self.path),
            "halo_source_description": "Externally provided two-line comparison only; not a physical halo",
            "window_response_source_description": "User-provided demonstrative power-response matrix",
            "coupling_grid_gev_inv": [0., 1e-10, 2e-10],
            "n_effective_trials": None,
            "physics": {
                "quadrature_power_fraction": [.5, .5],
                "reference_coupling_gev_inv": 1e-10,
                "finesse_s": 1510., "tau_s": 2.5e-9, "tau1_s": 1e-8,
                "rho_dm_gev4": 3.1e-42, "unit_conversion": 1.519e24,
            },
        }
        self.fits = [
            {"background_psd": np.ones(5)*1e-19,
             "skewnormal_loc_scale_alpha": [0., 1., 0.]}
            for _ in range(2)
        ]

    def make_input(self, **overrides):
        arrays = {
            "frequency_support_hz": np.array([4999.999, 5000.001]),
            "normalized_halo_weights": np.array([.5, .5]),
            "observed_frequency_windows_hz": self.frequency[1:4][None,:],
            "window_power_response": self.W,
        }
        arrays.update(overrides)
        np.savez(self.path, **arrays)

    def evaluate(self):
        return evaluate_multibin_windows(
            frequency_hz=self.frequency,
            observed_ps=np.ones((2, 5))*1.1e-19,
            fitted_backgrounds=self.fits,
            candidate_start_indices=[1],
            config=self.config,
            external_vendor_root="/unused",
        )

    def test_exact_forward_template_and_original_vendor_calls(self):
        with patch("ape_c_lpsd.analysis.ape_likelihood.log_likelihood",
                   side_effect=lambda mu, *args: -float(mu)) as vendor:
            result = self.evaluate()
        self.assertEqual(vendor.call_count, 3)
        self.assertEqual(result["status"],
                         "provisional_multibin_forward_likelihood_diagnostic_only")
        self.assertIsNone(result["expected_upper_limit"])
        self.assertIsNone(result["candidate_significance"])
        self.assertIsNone(result["additional_normalization_factor"])
        np.testing.assert_allclose(result["windows"][0]["log_likelihood_by_coupling"],
                                   [0., -1., -4.])

    def test_reject_wrong_c_grid(self):
        self.make_input(observed_frequency_windows_hz=np.array([[1., 2., 3.]]))
        with self.assertRaisesRegex(ValueError, "match actual native C"):
            self.evaluate()

    def test_reject_bad_halo_normalization(self):
        self.make_input(normalized_halo_weights=np.array([.8, .8]))
        with self.assertRaisesRegex(ValueError, "halo probability"):
            self.evaluate()

    def test_reject_missing_source(self):
        self.config["halo_source_description"] = ""
        with self.assertRaisesRegex(ValueError, "provenance"):
            self.evaluate()

    def test_reject_global_trials(self):
        self.config["n_effective_trials"] = 10
        with self.assertRaisesRegex(ValueError, "global trials"):
            self.evaluate()

    def test_reject_incorrect_window_shape(self):
        self.make_input(window_power_response=np.ones((1, 3, 3, 2)))
        with self.assertRaisesRegex(ValueError, "Window power response"):
            self.evaluate()

    def test_reject_uncontrolled_coupling_grid(self):
        self.config["coupling_grid_gev_inv"] = [1e-10, 2e-10]
        with self.assertRaisesRegex(ValueError, "starting at zero"):
            self.evaluate()


if __name__ == "__main__":
    unittest.main()
