"""Synthetic tests for the non-numerical APE response bridge (vendor untouched)."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from ape_post_lpsd import build_ape_likelihood_response as bridge


class ResponseBridgeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.h5 = self.root / "APE_0_1_X.h5"
        self.kernel = self.root / "kernel.npz"
        self.cfg = self.root / "params.json"
        self.out = self.root / "response.npz"
        self.grid = np.array([100., 100.1, 100.2])
        with h5py.File(self.h5, "w") as h:
            h["frequency"] = self.grid
        self.indices = np.arange(3)
        self.meta = {
            "kernel_units": "Hz^-1 per unit time-domain rad^2 variance",
            "source_format": "native C method-2 D = PSD rad^2/Hz",
            "c_configuration_source": "synthetic test",
        }
        self.write_kernel()
        self.parameters = {
            "normalized_halo_weights": [1.0],
            "quadrature_power_fraction": [0.5, 0.5],
            "reference_coupling_gev_inv": 1.e-10,
            "finesse_s": 1000., "tau_s": 1.e-8, "tau1_s": 1.e-8,
            "rho_dm_gev4": 1.e-42, "unit_conversion": 1.,
            "halo_source": "synthetic-only", "calibration_version": "synthetic-only",
        }
        self.cfg.write_text(json.dumps(self.parameters))

    def write_kernel(self, **changes):
        data = dict(
            native_frequency_hz=self.grid, output_bin_indices=self.indices,
            halo_frequency_hz=np.array([100.]),
            native_c_psd_kernel_per_variance_hz=np.ones((2, len(self.indices), 1)),
            metadata_json=json.dumps(self.meta),
        )
        data.update(changes)
        np.savez(self.kernel, **data)

    def build(self, **overrides):
        args = dict(spectrum_hdf5=self.h5, kernel_npz=self.kernel,
                    physics_json=self.cfg, output_npz=self.out)
        args.update(overrides)
        return bridge.build_from_existing_kernel(**args)

    def test_passes_arrays_and_physics_unchanged_to_existing_builder(self):
        seen = {}
        def fake_builder(path, **kwargs):
            seen.update(kwargs)
            np.savez(path, frequency_hz=kwargs["native_frequency_hz"])
            return {"optical_response_status": "PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS"}
        with patch.object(bridge, "save_provisional_response", fake_builder):
            report = self.build()
        np.testing.assert_array_equal(seen["native_frequency_hz"], self.grid)
        np.testing.assert_array_equal(seen["native_c_psd_kernel_per_variance_hz"],
                                      np.ones((2, 3, 1)))
        self.assertEqual(seen["reference_coupling_gev_inv"], 1.e-10)
        self.assertEqual(report["native_bins"], 3)
        self.assertIn("PROVISIONAL", report["physics_status"])

    def test_reject_incomplete_full_grid(self):
        self.indices = np.array([0, 2])
        self.write_kernel()
        with self.assertRaisesRegex(bridge.ResponseBuildError, "every native bin"):
            self.build()

    def test_reject_grid_mismatch(self):
        self.write_kernel(native_frequency_hz=[100., 100.15, 100.2])
        with self.assertRaisesRegex(bridge.ResponseBuildError, "frequency grids differ"):
            self.build()

    def test_reject_wrong_units(self):
        self.meta["kernel_units"] = "rad^2"
        self.write_kernel()
        with self.assertRaisesRegex(bridge.ResponseBuildError, "D/PSD"):
            self.build()

    def test_reject_missing_physics(self):
        del self.parameters["rho_dm_gev4"]
        self.cfg.write_text(json.dumps(self.parameters))
        with self.assertRaisesRegex(bridge.ResponseBuildError, "exactly"):
            self.build()

    def test_refuse_overwriting(self):
        self.out.write_text("do-not-change")
        with self.assertRaisesRegex(bridge.ResponseBuildError, "no overwrites"):
            self.build()
        self.assertEqual(self.out.read_text(), "do-not-change")


if __name__ == "__main__":
    unittest.main()
