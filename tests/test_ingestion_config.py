"""Tests for explicit APE DXD channel and calibration configuration."""
from __future__ import annotations

from pathlib import Path
import tempfile
import textwrap
import unittest

from ape_c_lpsd.ingestion.config import (
    IngestionConfigError,
    load_dxd_ingestion_config,
)


BASE = """
input_file: /tmp/run.dxd
reference_project: /tmp/reference
output_directory: /tmp/output
channels:
  x_1f: pem1_x
  y_1f: pem1_Y
  x_2f: pem2_X
  y_2f: pem2_Y
  cavity_transmission: cavity_trans
cleaning:
  block_size_seconds: 0.05
  reference_window_seconds: 60.0
  cavity_lower_fraction: 0.90
  eta0_lower_fraction: 0.90
cleaning:
  block_size_seconds: 0.05
  reference_window_seconds: 60.0
  cavity_lower_fraction: 0.90
  eta0_lower_fraction: 0.90
calibration:
  cavity_transmission_w_per_v: 0.020
  responsivity_a_per_w: 0.7
  transimpedance_gain_v_per_a: 1000000.0
  second_harmonic_pd_response: 0.7
"""


def _load(text: str):
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "config.yaml"
        path.write_text(textwrap.dedent(text))
        return load_dxd_ingestion_config(path)


class DXDConfigTests(unittest.TestCase):
    def test_explicit_values_are_loaded_without_hidden_defaults(self):
        config = _load(BASE)
        self.assertEqual(config.channels.x_1f, "pem1_x")
        self.assertEqual(config.channels.y_1f, "pem1_Y")
        self.assertAlmostEqual(config.cleaning.block_size_seconds, 0.05)
        self.assertAlmostEqual(config.cleaning.reference_window_seconds, 60.0)
        self.assertAlmostEqual(config.cleaning.cavity_lower_fraction, 0.90)
        self.assertAlmostEqual(config.cleaning.eta0_lower_fraction, 0.90)
        self.assertAlmostEqual(config.cleaning.block_size_seconds, 0.05)
        self.assertAlmostEqual(config.cleaning.reference_window_seconds, 60.0)
        self.assertAlmostEqual(config.cleaning.cavity_lower_fraction, 0.90)
        self.assertAlmostEqual(config.cleaning.eta0_lower_fraction, 0.90)
        self.assertAlmostEqual(config.calibration.cavity_transmission_w_per_v, 0.020)
        self.assertAlmostEqual(config.calibration.second_harmonic_pd_response, 0.7)

    def test_missing_calibration_value_is_refused(self):
        text = BASE.replace("  second_harmonic_pd_response: 0.7\n", "")
        with self.assertRaisesRegex(
            IngestionConfigError, "second_harmonic_pd_response is required"
        ):
            _load(text)

    def test_duplicate_physical_channel_mapping_is_refused(self):
        text = BASE.replace("  y_1f: pem1_Y", "  y_1f: pem1_x")
        with self.assertRaisesRegex(IngestionConfigError, "five distinct"):
            _load(text)


if __name__ == "__main__":
    unittest.main()
