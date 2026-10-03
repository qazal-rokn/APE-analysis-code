"""Synthetic-only tests for the canonical five-channel DXD contract."""
from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np

from ape_c_lpsd.ingestion.dxd_reader import DXDReaderError, _assemble_ape_data


def _inputs(n: int = 1000, fs: float = 20000.0):
    t = np.arange(n, dtype=np.float64) / fs
    arrays = {
        "x_1f": np.arange(n, dtype=np.float64),
        "y_1f": np.arange(n, dtype=np.float64) + 1,
        "x_2f": np.arange(n, dtype=np.float64) + 2,
        "y_2f": np.arange(n, dtype=np.float64) + 3,
        "cavity_transmission": np.arange(n, dtype=np.float64) + 4,
    }
    timestamps = {key: t.copy() for key in arrays}
    units = {key: "V" for key in arrays}
    names = {
        "x_1f": "pem1_x",
        "y_1f": "pem1_Y",
        "x_2f": "pem2_X",
        "y_2f": "pem2_Y",
        "cavity_transmission": "cavity_trans",
    }
    return arrays, timestamps, units, names


class DXDAssemblyTests(unittest.TestCase):
    def test_uniform_synchronized_channels_keep_no_redundant_time_array(self):
        arrays, timestamps, units, names = _inputs()
        data = _assemble_ape_data(
            source_file=Path("synthetic.dxd"),
            arrays=arrays,
            timestamps=timestamps,
            units=units,
            channel_names=names,
            start_time=None,
            metadata={"synthetic": True},
            declared_sample_rate_hz=20000.0,
        )
        self.assertEqual(data.n_samples, 1000)
        self.assertTrue(data.regular_sampling)
        self.assertEqual(data.sample_rate_hz, 20000.0)
        self.assertIsNone(data.relative_timestamps_s)
        self.assertEqual(data.channel_names["x_1f"], "pem1_x")


    def test_declared_rate_disagreement_is_refused(self):
        arrays, timestamps, units, names = _inputs()
        with self.assertRaisesRegex(
            DXDReaderError, "inconsistent with the Dewesoft declared sample rate"
        ):
            _assemble_ape_data(
                source_file=Path("synthetic.dxd"),
                arrays=arrays,
                timestamps=timestamps,
                units=units,
                channel_names=names,
                start_time=None,
                metadata={},
                declared_sample_rate_hz=19900.0,
            )

    def test_length_disagreement_is_refused_not_trimmed(self):
        arrays, timestamps, units, names = _inputs()
        arrays["y_1f"] = arrays["y_1f"][:-1]
        timestamps["y_1f"] = timestamps["y_1f"][:-1]
        with self.assertRaisesRegex(DXDReaderError, "lengths disagree"):
            _assemble_ape_data(
                source_file=Path("synthetic.dxd"),
                arrays=arrays,
                timestamps=timestamps,
                units=units,
                channel_names=names,
                start_time=None,
                metadata={},
            )

    def test_timestamp_disagreement_is_refused_not_resampled(self):
        arrays, timestamps, units, names = _inputs()
        timestamps["y_1f"] = timestamps["y_1f"] + 1e-4
        with self.assertRaisesRegex(DXDReaderError, "not synchronized"):
            _assemble_ape_data(
                source_file=Path("synthetic.dxd"),
                arrays=arrays,
                timestamps=timestamps,
                units=units,
                channel_names=names,
                start_time=None,
                metadata={},
            )

    def test_irregular_timing_is_preserved_and_not_given_fake_sample_rate(self):
        arrays, timestamps, units, names = _inputs(n=10)
        irregular = timestamps["x_1f"].copy()
        irregular[5:] += 2e-5
        timestamps = {key: irregular.copy() for key in timestamps}
        data = _assemble_ape_data(
            source_file=Path("synthetic.dxd"),
            arrays=arrays,
            timestamps=timestamps,
            units=units,
            channel_names=names,
            start_time=None,
            metadata={},
        )
        self.assertFalse(data.regular_sampling)
        self.assertIsNone(data.sample_rate_hz)
        self.assertIsNotNone(data.relative_timestamps_s)
        np.testing.assert_allclose(
            data.relative_timestamps_s, irregular - irregular[0]
        )


if __name__ == "__main__":
    unittest.main()
