"""Synthetic-only tests for the canonical five-channel DXD contract."""
from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np

from ape_c_lpsd.ingestion.dxd_reader import (
    DXDReaderError,
    _assemble_chunk,
    _validate_declared_rate,
    _validate_mapping,
)


NAMES = {
    "x_1f": "pem1_x",
    "y_1f": "pem1_Y",
    "x_2f": "pem2_X",
    "y_2f": "pem2_Y",
    "cavity_transmission": "cavity_trans",
}


class _FakeSession:
    def __init__(self, *, fs=20000.0, total=20, timestamp_offset_by_role=None):
        self.sample_rate_hz = fs
        self.source = Path("synthetic.dxd")
        self.channel_names = dict(NAMES)
        self.total = total
        self.timestamp_offset_by_role = timestamp_offset_by_role or {}

    def read_role_chunk(self, role, raw_start, count):
        idx = np.arange(raw_start, raw_start + count, dtype=np.float64)
        data = idx + {
            "x_1f": 0.0,
            "y_1f": 1.0,
            "x_2f": 2.0,
            "y_2f": 3.0,
            "cavity_transmission": 4.0,
        }[role]
        t = idx / self.sample_rate_hz
        t = t + float(self.timestamp_offset_by_role.get(role, 0.0))
        return data, t


class DXDChunkTests(unittest.TestCase):
    def test_mapping_requires_exact_five_distinct_channels(self):
        self.assertEqual(_validate_mapping(NAMES), list(NAMES.values()))
        duplicate = dict(NAMES)
        duplicate["y_1f"] = "pem1_x"
        with self.assertRaisesRegex(DXDReaderError, "five distinct"):
            _validate_mapping(duplicate)

    def test_declared_rate_is_validated_without_reestimating_it(self):
        t = np.arange(1000, dtype=np.float64) / 20000.0
        _validate_declared_rate(t, 20000.0)
        with self.assertRaisesRegex(DXDReaderError, "declared sample rate"):
            _validate_declared_rate(t, 19900.0)

    def test_chunk_contains_five_synchronized_arrays(self):
        session = _FakeSession()
        chunk, last = _assemble_chunk(session, 0, 10, None)
        self.assertEqual(chunk.raw_start, 0)
        self.assertEqual(chunk.raw_stop, 10)
        self.assertEqual(chunk.n_samples, 10)
        self.assertEqual(chunk.sample_rate_hz, 20000.0)
        self.assertAlmostEqual(last, 9 / 20000.0)
        np.testing.assert_allclose(chunk.y_1f, chunk.x_1f + 1.0)

    def test_cross_channel_timestamp_disagreement_is_refused(self):
        session = _FakeSession(timestamp_offset_by_role={"y_1f": 1e-4})
        with self.assertRaisesRegex(DXDReaderError, "not synchronized"):
            _assemble_chunk(session, 0, 10, None)

    def test_chunk_boundary_discontinuity_is_refused(self):
        session = _FakeSession()
        first, last = _assemble_chunk(session, 0, 10, None)
        self.assertEqual(first.raw_stop, 10)
        with self.assertRaisesRegex(DXDReaderError, "discontinuity"):
            _assemble_chunk(session, 11, 9, last)

    def test_adjacent_chunks_are_continuous(self):
        session = _FakeSession()
        first, last = _assemble_chunk(session, 0, 10, None)
        second, last2 = _assemble_chunk(session, 10, 10, last)
        self.assertEqual(first.raw_stop, second.raw_start)
        self.assertAlmostEqual(last2, 19 / 20000.0)


if __name__ == "__main__":
    unittest.main()
