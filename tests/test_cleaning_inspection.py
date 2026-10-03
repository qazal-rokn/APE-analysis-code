"""Unit tests for streaming cleaning-interval diagnostics."""
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "tools" / "inspect_cleaning.py"
SPEC = importlib.util.spec_from_file_location("inspect_cleaning", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
IntervalAccumulator = MODULE._IntervalAccumulator


class CleaningInspectionIntervalTests(unittest.TestCase):
    def test_interval_crossing_chunk_boundary_is_one_run(self):
        acc = IntervalAccumulator(target_value=False)
        acc.feed(np.array([True, True, False, False]), 0)
        acc.feed(np.array([False, False, True, True]), 4)
        summary = acc.summary(sample_rate_hz=2.0)
        self.assertEqual(summary["count"], 1)
        self.assertEqual(summary["total_duration_s"], 2.0)
        self.assertEqual(summary["longest_start_index"], 2)
        self.assertEqual(summary["longest_stop_index_exclusive"], 6)

    def test_good_segments_and_quantiles(self):
        acc = IntervalAccumulator(target_value=True)
        acc.feed(np.array([True, True, False, True]), 0)
        acc.feed(np.array([True, True, False, True]), 4)
        summary = acc.summary(sample_rate_hz=1.0)
        self.assertEqual(summary["count"], 3)
        self.assertEqual(summary["max_duration_s"], 3.0)
        self.assertEqual(summary["longest_start_index"], 3)
        self.assertEqual(summary["longest_stop_index_exclusive"], 6)


if __name__ == "__main__":
    unittest.main()
