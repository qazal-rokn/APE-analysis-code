"""Focused candidate tests, no real C execution."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from ape_c_lpsd.analysis.ape_candidates import select_candidates
from ape_c_lpsd.analysis.ape_pipeline_controller import read_candidate_ps_window


class Tests(unittest.TestCase):
    def test_terminal_cluster_is_preserved(self):
        result=select_candidates([1,2,3,4],[0]*4,[0,20,20,30],
                                 [0,1,1,1],local_sigma_threshold=5)
        self.assertEqual(len(result["candidates"]),1)
        self.assertEqual(result["candidates"][0]["frequency_hz"],4)
        self.assertEqual(result["candidates"][0]["cluster_end_index_exclusive"],4)

    def test_native_ps_alignment_rejects_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            a,b=Path(d)/"x.ps",Path(d)/"y.ps"
            a.write_text("1 2\n2 3\n")
            b.write_text("1 4\n2 5\n")
            f,p=read_candidate_ps_window({"X":a,"Y":b},start_index=0,bin_count=2)
            self.assertEqual(p.shape,(2,2))
            b.write_text("1 4\n2.1 5\n")
            with self.assertRaises(ValueError):
                read_candidate_ps_window({"X":a,"Y":b},start_index=0,bin_count=2)


if __name__=="__main__":
    unittest.main()
