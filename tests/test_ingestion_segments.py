"""Synthetic-only ingestion tests; intentionally not run in this implementation stage."""
import unittest
from pathlib import Path
from ape_c_lpsd.ingestion.models import ValidRun
from ape_c_lpsd.ingestion.segments import build_global_segments, choose_longest_global_segment

class SegmentTests(unittest.TestCase):
    def test_verified_file_endpoint_can_join(self):
        a=ValidRun(0,4,10,1.,False,Path('a')); b=ValidRun(1,0,8,1.,True,Path('b'))
        chosen=choose_longest_global_segment(build_global_segments((a,b),source_sample_counts={0:10,1:8}))
        self.assertEqual(chosen.sample_count,14)
    def test_unverified_boundary_never_joins(self):
        a=ValidRun(0,4,10,1.,False,Path('a')); b=ValidRun(1,0,8,1.,False,Path('b'))
        self.assertEqual(len(build_global_segments((a,b),source_sample_counts={0:10,1:8})),2)
