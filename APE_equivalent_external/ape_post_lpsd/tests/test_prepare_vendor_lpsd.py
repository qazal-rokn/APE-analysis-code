"""Synthetic contract tests for the separate post-LPSD handoff; vendor untouched."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_vendor_lpsd import HandoffError, prepare_pair


class HandoffTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.x = self.root / 'rho_x_combined.lpsd'
        self.y = self.root / 'rho_y_combined.lpsd'
        self.x.write_text('# original C PSD\n1 0.0001 20 2\n2 0.0002 20 4\n3 0.0003 30 5\n')
        self.y.write_text('1 0.0004 20 2\n2 0.0005 20 4\n3 0.0006 30 5\n')

    def run_handoff(self, **kw):
        return prepare_pair(self.x, self.y, self.root / 'product', fs=20., t0=100, t1=120,
                            expected_bins=3, chunk_size=2, **kw)

    def test_preserves_native_grid_and_independent_psd(self):
        summary = self.run_handoff()
        self.assertEqual(summary['bin_count_per_quadrature'], 3)
        for label, psd in [('X', [0.0001, 0.0002, 0.0003]),
                           ('Y', [0.0004, 0.0005, 0.0006])]:
            with h5py.File(summary['output'][label]) as f:
                np.testing.assert_array_equal(f['frequency'][:], [1., 2., 3.])
                np.testing.assert_array_equal(f['PSD'][:], psd)
                np.testing.assert_array_equal(f['logPSD'][:], np.log(psd))
                np.testing.assert_array_equal(f['N_avg'][:], [20., 20., 30.])
                np.testing.assert_array_equal(f['N_kernel'][:], [2., 4., 5.])
                self.assertEqual(f['logPSD'].attrs['Jdes'], 3)
                self.assertEqual(f['logPSD'].attrs['fs'], 20.)
                self.assertEqual(f.attrs['channel'], label)
        self.assertTrue((self.root/'product/APE_100_120_handoff.json').is_file())

    def test_unchanged_external_io_reader_accepts_native_arrays(self):
        # The repository checkout exposes the unchanged external package at
        # ../../scalar_dark_matter_lpsd relative to this test file.
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from scalar_dark_matter_lpsd.LPSDIO import LPSDOutput, LPSDDataGroup
        summary = self.run_handoff()
        x = LPSDOutput(summary['output']['X'])
        np.testing.assert_array_equal(x.freq, [1., 2., 3.])
        np.testing.assert_array_equal(x.logPSD, np.log([.0001, .0002, .0003]))
        group = LPSDDataGroup(summary['output']['X'])
        self.assertEqual(group.metadata['ifo'], ['X'])
        self.assertEqual(group.logPSD.shape, (1, 3))

    def test_rejects_xy_mismatch_and_removes_partial_products(self):
        self.y.write_text('1 .0004 2 3\n2.1 .0005 2 3\n3 .0006 2 3\n')
        with self.assertRaisesRegex(HandoffError, 'native frequency mismatch'):
            self.run_handoff()
        self.assertEqual(list((self.root / 'product').iterdir()), [])

    def test_rejects_invalid_psd(self):
        self.x.write_text('1 0 2 3\n2 1 2 3\n3 1 2 3\n')
        with self.assertRaisesRegex(HandoffError, 'invalid frequency'):
            self.run_handoff()
        self.assertEqual(list((self.root / 'product').iterdir()), [])

    def test_rejects_mismatched_count(self):
        self.y.write_text('1 .0004 2 3\n2 .0005 2 3\n')
        with self.assertRaisesRegex(HandoffError, 'count differs'):
            self.run_handoff()
        self.assertEqual(list((self.root / 'product').iterdir()), [])

    def test_refuses_overwrite(self):
        self.run_handoff()
        with self.assertRaisesRegex(HandoffError, 'overwrite'):
            self.run_handoff()

    def test_rejects_missing_columns(self):
        self.x.write_text('1 2 3\n2 3 4\n3 4 5\n')
        with self.assertRaisesRegex(HandoffError, 'four columns'):
            self.run_handoff()


if __name__ == '__main__':
    unittest.main()
