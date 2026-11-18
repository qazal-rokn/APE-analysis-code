"""Contract tests for APE adapter that invokes unchanged vendor LPSDOutput."""
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from ape_post_lpsd.prepare_vendor_lpsd_original import (
    VendorConversionError, convert_pair_original,
)


class OriginalVendorConverterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.x = self.root / 'x.lpsd'
        self.y = self.root / 'y.lpsd'

    def fill(self, freqs):
        self.x.write_text(''.join(
            f'{f:.17g} {(i + 1) * 1e-9:.17g} 10 3\n'
            for i, f in enumerate(freqs)))
        self.y.write_text(''.join(
            f'{f:.17g} {(i + 1) * 2e-9:.17g} 10 3\n'
            for i, f in enumerate(freqs)))

    def convert(self):
        return convert_pair_original(
            self.x, self.y, self.root / 'result',
            fs=20000.0, t0=100, t1=200, expected_bins=5)

    def test_original_vendor_log_grid_and_psd(self):
        frequencies = np.geomspace(100.0, 110.0, 5)
        self.fill(frequencies)
        summary = self.convert()
        for channel, factor in [('X', 1e-9), ('Y', 2e-9)]:
            with h5py.File(summary['outputs'][channel], 'r') as h:
                self.assertEqual(set(h.keys()), {'frequency', 'logPSD'})
                np.testing.assert_allclose(h['frequency'][:], frequencies, rtol=1e-13)
                np.testing.assert_allclose(
                    h['logPSD'][:], np.log(np.arange(1, 6) * factor))
                self.assertEqual(h['logPSD'].attrs['Jdes'], 5)
        self.assertEqual(summary['grid_comparison']['X']['bins_differing_above_1e-10_relative'], 0)

    def test_original_grid_reconstruction_is_not_overridden(self):
        self.fill([100, 101, 103, 106, 110])
        summary = self.convert()
        with h5py.File(summary['outputs']['X'], 'r') as h:
            np.testing.assert_allclose(h['frequency'][:], np.geomspace(100, 110, 5))
            self.assertNotEqual(float(h['frequency'][1]), 101.0)
        self.assertGreater(summary['grid_comparison']['X']['max_relative_difference'], 0)

    def test_x_y_grid_mismatch_rejected(self):
        self.fill([100, 101, 103, 106, 110])
        self.y.write_text(self.y.read_text().replace('103 ', '104 '))
        with self.assertRaisesRegex(VendorConversionError, 'frequency mismatch'):
            self.convert()

    def test_existing_outputs_not_overwritten(self):
        self.fill(np.geomspace(100, 110, 5))
        self.convert()
        with self.assertRaises(FileExistsError):
            self.convert()


if __name__ == '__main__':
    unittest.main()
