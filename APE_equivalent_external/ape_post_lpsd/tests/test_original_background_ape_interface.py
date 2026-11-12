"""Execute the original vendor background fitter on independent APE X/Y inputs.

Original vendor scientific functions are imported and used unchanged.
Requires the upstream project's normal Python dependencies.
"""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ape_post_lpsd.prepare_vendor_lpsd import prepare_pair
from scalar_dark_matter_lpsd.LPSDIO import LPSDOutput, LPSDJSONIO
from scalar_dark_matter_lpsd.lvk import fit_background


class OriginalBackgroundAPEInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.folder = Path(self.workspace.name)
        self.freq = np.geomspace(65, 220, 750)
        rng = np.random.default_rng(821)
        x = self.folder / "rho_x.lpsd"
        y = self.folder / "rho_y.lpsd"
        self.psds = {}
        for channel, path, offset, slope in (("X", x, -24.0, 0.25), ("Y", y, -24.1, 0.23)):
            psd = np.exp(offset + slope * np.log(self.freq / 65) + rng.normal(0, .2, len(self.freq)))
            self.psds[channel] = psd
            with path.open("x", encoding="utf-8") as stream:
                for f, p in zip(self.freq, psd):
                    stream.write(f"{f:.17g} {p:.17g} 12 3\n")
        self.result = prepare_pair(x, y, self.folder / "handoff", fs=1000.,
                                   t0=100, t1=200, expected_bins=750)

    def test_vendor_reader_preserves_native_spectra(self):
        for channel in ("X", "Y"):
            data = LPSDOutput(self.result["output"][channel])
            np.testing.assert_array_equal(data.freq, self.freq)
            np.testing.assert_allclose(np.exp(data.logPSD), self.psds[channel], rtol=2e-14, atol=0)

    def test_vendor_background_fitter_produces_two_original_json_records(self):
        out = self.folder / "original_background.json"
        for channel in ("X", "Y"):
            fit_background.fit_background_in_file(
                self.result["output"][channel], str(out),
                n_processes=1, buffer=20, pruning=1, segment_size=1000,
                ana_fmin=65., ana_fmax=220., k_min=4, k_pruning=1,
                k_plateau=0, bin_scaling=.2, cfd_fraction=.1,
                kernel_size=30, verbose=False, plot_path=None)
        data = LPSDJSONIO(str(out))
        self.assertEqual(len(data.data), 2)
        names = {"fmin", "fmax", "x_knots", "y_knots", "alpha_skew",
                 "loc_skew", "sigma_skew", "chi_sqr"}
        for channel in ("X", "Y"):
            label = data.get_label(self.result["output"][channel])
            self.assertEqual(label, f"bkginfo_100_200_{channel}")
            fit = data.get_df(label)
            self.assertEqual(len(fit), 1)
            self.assertTrue(names.issubset(fit.columns))
            self.assertTrue(np.isfinite(float(fit.iloc[0]["chi_sqr"])))
            self.assertGreater(float(fit.iloc[0]["sigma_skew"]), 0)


if __name__ == "__main__":
    unittest.main()
