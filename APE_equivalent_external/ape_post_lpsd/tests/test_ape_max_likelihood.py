"""Small synthetic integration checks using original vendor optimizer/IO."""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, mock
import json
import sys
import unittest
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO, MaxLKLIO
from scalar_dark_matter_lpsd.lvk import max_lkl_executable as vendor
from ape_post_lpsd.ape_max_likelihood import run_single_quadrature, APELikelihoodError


class TestAPEOriginalLikelihood(TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.freq = np.array([100.0, 101.0, 102.0])
        self.psd = {"X": np.array([1.12, 1.08, 1.07]), "Y": np.array([1.05, 1.04, 1.04])}
        self.spectrum = {}
        for channel in ("X", "Y"):
            name = self.root / f"APE_100_200_{channel}.h5"
            self.spectrum[channel] = name
            import h5py
            with h5py.File(name, 'x') as h5:
                h5.create_dataset('frequency', data=self.freq)
                ds=h5.create_dataset('logPSD',data=np.log(self.psd[channel]))
                ds.attrs['fmin']=self.freq[0]
                ds.attrs['fmax']=self.freq[-1]
                ds.attrs['fs']=1000.
                ds.attrs['Jdes']=len(self.freq)
        self.bkg = self.root/'bkg.json'
        records = LPSDJSONIO(str(self.bkg))
        df = pd.DataFrame([{ 'fmin':99.,'fmax':103.,'x_knots':list(np.log([99.,100.,102.,103.])),
                             'y_knots':[0.,0.,0.,0.], 'alpha_skew':0.,'loc_skew':0.,
                             'sigma_skew':0.15,'chi_sqr':1.0}])
        for channel in ('X','Y'):
            records.update_file(str(self.spectrum[channel]),df)
        self.response=self.root/'response.npz'
        self.meta={'mu_definition':'synthetic amplitude (not axion coupling)',
                   'response_source':'synthetic unit test', 'calibration_version':'synthetic',
                   'spectral_shape_convention':'unit sum applied in vendor likelihood',
                   'response_units':'rad^2/Hz per unit mu'}
        np.savez(self.response, frequency_hz=self.freq,
                 peak_norm_x=np.array([0.8,0.8,0.8]),peak_norm_y=np.array([0.4,0.4,0.4]),
                 metadata_json=json.dumps(self.meta))

    def invoke(self, channel='X', out=None, **kwargs):
        return run_single_quadrature(self.spectrum.get(channel, self.spectrum['X']), self.bkg, self.response,
            out if out is not None else self.root / f'{channel}_fit.npy', channel=channel,
            ana_fmin=100., ana_fmax=103., min_log10mu=-4, max_log10mu=1,
            max_frequencies=1, n_seed_lkl=80, **kwargs)

    def test_x_unchanged_vendor_optimizer_exact_parity(self):
        # Patch LIGO-specific function to fail if accidentally called.
        with mock.patch.object(vendor.stats, 'get_beta',side_effect=AssertionError('LIGO get_beta forbidden')):
            manifest=self.invoke()
        self.assertEqual(manifest['joint_xy_likelihood_used'],False)
        saved=MaxLKLIO(infile=str(self.root/'X_fit.npy')).data
        self.assertEqual(saved.shape,(1,5))
        expected=vendor.calc_max_lkl(np.log([[self.psd['X'][0]]]),np.zeros((1,1)),
             np.array([[0.,0.15,0.]]), np.array([[0.8]]), np.array([1.]),
             min_log10mu=-4,max_log10mu=1,n_seed_lkl=80)
        np.testing.assert_allclose(saved[0,1:],expected, equal_nan=True,atol=1e-12,rtol=1e-10)
        self.assertEqual(saved[0,0],100.)

    def test_y_runs_separately_without_x_as_second_detector(self):
        self.invoke('Y')
        saved=MaxLKLIO(infile=str(self.root/'Y_fit.npy')).data
        self.assertEqual(saved.shape,(1,5))
        self.assertEqual(saved[0,0],100.)
        self.assertTrue(np.isfinite(saved[0,2]))

    def test_refuses_nonmatching_response_grid(self):
        np.savez(self.response,frequency_hz=self.freq+0.1,
            peak_norm_x=np.array([.8]*3),peak_norm_y=np.array([.4]*3),
            metadata_json=json.dumps(self.meta))
        with self.assertRaisesRegex(ValueError,'frequency grid'):
            self.invoke()
        self.assertFalse((self.root/'X_fit.npy').exists())

    def test_refuses_bad_channel(self):
        with self.assertRaisesRegex(APELikelihoodError,'one quadrature'):
            self.invoke('H1')

    def test_refuses_output_overwrite(self):
        self.invoke()
        with self.assertRaisesRegex(APELikelihoodError,'overwrite'):
            self.invoke()

    def test_rejects_no_background_coverage(self):
        records=LPSDJSONIO(str(self.bkg))
        df=records.get_df(records.get_label(filepath=str(self.spectrum['X'])))
        df['fmin']=200.
        records.update_file(str(self.spectrum['X']),df)
        with self.assertRaisesRegex(APELikelihoodError,'No frequencies survived'):
            self.invoke()

if __name__ == '__main__':
    unittest.main()
