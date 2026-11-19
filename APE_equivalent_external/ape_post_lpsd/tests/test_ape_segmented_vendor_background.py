"""Real-vendor likelihood interface tests for background segment boundaries."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ape_provisional_qwp_multibin_vendor import build_vendor_background_inputs
from scalar_dark_matter_lpsd.lvk.stats import log_likelihood


def table():
    return pd.DataFrame([
        dict(fmin=87.0,fmax=88.0,x_knots=[-1.],y_knots=[1.],
             alpha_skew=0.0,loc_skew=0.0,sigma_skew=0.5,chi_sqr=1.),
        dict(fmin=88.0,fmax=89.0,x_knots=[-2.],y_knots=[2.],
             alpha_skew=0.5,loc_skew=0.01,sigma_skew=0.7,chi_sqr=2.)
    ])


def fake_bkg(logf, knots):
    return np.full_like(logf, knots[0], dtype=float)


class SegmentedVendorTests(unittest.TestCase):
    @patch('ape_provisional_qwp_multibin_vendor.bkg_model',side_effect=fake_bkg)
    def test_cross_boundary_all_bins_preserved(self,_):
        f=np.array([87.99,88.0,88.01])
        x=build_vendor_background_inputs(table(),f,10.)
        self.assertEqual(x['mode'],'per_bin_segment')
        self.assertEqual(x['pars'].shape,(3,3))
        self.assertEqual(x['logB_matrix'].shape,(3,1))
        np.testing.assert_array_equal(x['logB_vector'],[-1.,-2.,-2.])
        T=np.ones(3)*0.01
        observed=np.log(np.exp(x['logB_vector'])+0.2*T)
        result=log_likelihood(0.2,observed[:,None],
             x['logB_matrix'],T[:,None],np.ones(1),x['pars'])
        self.assertTrue(np.isfinite(result))

    @patch('ape_provisional_qwp_multibin_vendor.bkg_model',side_effect=fake_bkg)
    def test_same_segment_keeps_original_interface(self,_):
        f=np.array([87.1,87.2,87.3])
        x=build_vendor_background_inputs(table(),f,10.)
        self.assertEqual(x['mode'],'single_segment')
        self.assertEqual(x['logB_matrix'].shape,(1,3))
        self.assertEqual(x['pars'].shape,(1,3))
        T=np.array([.01,.02,.03])
        observed=np.log(np.exp(x['logB_vector'])+0.2*T)
        original=log_likelihood(0.2,observed[None,:],
              x['logB_matrix'],T[None,:],np.ones(3),x['pars'])
        split=log_likelihood(0.2,observed[:,None],
              x['logB_vector'][:,None],T[:,None],np.ones(1),
              np.repeat(x['pars'],3,axis=0))
        self.assertAlmostEqual(float(original),float(split),places=11)

    @patch('ape_provisional_qwp_multibin_vendor.bkg_model',side_effect=fake_bkg)
    def test_missing_coverage_rejected(self,_):
        df=table()
        df.loc[0,'fmax']=87.8
        with self.assertRaisesRegex(RuntimeError,'uncovered'):
            build_vendor_background_inputs(df,np.array([87.7,87.9,88.1]),10.)

    @patch('ape_provisional_qwp_multibin_vendor.bkg_model',side_effect=fake_bkg)
    def test_poor_fit_rejected(self,_):
        df=table()
        df.loc[1,'chi_sqr']=50.
        with self.assertRaisesRegex(RuntimeError,'chi-square'):
            build_vendor_background_inputs(df,np.array([87.99,88.01]),10.)

if __name__=='__main__':
    unittest.main()
