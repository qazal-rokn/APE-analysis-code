"""Dense four-hour native-grid scan unit tests (no real APE files needed)."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ape_post_lpsd.ape_full_native_likelihood_scan import (
    fitted_background_on_native_grid,fit_signed_vendor_likelihood,
    eligible_native_indices,halo_benchmark)
from ape_post_lpsd.ape_native_kaiser_kernel import (
    model_kernel,_cached_kaiser_parameters)

class FullNativeLikelihoodTests(unittest.TestCase):
    def test_selected_native_bins(self):
        f=np.arange(65.,125.,1.)
        np.testing.assert_array_equal(eligible_native_indices(f,70.,115.,3),
                                      np.arange(5,51))

    def test_cached_kaiser_original_tones(self):
        p=dict(fs=20000.,fmin=65.,fmax=10000.,jdes=5035956,psll=238.13)
        for center,expected in [(5000.0022089,63.966461226),
                                (7999.9970833,38.36500839)]:
            got=model_kernel(np.array([center]),np.array([center]),**p)
            self.assertAlmostEqual(float(got[0,0])/expected,1.,delta=2e-6)
        self.assertGreater(_cached_kaiser_parameters.cache_info().hits,0)

    def test_fixed_halo_reproducible(self):
        a,w,_=halo_benchmark(1000,16,20261119)
        b,v,_=halo_benchmark(1000,16,20261119)
        np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(w,v)

    @patch("ape_post_lpsd.ape_full_native_likelihood_scan.bkg_model")
    def test_mixed_background(self,model):
        model.side_effect=lambda logf,k:np.full(len(logf),k[0])
        t=pd.DataFrame([dict(fmin=100.,fmax=100.2,chi_sqr=1,
              x_knots=[-40.],y_knots=[0.],alpha_skew=0.,
              loc_skew=0.,sigma_skew=0.5),
              dict(fmin=100.2,fmax=100.5,chi_sqr=1,
              x_knots=[-39.],y_knots=[0.],alpha_skew=1.,
              loc_skew=0.,sigma_skew=0.4)])
        f=np.array([100.1,100.2,100.3])
        b,pars,own,valid=fitted_background_on_native_grid(t,f)
        np.testing.assert_array_equal(b,[-40.,-39.,-39.])
        self.assertTrue(valid.all())
        self.assertEqual(len(set(own.tolist())),2)

    def test_signed_vendor_likelihood_finite(self):
        n=15
        b=np.log(np.ones(n)*3e-16)
        T=np.ones(n)*2e-19
        pars=np.c_[np.zeros(n),np.ones(n)*0.4,np.zeros(n)]
        observed=b.copy()
        z=fit_signed_vendor_likelihood(observed,b,pars,T)
        self.assertTrue(np.isfinite(z["lnL0"]))
        self.assertTrue(np.isfinite(z["lnLmax"]))
        self.assertGreaterEqual(z["lnLmax"],z["lnL0"])
        self.assertIsNotNone(z["mu_hat_signed"])
        self.assertTrue(np.isfinite(z["mu_hat_signed"]))

if __name__=="__main__":unittest.main()
