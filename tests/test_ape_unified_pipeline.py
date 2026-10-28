"""Opt-in C-to-original-vendor orchestration checks; no C jobs or real data."""
import unittest
from unittest.mock import patch
import numpy as np
from tools.run_ape_pipeline import run_vendor_mle_stage


class UnifiedPipelineTests(unittest.TestCase):
    def setUp(self):
        self.f = np.linspace(4995, 5008, 1301)
        self.obs = np.ones((2, 1301))*1e-19
        self.fits = [
            {"background_psd": np.ones(1301)*1e-19,
             "skewnormal_loc_scale_alpha": np.array([0., 0.7, 0.])}
            for _ in range(2)]
        self.conf = {
            "peak_shape_file": None,
            "n_effective_trials": None,
            "vendor_optimizer": {"min_log10mu": -8, "max_log10mu": 3,
                                 "n_seed_lkl": 32},
            "local_sigma_threshold": 3.,
            "physics": {
                "quadrature_power_fraction": [0.5,0.5],
                "reference_coupling_gev_inv": 1e-10,
                "finesse_s": 1510., "tau_s": 2.5e-9, "tau1_s": 1e-8,
                "rho_dm_gev4": 3.1e-42, "unit_conversion": 1.519e24,
            }}
    def test_calls_original_vendor_bridge_without_recomputing_spectrum(self):
        result = {"selection": {"candidates": []}, "fits": []}
        with patch("tools.run_ape_pipeline.fit_and_report_single_bin_candidates",
                   return_value=result) as bridge:
            out=run_vendor_mle_stage(self.f, self.obs, self.fits,
                                     np.array([500,501]), self.conf, "/source")
        self.assertIs(out, result)
        kw=bridge.call_args.kwargs
        np.testing.assert_allclose(kw["candidate_frequencies_hz"],self.f[[500,501]])
        self.assertEqual(len(kw["prepared_windows"]),2)
        for p in kw["prepared_windows"]:
            self.assertEqual(p["Y"].shape,(2,1))
            self.assertEqual(p["peak_shape"].tolist(),[1.])
        self.assertEqual(kw["min_log10mu"],-8)
        self.assertEqual(kw["max_log10mu"],3)
        self.assertEqual(kw["n_seed_lkl"],32)

    def test_refuses_multibin_unreviewed_template(self):
        self.conf["peak_shape_file"]="halo.npy"
        with self.assertRaises(ValueError):
            run_vendor_mle_stage(self.f,self.obs,self.fits,[500],self.conf,"/source")

    def test_refuses_unvalidated_global_trial_factor(self):
        self.conf["n_effective_trials"]=10
        with self.assertRaises(ValueError):
            run_vendor_mle_stage(self.f,self.obs,self.fits,[500],self.conf,"/source")

    def test_refuses_invalid_candidate_indices(self):
        with self.assertRaises(ValueError):
            run_vendor_mle_stage(self.f,self.obs,self.fits,[500,500],self.conf,"/source")


if __name__=="__main__":
    unittest.main()
