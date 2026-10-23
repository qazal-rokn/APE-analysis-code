"""Finite stochastic halo fixture tests, no C or original vendor execution."""
import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np
from ape_c_lpsd.analysis.synthetic_axion_halo import inject_discrete_halo

class HaloTests(unittest.TestCase):
    def test_seed_and_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            args=dict(fs_hz=64.,duration_s=2.,seed=9,rest_frequency_hz=10.,
                      support_frequency_hz=[10.,10.5,11.],
                      normalized_vendor_shape=[.2,.5,.3],
                      coupling_gev_inv=1e-10,noise_rms_rad=1e-11,
                      quadrature_power_fraction=[.5,.5],
                      finesse_s=1510.,tau_s=2.5e-9,tau1_s=1e-8,
                      rho_dm_gev4=3.1e-42,unit_conversion=1.519e24)
            p=Path(d)/"test.h5"
            result=inject_discrete_halo(p,**args)
            self.assertEqual(result["n_samples"],128)
            self.assertTrue(result["stochastic"])
            with h5py.File(p) as h:
                self.assertTrue(h.attrs["synthetic_only"])
                self.assertEqual(h["rho_x"].shape,(128,))
                self.assertEqual(h.attrs["model"],"stochastic_gaussian_vendor_shape")
                np.testing.assert_array_equal(h["vendor_shape_weights"][:],[.2,.5,.3])
            with self.assertRaises(FileExistsError):
                inject_discrete_halo(p,**args)

    def test_unresolved_halo_bins_refused(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                inject_discrete_halo(Path(d)/"bad.h5",fs_hz=64.,duration_s=2.,
                    seed=9,rest_frequency_hz=10.,
                    support_frequency_hz=[10.,10.1],
                    normalized_vendor_shape=[.5,.5],
                    coupling_gev_inv=1e-10,noise_rms_rad=1e-11,
                    quadrature_power_fraction=[.5,.5],
                    finesse_s=1510.,tau_s=2.5e-9,tau1_s=1e-8,
                    rho_dm_gev4=3.1e-42,unit_conversion=1.519e24)
if __name__=="__main__":
    unittest.main()
