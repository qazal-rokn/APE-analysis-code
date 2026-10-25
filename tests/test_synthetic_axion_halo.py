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


    def test_peak_amplitude_power_matches_template(self):
        # Repeat the same random draw with zero coupling and subtract the
        # identical noise realization to isolate the injected signal.
        from physics.qwp_signal import qwp_amplitude_per_coupling
        with tempfile.TemporaryDirectory() as d:
            args=dict(fs_hz=64., duration_s=2., seed=9,
                      rest_frequency_hz=10., support_frequency_hz=[10., 10.5, 11.],
                      normalized_vendor_shape=[.2, .5, .3],
                      noise_rms_rad=1e-11, quadrature_power_fraction=[.5, .5],
                      finesse_s=1510., tau_s=2.5e-9, tau1_s=1e-8,
                      rho_dm_gev4=3.1e-42, unit_conversion=1.519e24)
            g=1e-10
            path0=Path(d)/"noise.h5"
            path1=Path(d)/"signal.h5"
            inject_discrete_halo(path0, coupling_gev_inv=0., **args)
            inject_discrete_halo(path1, coupling_gev_inv=g, **args)
            f=np.asarray(args["support_frequency_hz"])
            response=qwp_amplitude_per_coupling(
                f, finesse_s=args["finesse_s"], tau_s=args["tau_s"],
                tau1_s=args["tau1_s"], rho_dm_gev4=args["rho_dm_gev4"],
                unit_conversion=args["unit_conversion"])
            rng=np.random.default_rng(args["seed"])
            z=(rng.normal(size=len(f))+1j*rng.normal(size=len(f)))/np.sqrt(2)
            expected_total=np.sum(
                0.5*g**2*response**2*np.asarray(args["normalized_vendor_shape"])
                *np.abs(z)**2)
            with h5py.File(path0) as noise, h5py.File(path1) as injected:
                powers=[]
                for channel in ("rho_x", "rho_y"):
                    pure_signal=injected[channel][:]-noise[channel][:]
                    powers.append(np.mean(pure_signal**2))
            np.testing.assert_allclose(
                powers,
                expected_total*np.asarray(args["quadrature_power_fraction"]),
                rtol=1e-10, atol=0.)

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
