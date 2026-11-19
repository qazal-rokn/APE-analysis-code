"""Regression fixtures for APE diagnostic kernel/halo/QWP calculations.

Execute from APE_equivalent_external:
python -m unittest discover -s ape_post_lpsd/tests -p 'test_ape_diagnostic_integration.py' -v
"""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ape_native_kaiser_kernel import model_kernel
from ape_axion_halo_kernel_diagnostic import halo_modes
from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling

class DiagnosticKernelIntegrationTests(unittest.TestCase):
    def test_original_c_tone_peaks(self):
        p = dict(fs=20000., fmin=65., fmax=10000., jdes=5035956, psll=238.13)
        for tone, c_peak in ((5000.0022089, 63.966461226),
                             (7999.9970833, 38.36500839)):
            with self.subTest(tone=tone):
                k = model_kernel(np.array([tone]), np.array([tone]), **p)
                self.assertAlmostEqual(float(k[0, 0]) / c_peak, 1., delta=2e-6)

    def test_halo_weights_reproducible(self):
        p = dict(rest=5000., n_samples=5000, n_modes=64, seed=20261119,
                 sigma=229000*np.sqrt(3/2), escape=528000, lab_speed=229000)
        a, w, _ = halo_modes(**p)
        b, v, _ = halo_modes(**p)
        np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(w,v)
        self.assertAlmostEqual(float(w.sum()), 1.)
        self.assertTrue((a >= 5000.).all())

    def test_qwp_psd_peak_variance(self):
        f = np.array([5000.0022089, 5000.006])
        modes = np.array([5000.003, 5000.008])
        weights = np.array([0.4, 0.6])
        R = qwp_amplitude_per_coupling(modes, finesse_s=1500.,
             tau_s=2.5e-9, tau1_s=1e-8, rho_dm_gev4=3.1e-42,
             unit_conversion=1.519e24)
        K = model_kernel(f,modes,fs=20000.,fmin=65.,fmax=10000.,
                         jdes=5035956,psll=238.13)
        T = 0.5*(1e-10)**2*(K@(weights*R**2))
        self.assertTrue(np.isfinite(T).all())
        self.assertGreater(float(T.max()),0.)
