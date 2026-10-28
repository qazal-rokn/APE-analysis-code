"""Focused tests for paper-halo time-domain injection; no C execution."""
import unittest
import numpy as np
from tools.run_paper_axion_like_injection import (
    F_REST, Q, sample_modes, baseband_envelope, window_power_kernel,
)


class PaperAxionLikeInjectionTests(unittest.TestCase):
    def test_boosted_truncated_velocity_modes_and_seed(self):
        f1,p1=sample_modes(np.random.default_rng(17),128)
        f2,p2=sample_modes(np.random.default_rng(17),128)
        np.testing.assert_array_equal(f1,f2)
        np.testing.assert_array_equal(p1,p2)
        self.assertEqual(len(f1),128)
        self.assertTrue(np.all(f1>=F_REST))
        self.assertTrue(np.all(np.diff(f1)>=0))
        self.assertGreater(f1.max(),f1.min())

    def test_off_grid_coherent_superposition(self):
        f=np.array([F_REST+.001,F_REST+.004])
        phases=np.array([0.,np.pi/2])
        time,env=baseband_envelope(f,phases,np.ones(2),2.,f_rest=F_REST,duration=1.)
        self.assertEqual(len(time),9)
        np.testing.assert_allclose(env[0],np.sqrt(2)*(1+1j),atol=1e-12)
        self.assertGreater(np.abs(env[-1]),0)

    def test_kaiser_coherent_tone_on_and_off_bin(self):
        bins=np.array([4999.995,5000.,5000.005])
        w=window_power_kernel(bins,np.array([5000.]),q=Q,psll=238.13)
        self.assertEqual(w.shape,(3,1))
        self.assertAlmostEqual(w[1,0],1.,places=12)
        self.assertAlmostEqual(w[0,0],.72697,delta=5e-5)
        self.assertAlmostEqual(w[2,0],.72697,delta=5e-5)


if __name__=="__main__":
    unittest.main()
