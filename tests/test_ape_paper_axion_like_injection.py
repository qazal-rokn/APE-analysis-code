"""Focused scientific sanity checks. No C runs or generated data."""
import unittest
import numpy as np
from tools.run_paper_axion_like_injection import (
    sample_modes, baseband_envelope, window_power_kernel,
    F_REST, V_ESC, V_LAB, c,
)


class AxionHaloInjectionChecks(unittest.TestCase):
    def test_velocity_broadening_and_seed(self):
        a,phase=sample_modes(np.random.default_rng(20261028),1024)
        b,other=sample_modes(np.random.default_rng(20261028),1024)
        np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(phase,other)
        self.assertTrue(np.all(a >= F_REST))
        self.assertLessEqual(float(np.max(a-F_REST)),
            F_REST*(V_ESC+V_LAB)**2/(2*c*c)+1e-12)
        self.assertGreater(float(np.std(a)), 0)
        self.assertLess(float(np.std(a)),.03)

    def test_one_mode_field_near_rest_frequency(self):
        freq=np.array([5000.003,5000.004])
        phase=np.array([0.,1.])
        t,z=baseband_envelope(freq,phase,np.ones(2),1.,f_rest=F_REST,
                              duration=1.,envelope_rate=8.)
        self.assertEqual(len(t),9)
        self.assertTrue(np.isfinite(z).all())
        self.assertGreater(float(np.std(z.real)),0)

    def test_kernel_is_coherent_power_not_probability(self):
        f=np.array([F_REST,F_REST+0.005,F_REST+0.01])
        w=window_power_kernel(f,np.array([F_REST]),q=1e6,psll=238.13)
        self.assertEqual(w.shape,(3,1))
        self.assertAlmostEqual(w[0,0],1.,places=11)
        self.assertGreater(w[1,0],0)
        self.assertGreater(w.sum(),1.)  # overlap MUST NOT be renormalized


if __name__ == "__main__":
    unittest.main()
