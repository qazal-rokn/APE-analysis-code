"""Focused axion-like halo injection tests: mathematics only; no C and no outputs."""
import unittest
import numpy as np
from tools.run_paper_axion_like_injection import (
    F_REST, Q, PSLL, sample_modes, baseband_envelope, window_power_kernel,
)

class PaperHaloInjectionTests(unittest.TestCase):
    def test_source_grounded_velocity_field_and_reproducibility(self):
        first=sample_modes(np.random.default_rng(21),512)
        second=sample_modes(np.random.default_rng(21),512)
        for a,b in zip(first,second):
            np.testing.assert_array_equal(a,b)
        f,phase=first
        self.assertTrue(np.all(f>=F_REST))
        self.assertTrue(np.all(np.diff(f)>=0))
        self.assertTrue(np.all((phase>=0)&(phase<2*np.pi)))
        self.assertGreater(np.ptp(f), 0.001)

    def test_coherent_shared_mode_synthesis(self):
        f=np.array([F_REST+0.003,F_REST+0.005])
        phase=np.array([0.,np.pi/3])
        t,e=baseband_envelope(f,phase,np.ones(2),1.,
                              f_rest=F_REST,duration=2)
        self.assertEqual(len(t),len(e))
        self.assertAlmostEqual(e[0].real,
            np.cos(phase).sum()/np.sqrt(2), places=12)
        self.assertTrue(np.iscomplexobj(e))

    def test_c_kaiser_response_not_a_probability_shape(self):
        f=np.array([F_REST-0.005,F_REST,F_REST+0.005])
        power=window_power_kernel(f,np.array([F_REST]),q=Q,psll=PSLL)
        self.assertEqual(power.shape,(3,1))
        self.assertAlmostEqual(power[1,0],1.,places=12)
        self.assertGreater(power.sum(),1.)

if __name__=="__main__": unittest.main()
