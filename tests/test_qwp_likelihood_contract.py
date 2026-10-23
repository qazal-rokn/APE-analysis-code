"""Minimal synthetic verification of provisional QWP response."""
import unittest
import numpy as np
from physics.qwp_signal import qwp_amplitude_per_coupling, qwp_signal_power_template


class Tests(unittest.TestCase):
    def test_finesse_matches_original_losses(self):
        t1,t2,ts=600e-6,580e-6,745e-6
        f=np.geomspace(1,10000,25)
        config=dict(finesse_s=2*np.pi/(t1+t2+4*ts),
                    tau_s=2.5e-9,tau1_s=1e-8,
                    rho_dm_gev4=3.1e-42,unit_conversion=1.519e24)
        response=qwp_amplitude_per_coupling(f,**config)
        old=2.5e-9*1.519e24*np.sqrt(
            2*3.1e-42/(((t1+t2+4*ts)/2)**2+4*np.sin(np.pi*f*1e-8)**2))
        np.testing.assert_allclose(response,old,rtol=1e-13)

    def test_peak_power_scaling(self):
        config=dict(finesse_s=1510.,tau_s=2.5e-9,tau1_s=1e-8,
                    rho_dm_gev4=3.1e-42,unit_conversion=1.519e24)
        amplitude=qwp_amplitude_per_coupling([100.],**config)[0]
        ps=qwp_signal_power_template(
            signal_frequency_hz=[100.],normalized_halo_weights=[1.],
            window_power_response=np.ones((2,3,1)),
            quadrature_power_fraction=[1.,.25],
            reference_coupling_gev_inv=1e-10,**config)
        np.testing.assert_allclose(ps[0],np.full(3,.5*(1e-10*amplitude)**2))
        np.testing.assert_allclose(ps[1],.25*ps[0])

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            qwp_amplitude_per_coupling([-1.],finesse_s=1510.,
                tau_s=2.5e-9,tau1_s=1e-8,rho_dm_gev4=3.1e-42,
                unit_conversion=1.519e24)


if __name__=="__main__":
    unittest.main()
