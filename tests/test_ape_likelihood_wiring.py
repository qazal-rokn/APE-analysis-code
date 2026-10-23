"""Necessary contract checks; external vendor and C are not executed here."""
import unittest
from unittest.mock import patch
import numpy as np
from physics.axion_halo import maxwell_halo_bin_weights
from ape_c_lpsd.analysis import ape_likelihood


class IntegratedTests(unittest.TestCase):
    def test_halo_mass_and_inadequate_grid(self):
        edges=np.linspace(99.999,100.01,500)
        f,h,m=maxwell_halo_bin_weights(
            frequency_edges_hz=edges,rest_frequency_hz=100,
            velocity_scale_m_s=220000,light_speed_m_s=299792458,
            min_captured_probability=.999)
        self.assertAlmostEqual(h.sum(),1.,places=12)
        self.assertGreater(m,.999)
        with self.assertRaises(ValueError):
            maxwell_halo_bin_weights(
                frequency_edges_hz=np.linspace(100,100.00001,5),
                rest_frequency_hz=100,velocity_scale_m_s=220000,
                light_speed_m_s=299792458,min_captured_probability=.999)

    def test_ape_template_vendor_contract(self):
        inputs=dict(observed_ps=np.full((2,3),2.e-16),
                    background_ps=np.full((2,3),1.e-16),
                    skewnormal_parameters=np.array([[0.,1.,0.],[0.,1.,0.]]),
                    frequency_support_hz=[100.],
                    normalized_halo_weights=[1.],
                    window_power_response=np.ones((2,3,1)),
                    quadrature_power_fraction=[1.,.5],
                    reference_coupling_gev_inv=1.e-10,
                    finesse_s=1510.,tau_s=2.5e-9,tau1_s=1.e-8,
                    rho_dm_gev4=3.1e-42,unit_conversion=1.519e24)
        prepared=ape_likelihood.prepare_ape_likelihood(**inputs)
        np.testing.assert_allclose(
            prepared["peak_norm"]*prepared["peak_shape"][None,:],
            prepared["signal_template_ps"],rtol=1e-13)
        with patch.object(ape_likelihood,"log_likelihood",return_value=4.2) as call:
            value=ape_likelihood.log_likelihood_at_coupling(
                2e-10,prepared=prepared,external_root="/unused")
            self.assertEqual(value,4.2)
            self.assertAlmostEqual(call.call_args.args[0],4.)
        with self.assertRaises(ValueError):
            ape_likelihood.prepare_ape_likelihood(
                **dict(inputs,observed_ps=np.ones((2,2))))


if __name__=="__main__":
    unittest.main()
