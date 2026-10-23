"""Contract tests: vendor shape unchanged; only optical response replaced."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from ape_c_lpsd.analysis import vendor_shape_ape


class VendorShapeTests(unittest.TestCase):
    def setUp(self):
        self.f=np.array([100.,100.1,100.2])
        self.shape=np.array([.2,.5,.3])
        self.inputs=dict(observed_ps=np.full((2,3),2e-16),
                         background_ps=np.full((2,3),1e-16),
                         skewnormal_parameters=np.array([[0.,1.,0.],[0.,1.,0.]]),
                         frequency_bins_hz=self.f, vendor_peak_shape=self.shape,
                         quadrature_power_fraction=np.array([1.,.25]),
                         reference_coupling_gev_inv=1e-10,
                         finesse_s=1510.,tau_s=2.5e-9,tau1_s=1e-8,
                         rho_dm_gev4=3.1e-42,unit_conversion=1.519e24)

    def test_shape_exact_and_quadrature_response(self):
        actual=vendor_shape_ape.prepare_vendor_shape_ape_likelihood(**self.inputs)
        np.testing.assert_array_equal(actual["peak_shape"],self.shape)
        np.testing.assert_allclose(actual["peak_norm"][1],.25*actual["peak_norm"][0],
                                   rtol=1e-14)
        self.assertEqual(actual["Y"].shape,(2,3))
        with patch.object(vendor_shape_ape,"log_likelihood",return_value=2.) as stub:
            self.assertEqual(vendor_shape_ape.evaluate_vendor_shape_at_coupling(
                2e-10,prepared=actual,external_root="/unused"),2.)
            self.assertAlmostEqual(stub.call_args.args[0],4.)

    def test_load_original_file_and_reject_grid_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"peak_shape.npy"
            np.save(path,self.shape)
            got=vendor_shape_ape.load_vendor_peak_shape(
                path,frequency_bins_hz=self.f)
            np.testing.assert_array_equal(got,self.shape)
            with self.assertRaises(ValueError):
                vendor_shape_ape.load_vendor_peak_shape(
                    path,frequency_bins_hz=[100.,100.1])

    def test_invalid_normalization_rejected(self):
        with self.assertRaises(ValueError):
            vendor_shape_ape.prepare_vendor_shape_ape_likelihood(
                **dict(self.inputs,vendor_peak_shape=[.2,.5,.4]))


if __name__=="__main__":
    unittest.main()
