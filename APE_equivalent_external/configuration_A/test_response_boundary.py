"""Synthetic interface tests only; NOT optical/PSD scientific validation.

Run from APE_equivalent_external:
python -m unittest discover -s configuration_A -p 'test_*.py' -v
"""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from configuration_A.build_response import (
    ACalibrationError, make_singlebin_peak_norm, build_response_file
)
from ape_post_lpsd.ape_signal_response import load_peak_norm
from configuration_A.run import validate_a_response_file, ConfigurationAInputError
from configuration_A.derive_qwp_calibration import derive_provisional_qwp_calibration
from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling


class TestSingleBinAPEBoundary(unittest.TestCase):
    def test_peak_to_variance_and_channel_projection(self):
        f = np.array([100., 200.])
        x, y = make_singlebin_peak_norm(
            f, reference_coupling_gev_inv=2.,
            peak_rotation_rad_per_gev_inv_x=[3., 4.],
            peak_rotation_rad_per_gev_inv_y=[1.5, 2.],
            c_psd_hz_inv_per_unit_variance_x=[5., 6.],
            c_psd_hz_inv_per_unit_variance_y=[7., 8.]
        )
        np.testing.assert_allclose(x, [90., 192.])
        np.testing.assert_allclose(y, [31.5, 64.])

    def test_bad_grid_or_kernel_is_rejected(self):
        kw = {
            "reference_coupling_gev_inv": 1.,
            "peak_rotation_rad_per_gev_inv_x": [1., 2.],
            "peak_rotation_rad_per_gev_inv_y": [1., 2.],
            "c_psd_hz_inv_per_unit_variance_x": [1., 1.],
            "c_psd_hz_inv_per_unit_variance_y": [1., 1.],
        }
        with self.assertRaises(ACalibrationError):
            make_singlebin_peak_norm([200., 100.], **kw)
        kw["c_psd_hz_inv_per_unit_variance_x"] = [-1., 1.]
        with self.assertRaises(ACalibrationError):
            make_singlebin_peak_norm([100., 200.], **kw)

    def test_generated_file_loads_in_existing_unmodified_adapter(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            src, prov, out = (d / x for x in ("measured.npz", "source.json", "response.npz"))
            f = np.array([100., 200.])
            np.savez(
                src, frequency_hz=f,
                peak_rotation_rad_per_gev_inv_x=[3., 4.],
                peak_rotation_rad_per_gev_inv_y=[1.5, 2.],
                c_psd_hz_inv_per_unit_variance_x=[5., 6.],
                c_psd_hz_inv_per_unit_variance_y=[7., 8.],
            )
            provenance = {
                "response_source": "synthetic fixture",
                "calibration_version": "synthetic-v1",
                "native_c_kernel_source": "synthetic unit-variance PSD kernel",
                "single_bin_convention": "synthetic unit-variance tone",
                "peak_amplitude_calibration_source": "synthetic peak amplitude",
                "quadrature_projection_source": "synthetic X/Y projection",
                "response_status": "PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS",
                "reference_coupling_gev_inv": 2.,
            }
            prov.write_text(json.dumps(provenance), encoding="utf-8")
            meta = build_response_file(src, prov, out)
            self.assertFalse(meta["physical_upper_limits_authorized"])
            matrix, loaded = load_peak_norm(out, f)
            np.testing.assert_allclose(matrix, [[90., 192.], [31.5, 64.]])
            self.assertEqual(loaded["configuration"], "A_original_vendor_single_bin")
            self.assertEqual(validate_a_response_file(out)["configuration"], "A_original_vendor_single_bin")
            with self.assertRaises(ACalibrationError):
                build_response_file(src, prov, out)
            with self.assertRaises(Exception):
                load_peak_norm(out, np.array([101., 200.]))

    def test_reused_qwp_model_transfers_to_single_bin_interface(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            f = np.array([100., 200.])
            src, optics, output = (d / x for x in ("kernel.npz", "optics.json", "prepared.npz"))
            np.savez(src, frequency_hz=f,
                     quadrature_amplitude_projection_x=[1., 0.5],
                     quadrature_amplitude_projection_y=[0., 0.5],
                     c_psd_hz_inv_per_unit_variance_x=[5., 6.],
                     c_psd_hz_inv_per_unit_variance_y=[7., 8.])
            parameters = dict(finesse_s=1400., tau_s=1.,
                              tau1_s=2., rho_dm_gev4=1., unit_conversion=1.)
            optics.write_text(json.dumps(parameters), encoding="utf-8")
            derive_provisional_qwp_calibration(src, optics, output)
            theoretical = qwp_amplitude_per_coupling(f, **parameters)
            with np.load(output, allow_pickle=False) as z:
                np.testing.assert_array_equal(z["peak_rotation_rad_per_gev_inv_x"], theoretical*[1., 0.5])
                np.testing.assert_array_equal(z["peak_rotation_rad_per_gev_inv_y"], theoretical*[0., 0.5])

    def test_configuration_b_response_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "b_response.npz"
            np.savez(bad, frequency_hz=np.array([100., 200.]),
                     peak_norm_x=np.ones(2), peak_norm_y=np.ones(2),
                     metadata_json=json.dumps({"configuration": "B_multibin"}))
            with self.assertRaises(ConfigurationAInputError):
                validate_a_response_file(bad)


if __name__ == "__main__":
    unittest.main()
