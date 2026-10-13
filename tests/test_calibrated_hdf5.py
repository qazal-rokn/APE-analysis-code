"""Synthetic validation for calibrated HDF5 writing."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ape_c_lpsd.calibrated_hdf5 import (
    create_calibrated_hdf5,
    validate_completed_hdf5,
    write_calibrated_chunk,
)


class CalibratedHDF5Tests(unittest.TestCase):
    def test_writes_aligned_float64_datasets_and_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibrated.h5"
            h5 = create_calibrated_hdf5(
                path,
                sample_count=8,
                sample_rate_hz=20000.0,
                manifest_sha256="abc123",
                calibration_method="ape_rms_quadrature_calibration_v1",
                calibration_parameters={
                    "cavity_transmission_w_per_v": 0.020,
                    "responsivity_a_per_w": 0.7,
                    "transimpedance_gain_v_per_a": 1.0e6,
                    "second_harmonic_pd_response": 0.7,
                },
                source_parts=[
                    {"path": "a.dxd", "raw_start": 10, "raw_stop": 18}
                ],
                code_commit="deadbeef",
                chunk_samples=4,
            )
            try:
                next_offset = write_calibrated_chunk(
                    h5,
                    output_start=0,
                    i0_w=np.linspace(0.01, 0.011, 4),
                    eta0=np.linspace(0.02, 0.021, 4),
                    rho_x=np.linspace(1e-6, 2e-6, 4),
                    rho_y=np.linspace(-2e-6, -1e-6, 4),
                )
                self.assertEqual(next_offset, 4)

                next_offset = write_calibrated_chunk(
                    h5,
                    output_start=4,
                    i0_w=np.linspace(0.011, 0.012, 4),
                    eta0=np.linspace(0.021, 0.022, 4),
                    rho_x=np.linspace(2e-6, 3e-6, 4),
                    rho_y=np.linspace(-1e-6, 0.0, 4),
                )
                self.assertEqual(next_offset, 8)
                validate_completed_hdf5(h5, 8)

                self.assertEqual(h5["rho_x"].shape, (8,))
                self.assertEqual(h5["rho_x"].dtype, np.dtype(np.float64))
                self.assertEqual(h5["rho_x"].attrs["unit"], "rad")
                self.assertEqual(h5.attrs["manifest_sha256"], "abc123")
                self.assertEqual(h5.attrs["sample_count"], 8)
                self.assertFalse(bool(h5.attrs["x_y_combined"]))
            finally:
                h5.close()

    def test_refuses_nonfinite_calibrated_chunk(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibrated.h5"
            h5 = create_calibrated_hdf5(
                path,
                sample_count=2,
                sample_rate_hz=20000.0,
                manifest_sha256="abc123",
                calibration_method="ape_rms_quadrature_calibration_v1",
                calibration_parameters={},
                source_parts=[],
                code_commit="deadbeef",
                chunk_samples=2,
            )
            try:
                with self.assertRaisesRegex(ValueError, "rho_x contains non-finite"):
                    write_calibrated_chunk(
                        h5,
                        output_start=0,
                        i0_w=np.array([0.01, 0.01]),
                        eta0=np.array([0.02, 0.02]),
                        rho_x=np.array([1e-6, np.nan]),
                        rho_y=np.array([0.0, 0.0]),
                    )
            finally:
                h5.close()


if __name__ == "__main__":
    unittest.main()
