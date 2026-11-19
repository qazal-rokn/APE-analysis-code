"""Synthetic checks for linear-PSD sum before unchanged vendor background fitting."""
import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np

from configuration_A.combine_quadrature_psd import (
    CombinedPSDInputError, combine_psd_hdf5,
)
from scalar_dark_matter_lpsd.LPSDIO import LPSDOutput


class TestSummedQuadratures(unittest.TestCase):
    @staticmethod
    def make_h5(path, channel, psd, *, t0=0):
        freq = np.array([100., 101., 102.], dtype=float)
        with h5py.File(path, "w") as h:
            h.attrs["channel"] = channel
            h.attrs["spectrum_unit"] = "rad^2/Hz"
            h.attrs["sample_rate_hz"] = 20000.
            h.attrs["segment_t0"] = t0
            h.attrs["segment_t1"] = 15435
            h.create_dataset("frequency", data=freq)
            h.create_dataset("PSD", data=np.array(psd, dtype=float))
            d = h.create_dataset("logPSD", data=np.log(psd))
            for key, value in (("fmin",freq[0]),("fmax",freq[-1]),
                               ("fs",20000.),("Jdes",len(freq))):
                d.attrs[key] = value
            h.create_dataset("N_avg", data=[3.,3.,3.])
            h.create_dataset("N_kernel", data=[10.,10.,10.])

    def test_sum_is_linear_then_log_and_vendor_reads_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            x,y = d/"APE_0_15435_X.h5",d/"APE_0_15435_Y.h5"
            self.make_h5(x,"X",[2.,3.,5.])
            self.make_h5(y,"Y",[7.,11.,13.])
            info = combine_psd_hdf5(x,y,d/"summed",chunk_size=2)
            path = Path(info["output_file"])
            with h5py.File(path) as h:
                np.testing.assert_array_equal(h["PSD"][:],[9.,14.,18.])
                np.testing.assert_allclose(h["logPSD"][:],np.log([9.,14.,18.]))
                np.testing.assert_array_equal(h["N_avg"][:],[3.,3.,3.])
                self.assertEqual(h.attrs["channel"],"Q")
            vendor = LPSDOutput(str(path))
            np.testing.assert_array_equal(vendor.freq,[100.,101.,102.])
            np.testing.assert_allclose(vendor.logPSD,np.log([9.,14.,18.]))
            with self.assertRaises(CombinedPSDInputError):
                combine_psd_hdf5(x,y,d/"summed",chunk_size=2)

    def test_mismatched_segment_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            x,y = d/"X.h5",d/"Y.h5"
            self.make_h5(x,"X",[2.,3.,5.])
            self.make_h5(y,"Y",[7.,11.,13.],t0=1)
            with self.assertRaises(CombinedPSDInputError):
                combine_psd_hdf5(x,y,d/"sum")

    def test_log_psd_cannot_be_directly_summed(self):
        self.assertNotAlmostEqual(
            np.log(2.)+np.log(7.), np.log(2.+7.)
        )


if __name__=="__main__":
    unittest.main()
