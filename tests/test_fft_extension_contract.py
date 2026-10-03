"""Future-only unit tests for FFT extension; intentionally not executed in this stage."""
import unittest
import numpy as np
from ape_c_lpsd.fft_extension.config import FFTAnalysisOptions
from ape_c_lpsd.fft_extension.preflight import FFTPreflightError, preflight_full_record_fft
from ape_c_lpsd.fft_extension.spectrum import full_record_one_sided_psd

def _options(**changes):
    values = dict(requested_display_min_hz=1., requested_display_max_hz=40., fft_min_hz=1., fft_max_hz=20., c_lpsd_handoff_hz=20., c_lpsd_min_hz=20., c_lpsd_max_hz=40., detrend="none", window_policy="rectangular", kaiser_beta=None, one_sided=True, full_record_only=True, zero_padding_factor=1, memory_limit_bytes=10_000_000)
    values.update(changes)
    return FFTAnalysisOptions(**values)

class FFTContractTests(unittest.TestCase):
    def test_preflight_refuses_insufficient_memory(self):
        with self.assertRaises(FFTPreflightError):
            preflight_full_record_fft(sample_count=1_000_000, sample_rate_hz=100., options=_options(memory_limit_bytes=1))

    def test_one_sided_even_record_does_not_double_dc_or_nyquist(self):
        dc = full_record_one_sided_psd(np.ones(64), sample_rate_hz=64., options=_options())
        nyquist = full_record_one_sided_psd((-1.0) ** np.arange(64), sample_rate_hz=64., options=_options())
        self.assertEqual(dc.method_identity, "fft_full_record")
        self.assertAlmostEqual(dc.psd_v2_per_hz[0], 1.0)
        self.assertAlmostEqual(nyquist.psd_v2_per_hz[-1], 1.0)

    def test_kaiser_policy_is_explicitly_not_claimed_c_equivalent(self):
        result = full_record_one_sided_psd(np.arange(64, dtype=float), sample_rate_hz=64., options=_options(window_policy="numpy_kaiser_to_be_validated", kaiser_beta=20.))
        self.assertTrue(np.all(np.isfinite(result.psd_v2_per_hz)))
