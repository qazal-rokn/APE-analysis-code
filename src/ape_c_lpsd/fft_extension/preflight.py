"""No-data resource and range preflight for a full-record FFT."""
from __future__ import annotations
from dataclasses import dataclass
import math
from .config import FFTAnalysisOptions

class FFTPreflightError(RuntimeError):
    pass

@dataclass(frozen=True)
class FFTPreflight:
    sample_count: int
    sample_rate_hz: float
    nyquist_hz: float
    bin_width_hz: float
    rfft_bin_count: int
    expected_memory_bytes: int
    memory_limit_bytes: int
    fft_calculable_min_hz: float
    fft_calculable_max_hz: float
    c_lpsd_handoff_hz: float

def preflight_full_record_fft(*, sample_count: int, sample_rate_hz: float, options: FFTAnalysisOptions) -> FFTPreflight:
    """Refuse before any FFT allocation if full-record resources/ranges are invalid.

    Conservative memory model: float64 input, detrended work array, window,
    windowed work array, complex128 rFFT, PSD, and frequency vectors.
    """
    if sample_count < 2 or not math.isfinite(sample_rate_hz) or sample_rate_hz <= 0:
        raise FFTPreflightError("sample count and sample rate must be positive")
    nyquist = sample_rate_hz / 2.0
    if options.fft_max_hz > nyquist or options.requested_display_max_hz > nyquist:
        raise FFTPreflightError("FFT/display range exceeds Nyquist")
    if not (options.fft_min_hz <= options.c_lpsd_handoff_hz <= options.fft_max_hz):
        raise FFTPreflightError("handoff must lie within requested FFT range")
    bins = sample_count // 2 + 1
    expected = (4 * sample_count * 8) + (bins * 16) + (2 * bins * 8)
    if expected > options.memory_limit_bytes:
        raise FFTPreflightError(f"full-record FFT refused: expected_memory_bytes={expected}, memory_limit_bytes={options.memory_limit_bytes}")
    return FFTPreflight(sample_count, sample_rate_hz, nyquist, sample_rate_hz / sample_count, bins, expected, options.memory_limit_bytes, sample_rate_hz / sample_count, nyquist, options.c_lpsd_handoff_hz)
