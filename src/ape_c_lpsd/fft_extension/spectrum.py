"""Future execution kernel for a full-record, one-sided voltage PSD."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .config import FFTAnalysisOptions
from .preflight import FFTPreflight, preflight_full_record_fft

@dataclass(frozen=True)
class FullRecordFFTPSD:
    frequency_hz: np.ndarray
    psd_v2_per_hz: np.ndarray
    asd_v_per_sqrt_hz: np.ndarray
    preflight: FFTPreflight
    method_identity: str = "fft_full_record"

def _linear_detrend(values: np.ndarray) -> np.ndarray:
    index = np.arange(values.size, dtype=np.float64)
    centered = index - index.mean()
    slope = float(np.dot(centered, values - values.mean()) / np.dot(centered, centered))
    return values - (values.mean() + slope * centered)

def _window(sample_count: int, options: FFTAnalysisOptions) -> np.ndarray:
    if options.window_policy == "rectangular":
        return np.ones(sample_count, dtype=np.float64)
    # Explicitly not labelled C-equivalent. Numerical equivalence is a future test.
    return np.kaiser(sample_count, float(options.kaiser_beta))

def full_record_one_sided_psd(samples_v: np.ndarray, *, sample_rate_hz: float, options: FFTAnalysisOptions) -> FullRecordFFTPSD:
    """Calculate a future full-record one-sided density PSD for one validated V sequence.

    Callers must supply one uninterrupted, already-valid, correctly ordered channel.
    This routine neither cleans, fills, resamples, reorders, nor combines data.
    """
    values = np.asarray(samples_v, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
        raise ValueError("FFT input must be a finite one-dimensional validated voltage sequence")
    plan = preflight_full_record_fft(sample_count=values.size, sample_rate_hz=sample_rate_hz, options=options)
    if options.detrend == "none":
        detrended = values.copy()
    elif options.detrend == "mean":
        detrended = values - values.mean()
    else:
        detrended = _linear_detrend(values)
    window = _window(values.size, options)
    transform = np.fft.rfft(detrended * window)
    psd = (np.abs(transform) ** 2) / (sample_rate_hz * np.dot(window, window))
    if values.size % 2 == 0:
        psd[1:-1] *= 2.0
    else:
        psd[1:] *= 2.0
    frequency = np.fft.rfftfreq(values.size, d=1.0 / sample_rate_hz)
    return FullRecordFFTPSD(frequency, psd, np.sqrt(psd), plan)
