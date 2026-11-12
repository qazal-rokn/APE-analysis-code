"""Read-only calibrated quadrature review utilities.

This module extends the frozen reference scalar denominator to two separately
reported quadratures. It deliberately introduces no phase, gain-imbalance, or
leakage correction and never combines X/Y in its outputs.
"""
from __future__ import annotations

import numpy as np


def denominator_for_raw_samples(quality, raw_start: int, sample_count: int, parameters) -> np.ndarray:
    """Return the exact 1000:1 mapped frozen-reference denominator D(t)."""
    parent = np.arange(raw_start, raw_start + sample_count, dtype=np.int64) // parameters.samples_per_reduced_block
    i0 = np.asarray([record.i0_w for record in quality.records], dtype=np.float64)[parent]
    eta0 = np.asarray([record.eta0 for record in quality.records], dtype=np.float64)[parent]
    return (parameters.first_harmonic_rms_factor * i0 * eta0 *
            parameters.responsivity_a_per_w * parameters.transimpedance_gain_v_per_a)


def calibrated_quadratures(x_v, y_v, quality, raw_start: int, parameters, valid_mask):
    """Return D, X/D, and Y/D; invalid locations remain NaN, never filled."""
    x = np.asarray(x_v, dtype=np.float64)
    y = np.asarray(y_v, dtype=np.float64)
    valid = np.asarray(valid_mask, dtype=bool)
    if x.shape != y.shape or x.shape != valid.shape:
        raise ValueError("X, Y, and valid mask must have identical shapes")
    denominator = denominator_for_raw_samples(quality, raw_start, len(x), parameters)
    x_ellipticity = np.full(x.shape, np.nan, dtype=np.float64)
    y_ellipticity = np.full(y.shape, np.nan, dtype=np.float64)
    np.divide(x, denominator, out=x_ellipticity, where=valid)
    np.divide(y, denominator, out=y_ellipticity, where=valid)
    return denominator, x_ellipticity, y_ellipticity


def calibration_profile(parameters) -> dict[str, object]:
    return {
        "method": "reference_scalar_denominator_applied_separately_to_pem1_x_and_pem1_Y",
        "equation": "D=sqrt(2)*I0*eta0*responsivity*transimpedance_gain; X_ellipticity=X_voltage/D; Y_ellipticity=Y_voltage/D",
        "raw_to_reduced_mapping": f"raw_index // {parameters.samples_per_reduced_block}",
        "no_phase_gain_imbalance_or_leakage_correction": True,
        "parameters": {name: getattr(parameters, name) for name in (
            "sample_rate_hz", "reduced_block_size_seconds", "samples_per_reduced_block",
            "responsivity_a_per_w", "transimpedance_gain_v_per_a", "first_harmonic_rms_factor")},
        "input_unit": "V",
        "denominator_unit": "V",
        "output_unit": "1 (dimensionless ellipticity)",
    }
