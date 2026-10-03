"""Canonical APE chunk calibration using the approved RMS equations."""
from __future__ import annotations

import numpy as np

from .ingestion.config import CalibrationConfig
from .ingestion.models import APEDXDChunk, APECalibratedChunk


class CalibrationError(ValueError):
    pass


def _positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise CalibrationError(f"{name} must be positive and finite")
    return value


def calibrate_chunk(chunk: APEDXDChunk, calibration: CalibrationConfig) -> APECalibratedChunk:
    """Calibrate one synchronized raw chunk without smoothing or X/Y combination."""
    c_trans = _positive_finite("cavity_transmission_w_per_v", calibration.cavity_transmission_w_per_v)
    q = _positive_finite("responsivity_a_per_w", calibration.responsivity_a_per_w)
    gain = _positive_finite("transimpedance_gain_v_per_a", calibration.transimpedance_gain_v_per_a)
    h2 = _positive_finite("second_harmonic_pd_response", calibration.second_harmonic_pd_response)

    x1 = np.asarray(chunk.x_1f, dtype=np.float64)
    y1 = np.asarray(chunk.y_1f, dtype=np.float64)
    x2 = np.asarray(chunk.x_2f, dtype=np.float64)
    y2 = np.asarray(chunk.y_2f, dtype=np.float64)
    vtrans = np.asarray(chunk.cavity_transmission, dtype=np.float64)

    i0_w = c_trans * vtrans
    r2f = np.hypot(x2, y2)

    valid = (
        np.isfinite(x1)
        & np.isfinite(y1)
        & np.isfinite(x2)
        & np.isfinite(y2)
        & np.isfinite(vtrans)
        & np.isfinite(i0_w)
        & np.isfinite(r2f)
        & (i0_w > 0.0)
        & (r2f > 0.0)
    )

    eta0 = np.full(chunk.n_samples, np.nan, dtype=np.float64)
    arg = np.full(chunk.n_samples, np.nan, dtype=np.float64)
    arg[valid] = 2.0 * np.sqrt(2.0) * r2f[valid] / (h2 * q * gain * i0_w[valid])
    valid &= np.isfinite(arg) & (arg > 0.0)
    eta0[valid] = np.sqrt(arg[valid])

    denominator = np.full(chunk.n_samples, np.nan, dtype=np.float64)
    denominator[valid] = np.sqrt(2.0) * q * gain * i0_w[valid] * eta0[valid]
    valid &= np.isfinite(denominator) & (denominator > 0.0)

    rho_x = np.full(chunk.n_samples, np.nan, dtype=np.float64)
    rho_y = np.full(chunk.n_samples, np.nan, dtype=np.float64)
    np.divide(x1, denominator, out=rho_x, where=valid)
    np.divide(y1, denominator, out=rho_y, where=valid)

    valid &= np.isfinite(eta0) & np.isfinite(rho_x) & np.isfinite(rho_y)
    eta0[~valid] = np.nan
    rho_x[~valid] = np.nan
    rho_y[~valid] = np.nan

    i0_out = i0_w.astype(np.float64, copy=True)
    i0_out[~np.isfinite(i0_out)] = np.nan

    metadata = {
        "method": "ape_rms_quadrature_calibration_v1",
        "input_quadrature_convention": "RMS voltage",
        "equations": {
            "i0": "I0=C_trans*V_trans",
            "r2f": "R2f=hypot(X2f,Y2f)",
            "eta0": "sqrt(2*sqrt(2)*R2f/(H2*q*G*I0))",
            "rho_x": "X1f/(sqrt(2)*q*G*I0*eta0)",
            "rho_y": "Y1f/(sqrt(2)*q*G*I0*eta0)",
        },
        "parameters": {
            "cavity_transmission_w_per_v": c_trans,
            "responsivity_a_per_w": q,
            "transimpedance_gain_v_per_a": gain,
            "second_harmonic_pd_response": h2,
        },
        "invalid_policy": "preserve index; calibrated outputs NaN; valid_mask false",
        "x_y_combined": False,
    }

    return APECalibratedChunk(
        source_file=chunk.source_file,
        start_index=chunk.start_index,
        stop_index=chunk.stop_index,
        sample_rate_hz=chunk.sample_rate_hz,
        i0_w=i0_out,
        eta0=eta0,
        rho_x=rho_x,
        rho_y=rho_y,
        valid_mask=valid.astype(bool, copy=False),
        calibration_metadata=metadata,
    )
