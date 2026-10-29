"""Explicit unit boundary for archived PSD products and APE amplitude reporting.

The native C direct-PS path is already rad**2 and MUST NOT use this conversion.
No axion Q or frequency-dependent bandwidth is supplied by this module.
"""
from __future__ import annotations
import numpy as np


def archived_psd_to_ps(psd_rad2_per_hz, *, effective_bandwidth_hz):
    """Convert an independently identified PSD to band power, rad**2.

    The supplied bandwidth must come from documented estimator settings,
    NOT an assumed f/Q_a factor. Does not read or change source files.
    """
    p = np.asarray(psd_rad2_per_hz, dtype=float)
    bw = np.asarray(effective_bandwidth_hz, dtype=float)
    if (p.size == 0 or not np.isfinite(p).all() or np.any(p < 0)
            or not np.isfinite(bw).all() or np.any(bw <= 0)):
        raise ValueError("PSD must be finite nonnegative and bandwidth finite positive")
    try:
        out = p * np.broadcast_to(bw, p.shape)
    except ValueError as exc:
        raise ValueError("Bandwidth must be broadcastable to PSD shape") from exc
    return out


def amplitude_from_power(signal_power_rad2):
    """Signal amplitude (radians), not raw measured signal-plus-noise amplitude."""
    p = np.asarray(signal_power_rad2, dtype=float)
    if p.size == 0 or not np.isfinite(p).all() or np.any(p < 0):
        raise ValueError("Signal power must be finite nonnegative")
    return np.sqrt(p)


def combined_xy_amplitude(ps_x_rad2, ps_y_rad2):
    """Diagnostic RSS amplitude in rad; do NOT substitute for two-channel likelihood."""
    x, y = np.asarray(ps_x_rad2, float), np.asarray(ps_y_rad2, float)
    if x.shape != y.shape or x.size == 0:
        raise ValueError("X/Y shapes must match")
    return amplitude_from_power(x + y)


def coupling_from_signal_amplitude(amplitude_rad, *, response_rad_per_gev_inv):
    """Apply the caller-validated APE amplitude response: g = A / R.

    The amplitude must be inferred SIGNAL amplitude after the power fit, not
    sqrt(raw observed PS). The response includes the approved RMS/peak
    convention and any quadrature partition; not inferred here.
    """
    a = np.asarray(amplitude_rad, float)
    r = np.asarray(response_rad_per_gev_inv, float)
    if (a.size == 0 or not np.isfinite(a).all() or np.any(a < 0)
            or not np.isfinite(r).all() or np.any(r <= 0)):
        raise ValueError("Require finite nonnegative signal amplitude and positive response")
    try:
        return a / np.broadcast_to(r, a.shape)
    except ValueError as exc:
        raise ValueError("Response must be broadcastable to amplitude") from exc


def coupling_from_fitted_power(mu, *, g_ref_gev_inv):
    """Equivalent amplitude conversion for vendor mu=(g/g_ref)**2, mu>=0."""
    mu = np.asarray(mu, float)
    g_ref = float(g_ref_gev_inv)
    if (mu.size == 0 or not np.isfinite(mu).all() or np.any(mu < 0)
            or not np.isfinite(g_ref) or g_ref <= 0):
        raise ValueError("Invalid fitted signal power or reference coupling")
    return coupling_from_signal_amplitude(
        amplitude_from_power(mu), response_rad_per_gev_inv=1.0 / g_ref)
