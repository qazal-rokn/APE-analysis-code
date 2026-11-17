"""Provisional APE QWP Case-2 signal model. All physical inputs mandatory.

The peak-amplitude and X/Y signal response still need injection validation.
"""
import numpy as np


def qwp_amplitude_per_coupling(frequency_hz, *, finesse_s, tau_s,
                               tau1_s, rho_dm_gev4, unit_conversion):
    """Nominal resonant high-finesse response, rad per GeV^-1."""
    f = np.asarray(frequency_hz, dtype=float)
    values = (finesse_s, tau_s, tau1_s, rho_dm_gev4, unit_conversion)
    if f.ndim != 1 or not f.size or not np.isfinite(f).all() or np.any(f <= 0):
        raise ValueError("frequency_hz must be nonempty positive 1D")
    if not all(np.isscalar(x) and np.isfinite(x) and x > 0 for x in values):
        raise ValueError("All QWP parameters must be explicit positive finite scalars")
    d = (np.pi / finesse_s)**2 + 4*np.sin(np.pi*f*tau1_s)**2
    return tau_s*unit_conversion*np.sqrt(2*rho_dm_gev4/d)


def qwp_signal_power_template(*, signal_frequency_hz, normalized_halo_weights,
                              window_power_response, quadrature_power_fraction,
                              reference_coupling_gev_inv, finesse_s, tau_s,
                              tau1_s, rho_dm_gev4, unit_conversion):
    """Reference-coupling PS (M,N), rad^2, with explicit window kernel (M,N,K)."""
    f = np.asarray(signal_frequency_hz, dtype=float)
    h = np.asarray(normalized_halo_weights, dtype=float)
    w = np.asarray(window_power_response, dtype=float)
    a = np.asarray(quadrature_power_fraction, dtype=float)
    if f.ndim != 1 or not f.size or h.shape != f.shape:
        raise ValueError("Frequency and probability arrays must match (K,)")
    if w.ndim != 3 or w.shape[2] != f.size or a.shape != (w.shape[0],):
        raise ValueError("Require response (M,N,K) and quadrature fraction (M,)")
    if not all(np.isfinite(x).all() for x in (h, w, a)):
        raise ValueError("Nonfinite signal model inputs")
    if (np.any(h < 0) or np.any(w < 0) or np.any(a < 0) or np.any(a > 1)
            or not np.isclose(h.sum(), 1, rtol=1e-12, atol=1e-12)):
        raise ValueError("Require nonnegative kernel, valid fractions, unit-sum halo")
    if (not np.isscalar(reference_coupling_gev_inv)
            or not np.isfinite(reference_coupling_gev_inv)
            or reference_coupling_gev_inv <= 0):
        raise ValueError("Positive reference coupling required")
    r = qwp_amplitude_per_coupling(
        f, finesse_s=finesse_s, tau_s=tau_s, tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4, unit_conversion=unit_conversion)
    return (reference_coupling_gev_inv**2 / 2)*a[:, None]*np.einsum(
        "k,mnk->mn", h*r**2, w)
