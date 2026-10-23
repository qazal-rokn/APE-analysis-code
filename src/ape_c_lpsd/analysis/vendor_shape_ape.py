"""Adapter for the external scalar-DM peak-shape convention.

The vendor's peak_shape is a supplied normalized probability per C-LPSD
frequency bin. This module does not invent a halo or resample it. The shape
must already describe the final detector estimator binning (including any
window effects); if that condition is not met, use the general QWP kernel
interface instead and do not call this adapter.
"""
from pathlib import Path
import numpy as np
from physics.qwp_signal import qwp_amplitude_per_coupling
from ape_c_lpsd.analysis.vendor_statistics import log_likelihood


def load_vendor_peak_shape(path, *, frequency_bins_hz):
    """Read the vendor .npy peak shape unchanged, with exact grid-length checks."""
    frequency = np.asarray(frequency_bins_hz, dtype=float)
    if (frequency.ndim != 1 or frequency.size == 0 or not np.isfinite(frequency).all()
            or np.any(frequency <= 0) or np.any(np.diff(frequency) <= 0)):
        raise ValueError("Require strictly increasing positive C-LPSD bin frequencies")
    path = Path(path)
    if path.suffix != ".npy" or not path.is_file():
        raise ValueError("Supply an existing original vendor peak-shape .npy")
    shape = np.load(path, allow_pickle=False)
    return validate_vendor_peak_shape(shape, frequency)


def validate_vendor_peak_shape(peak_shape, frequency_bins_hz):
    """Reject non-matching bins rather than interpolate or re-normalize silently."""
    frequency = np.asarray(frequency_bins_hz, dtype=float)
    shape = np.asarray(peak_shape, dtype=float)
    if (frequency.ndim != 1 or not frequency.size
            or not np.isfinite(frequency).all()
            or np.any(frequency <= 0) or np.any(np.diff(frequency) <= 0)
            or shape.shape != frequency.shape or not np.isfinite(shape).all()
            or np.any(shape < 0) or not np.isclose(shape.sum(), 1., rtol=1e-12, atol=1e-12)):
        raise ValueError("Vendor peak shape must be nonnegative, unit sum and exactly grid-aligned")
    return shape.copy()


def prepare_vendor_shape_ape_likelihood(*, observed_ps, background_ps,
                                        skewnormal_parameters, frequency_bins_hz,
                                        vendor_peak_shape, quadrature_power_fraction,
                                        reference_coupling_gev_inv,
                                        finesse_s, tau_s, tau1_s,
                                        rho_dm_gev4, unit_conversion):
    """Return the exact vendor arrays, with only detector power replaced by APE.

    The input peak shape must already encode observed estimator-bin response.
    APE reference signal power: peak_norm[m,j] * peak_shape[j].
    This is a provisional peak/RMS and X/Y calibration convention.
    """
    f = np.asarray(frequency_bins_hz, dtype=float)
    shape = validate_vendor_peak_shape(vendor_peak_shape, f)
    data = np.asarray(observed_ps, dtype=float)
    bg = np.asarray(background_ps, dtype=float)
    pars = np.asarray(skewnormal_parameters, dtype=float)
    frac = np.asarray(quadrature_power_fraction, dtype=float)
    if (data.ndim != 2 or data.shape[1] != f.size or not data.size
            or bg.shape != data.shape or pars.shape != (data.shape[0], 3)
            or frac.shape != (data.shape[0],)
            or not all(np.isfinite(a).all() for a in (data, bg, pars, frac))
            or np.any(data <= 0) or np.any(bg <= 0)
            or np.any(pars[:, 1] <= 0) or np.any(frac < 0) or np.any(frac > 1)):
        raise ValueError("Observed and fitted PS must be positive finite (M,N); X/Y fractions explicit")
    g_ref = float(reference_coupling_gev_inv)
    if not np.isfinite(g_ref) or g_ref <= 0:
        raise ValueError("reference coupling must be positive")
    response = qwp_amplitude_per_coupling(
        f, finesse_s=finesse_s, tau_s=tau_s, tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4, unit_conversion=unit_conversion)
    norm = 0.5*(g_ref*response[None, :])**2 * frac[:, None]
    return dict(Y=np.log(data), bkg=np.log(bg), peak_norm=norm,
                peak_shape=shape, model_args=pars,
                reference_coupling_gev_inv=g_ref)


def evaluate_vendor_shape_at_coupling(coupling_gev_inv, *, prepared, external_root):
    """Call the original unchanged vendor likelihood at nonnegative coupling."""
    g = float(coupling_gev_inv)
    if not np.isfinite(g) or g < 0:
        raise ValueError("Physical axion coupling must be nonnegative finite")
    mu = (g / prepared["reference_coupling_gev_inv"])**2
    return log_likelihood(mu, prepared["Y"], prepared["bkg"],
                          prepared["peak_norm"], prepared["peak_shape"],
                          prepared["model_args"], external_root)
