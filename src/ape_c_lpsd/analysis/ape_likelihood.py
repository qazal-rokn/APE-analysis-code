"""APE-specific signal-power -> unchanged vendor likelihood connection.

All detector/halo/window choices must be supplied by the caller. This does
not establish the physical validity of those inputs or statistical coverage.
"""
import numpy as np
from physics.qwp_signal import qwp_signal_power_template
from ape_c_lpsd.analysis.vendor_statistics import factor_power_template, log_likelihood


def prepare_ape_likelihood(*, observed_ps, background_ps,
                           skewnormal_parameters, frequency_support_hz,
                           normalized_halo_weights, window_power_response,
                           quadrature_power_fraction, reference_coupling_gev_inv,
                           finesse_s, tau_s, tau1_s, rho_dm_gev4,
                           unit_conversion):
    """Build vendor inputs for user-provided PS (M,N), all in rad^2."""
    observed = np.asarray(observed_ps, dtype=float)
    background = np.asarray(background_ps, dtype=float)
    pars = np.asarray(skewnormal_parameters, dtype=float)
    if (observed.ndim != 2 or not observed.size
            or background.shape != observed.shape
            or pars.shape != (observed.shape[0], 3)
            or not np.isfinite(observed).all()
            or not np.isfinite(background).all()
            or np.any(observed <= 0) or np.any(background <= 0)):
        raise ValueError("Observed/background must be positive finite PS (M,N); parameters (M,3)")
    template = qwp_signal_power_template(
        signal_frequency_hz=frequency_support_hz,
        normalized_halo_weights=normalized_halo_weights,
        window_power_response=window_power_response,
        quadrature_power_fraction=quadrature_power_fraction,
        reference_coupling_gev_inv=reference_coupling_gev_inv,
        finesse_s=finesse_s, tau_s=tau_s, tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4, unit_conversion=unit_conversion)
    if template.shape != observed.shape:
        raise ValueError("Physical signal template must match observed PS (M,N)")
    peak_norm, peak_shape = factor_power_template(template)
    return dict(Y=np.log(observed), bkg=np.log(background),
                peak_norm=peak_norm, peak_shape=peak_shape,
                model_args=pars, signal_template_ps=template,
                reference_coupling_gev_inv=float(reference_coupling_gev_inv))


def log_likelihood_at_coupling(coupling_gev_inv, *, prepared, external_root):
    """Delegate exactly to the original vendor likelihood, for g >= 0."""
    g = float(coupling_gev_inv)
    if not np.isfinite(g) or g < 0:
        raise ValueError("Physical coupling must be finite nonnegative")
    g_ref = prepared["reference_coupling_gev_inv"]
    mu = (g/g_ref)**2
    return log_likelihood(
        mu, prepared["Y"], prepared["bkg"],
        prepared["peak_norm"], prepared["peak_shape"],
        prepared["model_args"], external_root)
