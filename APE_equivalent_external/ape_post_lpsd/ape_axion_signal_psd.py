"""APE-specific deterministic QWP/halo response on the native C-LPSD D (PSD) grid.

Preserves the existing provisional APE QWP amplitude formula verbatim by import.
No LIGO response, vendor statistics, inferred axion phase, or guessed ENBW.
Inputs include an externally validated C *PSD* kernel (per unit signal variance).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling


class APEModelError(ValueError):
    """An APE optical or native C response input fails the model contract."""


def build_response_psd(
    *, native_frequency_hz, output_bin_indices, halo_frequency_hz,
    normalized_halo_weights, native_c_psd_kernel_per_variance_hz,
    quadrature_power_fraction, reference_coupling_gev_inv,
    finesse_s, tau_s, tau1_s, rho_dm_gev4, unit_conversion,
):
    """Return full-grid X/Y PSD per unit mu=(g/g_ref)^2, without interpolation.

    Native C kernel K[m,j,k] is the output *D* [1/Hz] for a unit-variance
    sinusoid at halo frequency k, with original APE C window/averaging flags.
    It is NOT the unit-peak PS kernel used in earlier APE demonstrations.

    T[m,j] = 0.5*g_ref^2*a[m]*sum_k h[k]*R(f_k)^2*K[m,j,k].
    R is the existing provisional peak amplitude per axion coupling.
    """
    native = np.asarray(native_frequency_hz, dtype=np.float64)
    indices = np.asarray(output_bin_indices)
    halo = np.asarray(halo_frequency_hz, dtype=np.float64)
    weights = np.asarray(normalized_halo_weights, dtype=np.float64)
    kernel = np.asarray(native_c_psd_kernel_per_variance_hz, dtype=np.float64)
    fractions = np.asarray(quadrature_power_fraction, dtype=np.float64)
    if (native.ndim != 1 or native.size < 2 or not np.isfinite(native).all()
            or np.any(native <= 0) or np.any(np.diff(native) <= 0)):
        raise APEModelError('Native C frequency grid must be finite, positive, strictly increasing')
    if (indices.ndim != 1 or indices.size == 0 or indices.dtype.kind not in 'iu'
            or np.any(indices < 0) or np.any(indices >= native.size)
            or np.any(np.diff(indices) <= 0)):
        raise APEModelError('Selected output-bin indices must be unique, increasing native-grid indices')
    if (halo.ndim != 1 or halo.size == 0 or not np.isfinite(halo).all()
            or np.any(halo <= 0) or np.any(np.diff(halo) <= 0)):
        raise APEModelError('Halo support frequencies must be increasing, positive and finite')
    if (weights.shape != halo.shape or not np.isfinite(weights).all()
            or np.any(weights < 0) or not np.isclose(weights.sum(), 1., atol=1e-12, rtol=1e-12)):
        raise APEModelError('Supplied halo probability weights must be nonnegative and sum to one')
    if (kernel.shape != (2, indices.size, halo.size)
            or not np.isfinite(kernel).all() or np.any(kernel < 0)):
        raise APEModelError('Native C PSD kernel must have shape (2, selected bins, halo modes), nonnegative')
    if (fractions.shape != (2,) or not np.isfinite(fractions).all()
            or np.any(fractions < 0) or np.any(fractions > 1) or not np.any(fractions > 0)):
        raise APEModelError('Explicit X/Y power fractions must each be in [0,1] and not both zero')
    gref = float(reference_coupling_gev_inv)
    if not np.isfinite(gref) or gref <= 0:
        raise APEModelError('Reference axion coupling must be positive and finite')
    response = qwp_amplitude_per_coupling(
        halo, finesse_s=finesse_s, tau_s=tau_s, tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4, unit_conversion=unit_conversion)
    # The *existing* physical response is an amplitude. Convert peak-square
    # into time-averaged variance via /2; K already supplies PSD per variance.
    selected = (0.5*gref**2)*fractions[:, None]*np.einsum(
        'k,mjk->mj', weights*response**2, kernel, optimize=True)
    if not np.isfinite(selected).all() or np.any(selected < 0) or not np.any(selected > 0):
        raise APEModelError('Model PSD is zero, nonfinite or negative')
    result = np.zeros((2, native.size), dtype=np.float64)
    result[:, indices] = selected
    return result


def save_provisional_response(
    path: str | Path, *, native_frequency_hz, output_bin_indices,
    halo_frequency_hz, normalized_halo_weights,
    native_c_psd_kernel_per_variance_hz, quadrature_power_fraction,
    reference_coupling_gev_inv, finesse_s, tau_s, tau1_s,
    rho_dm_gev4, unit_conversion, halo_source: str,
    native_c_kernel_source: str, calibration_version: str,
) -> dict:
    """Write Step-3-compatible NPZ with explicit non-approved science status.

    This file is not a physical axion upper-limit calibration. The caller must
    supply the halo shape, C PSD normalization and quadrature response.
    """
    for name, value in (('halo_source', halo_source),
                        ('native_c_kernel_source', native_c_kernel_source),
                        ('calibration_version', calibration_version)):
        if not isinstance(value, str) or not value.strip():
            raise APEModelError(f'Missing provenance: {name}')
    template = build_response_psd(
        native_frequency_hz=native_frequency_hz,
        output_bin_indices=output_bin_indices,
        halo_frequency_hz=halo_frequency_hz,
        normalized_halo_weights=normalized_halo_weights,
        native_c_psd_kernel_per_variance_hz=native_c_psd_kernel_per_variance_hz,
        quadrature_power_fraction=quadrature_power_fraction,
        reference_coupling_gev_inv=reference_coupling_gev_inv,
        finesse_s=finesse_s, tau_s=tau_s, tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4, unit_conversion=unit_conversion)
    output = Path(path).expanduser().resolve()
    if output.suffix != '.npz' or output.exists():
        raise APEModelError('Output must be a new .npz file (no overwrites)')
    metadata = {
        'mu_definition': 'mu=(g_a_gamma_gamma/reference_coupling_gev_inv)^2',
        'reference_coupling_gev_inv': float(reference_coupling_gev_inv),
        'response_units': 'rad^2/Hz per unit mu',
        'response_source': 'provisional APE QWP Case-2 from development/likelihood/src/physics/qwp_signal.py',
        'source_qwp_code_sha': '7dcee0f86c8bc314cd46c525b37724ee7c8cb189',
        'calibration_version': calibration_version,
        'spectral_shape_convention': 'halo and native C window already in peak_norm; use peak_shape=[1] for current single-bin bridge',
        'halo_source': halo_source,
        'native_c_kernel_source': native_c_kernel_source,
        'native_c_output_column': 'FORMAT fDNbri: D = one-sided PSD',
        'native_c_kernel_units': 'Hz^-1 per unit time-domain signal variance',
        'optical_response_status': 'PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS',
        'halo_and_kernel_validation_status': 'CALLER_SUPPLIED_NOT_CERTIFIED_BY_THIS_MODULE',
        'random_DM_phase': False,
        'random_DM_amplitude': False,
        'stochasticity_of_DM_signal_neglected': True,
        'quadratures_combined_in_likelihood': False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output.open('xb') as stream:
            np.savez(stream, frequency_hz=np.asarray(native_frequency_hz, dtype=np.float64),
                     peak_norm_x=template[0], peak_norm_y=template[1],
                     metadata_json=json.dumps(metadata, sort_keys=True))
    except BaseException:
        output.unlink(missing_ok=True)
        raise
    return metadata


def inject_expected_psd(background_psd, signal_psd_per_mu, *, mu: float):
    """Deterministic POWER-domain injection; no random DM phases or C execution."""
    b = np.asarray(background_psd, dtype=np.float64)
    s = np.asarray(signal_psd_per_mu, dtype=np.float64)
    if (b.shape != s.shape or not b.size or not np.isfinite(b).all()
            or not np.isfinite(s).all() or np.any(b <= 0) or np.any(s < 0)
            or not np.isfinite(mu) or mu < 0):
        raise APEModelError('Invalid background, signal template or nonnegative coupling-squared mu')
    out = b + float(mu)*s
    if not np.isfinite(out).all() or np.any(out <= 0):
        raise APEModelError('Invalid injected PSD')
    return out
