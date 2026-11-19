"""APE-only signal-power response boundary for the unchanged vendor likelihood.

No optical transfer function or axion coupling convention is assumed here.
The input must already be a scientifically approved expected PSD *per unit mu*.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


class APEResponseError(ValueError):
    """An APE response cannot safely be supplied to the vendor likelihood."""


def load_peak_norm(path: str | Path, native_frequency_hz: np.ndarray) -> tuple[np.ndarray, dict]:
    """Return vendor peak_norm[2, N] for X/Y, with no interpolation or rescaling.

    NPZ contract:
        frequency_hz: native positive monotonically increasing C-LPSD grid
        peak_norm_x: additive PSD for X *per unit mu* [rad^2/Hz]
        peak_norm_y: additive PSD for Y *per unit mu* [rad^2/Hz]
        metadata_json: JSON scalar string with approved normalization provenance.

    The physical meaning of mu is defined by the supplied metadata, NOT by
    this adapter. In particular, g_{a gamma gamma} conversion is NOT assumed.
    """
    grid = np.asarray(native_frequency_hz, dtype=np.float64)
    if grid.ndim != 1 or grid.size == 0 or not np.isfinite(grid).all() or np.any(grid <= 0) or np.any(np.diff(grid) <= 0):
        raise APEResponseError("Native frequency grid must be finite, positive, 1D and strictly increasing")
    try:
        with np.load(Path(path), allow_pickle=False) as data:
            if set(data.files) != {"frequency_hz", "peak_norm_x", "peak_norm_y", "metadata_json"}:
                raise APEResponseError("Response NPZ requires exactly frequency_hz, peak_norm_x, peak_norm_y, metadata_json")
            frequencies = np.asarray(data["frequency_hz"], dtype=np.float64)
            x = np.asarray(data["peak_norm_x"], dtype=np.float64)
            y = np.asarray(data["peak_norm_y"], dtype=np.float64)
            raw = data["metadata_json"]
            if raw.ndim != 0:
                raise APEResponseError("metadata_json must be a scalar JSON string")
            metadata = json.loads(str(raw.item()))
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise APEResponseError(f"Cannot read APE response: {exc}") from exc
    if frequencies.shape != grid.shape or not np.array_equal(frequencies, grid):
        raise APEResponseError("Response frequency grid differs from native C-LPSD frequency grid")
    if any(a.shape != grid.shape or not np.isfinite(a).all() or np.any(a < 0) for a in (x, y)):
        raise APEResponseError("APE X/Y PSD responses must be finite and nonnegative at every native bin")
    if not np.any(x > 0) and not np.any(y > 0):
        raise APEResponseError("Both APE PSD responses are identically zero")
    if not isinstance(metadata, dict) or metadata.get("response_units") != "rad^2/Hz per unit mu":
        raise APEResponseError("Response units are not declared as rad^2/Hz per unit mu")
    for field in ("mu_definition", "response_source", "calibration_version", "spectral_shape_convention"):
        if not isinstance(metadata.get(field), str) or not metadata[field].strip():
            raise APEResponseError(f"Missing approved response metadata: {field}")
    return np.stack((x, y)), metadata


def vendor_likelihood_inputs(
    log_psd_xy: np.ndarray,
    log_background_xy: np.ndarray,
    skewnormal_parameters_xy: np.ndarray,
    peak_norm_xy: np.ndarray,
    peak_shape: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Validate and pass arrays through unchanged to vendor stats.log_likelihood.

    In the original model, peak_norm[:, j] * peak_shape[j] * mu is
    an *additive one-sided PSD*. Do not put a peak shape into the supplied
    peak_norm twice. Do not sum X/Y here or modify the skew-normal model.
    """
    y = np.asarray(log_psd_xy, dtype=np.float64)
    b = np.asarray(log_background_xy, dtype=np.float64)
    pars = np.asarray(skewnormal_parameters_xy, dtype=np.float64)
    response = np.asarray(peak_norm_xy, dtype=np.float64)
    shape = np.asarray(peak_shape, dtype=np.float64)
    if y.ndim != 2 or y.shape[0] != 2 or y.shape[1] == 0:
        raise APEResponseError("APE likelihood data must have independent X/Y rows")
    if b.shape != y.shape or response.shape != y.shape or pars.shape != (2, 3) or shape.shape != (y.shape[1],):
        raise APEResponseError("Original likelihood input array shapes do not agree")
    if not all(np.isfinite(a).all() for a in (y, b, pars, response, shape)):
        raise APEResponseError("Nonfinite likelihood input")
    if np.any(response < 0) or np.any(shape < 0) or not np.isclose(shape.sum(), 1, atol=1e-14, rtol=0):
        raise APEResponseError("Response and normalized peak shape must be nonnegative")
    if np.any(pars[:, 1] <= 0):
        raise APEResponseError("Skewnormal scales must be positive")
    return y, b, response, shape, pars


# Previously used APE calibration constants, preserved as provenance only:
# the calibrated input HDF5 already contains rho_x/rho_y in radians.
# These numbers MUST NOT be reapplied to its PSD.
APE_MEASURED_CALIBRATION = {
    "cavity_transmission_w_per_v": 0.0284,
    "responsivity_a_per_w": 0.7,
    "second_harmonic_pd_response": 0.7,
    "transimpedance_gain_v_per_a": 1000000.0,
    "input_calibration_method": "ape_rms_quadrature_calibration_v1",
}
# Original APE *provisional working values*, not verified cavity measurements.
APE_QWP_WORKING_VALUES = {
    "finesse_s": 1500.0,
    "tau_s": 2.5e-9,
    "tau1_s": 1.0e-8,
    "rho_dm_gev4": 3.1e-42,
    "unit_conversion": 1.519e24,
}
APE_REFERENCE_COUPLING_GEV_INV = 1e-10
# Retain original LIGO get_beta() linewidth convention, not a new halo model.
ORIGINAL_VENDOR_FRACTIONAL_LINEWIDTH = 1e-6


def make_ape_q_peak_norm(frequency_hz, *,
                         finesse_s=1500.0, tau_s=2.5e-9,
                         tau1_s=1.0e-8, rho_dm_gev4=3.1e-42,
                         unit_conversion=1.519e24,
                         reference_coupling_gev_inv=1e-10):
    """APE replacement for LIGO get_beta: original single-bin linewidth.

    Assume the already calibrated X,Y readout measures orthogonal components
    of the same angle: theta^2 = theta_x^2 + theta_y^2. Only the SUM Q is
    calibrated by this routine. R is existing nominal peak angular response
    per GeV^-1. Peak-to-variance factor 1/2; original scalar get_beta uses
    1 / (pi*f*10^-6) in its peak PSD. No additional halo or C window kernel.
    """
    from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling
    f = np.asarray(frequency_hz, dtype=np.float64)
    if (f.ndim != 1 or f.size == 0 or not np.isfinite(f).all()
            or np.any(f <= 0) or np.any(np.diff(f) <= 0)):
        raise APEResponseError("Need full positive increasing native frequency grid")
    gref = float(reference_coupling_gev_inv)
    if not np.isfinite(gref) or gref <= 0:
        raise APEResponseError("Reference coupling must be finite and positive")
    response = qwp_amplitude_per_coupling(
        f, finesse_s=finesse_s, tau_s=tau_s, tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4, unit_conversion=unit_conversion)
    # Original single-bin statistical signal-shape/linewidth convention.
    beta_q = 0.5 * np.square(gref * response) / (
        np.pi * f * ORIGINAL_VENDOR_FRACTIONAL_LINEWIDTH
    )
    if not np.isfinite(beta_q).all() or np.any(beta_q <= 0):
        raise APEResponseError("APE Q response is nonfinite or nonpositive")
    return beta_q


def save_ape_q_response(spectrum_h5, output_npz, *,
                        finesse_s=1500.0, tau_s=2.5e-9,
                        tau1_s=1.0e-8, rho_dm_gev4=3.1e-42,
                        unit_conversion=1.519e24,
                        reference_coupling_gev_inv=1e-10):
    """Save APE-only Q beta on original C frequency grid; do not rescale PSD."""
    from scalar_dark_matter_lpsd.LPSDIO import LPSDOutput
    source = Path(spectrum_h5).expanduser().resolve()
    dst = Path(output_npz).expanduser().resolve()
    if not source.name.endswith("_Q.h5") or not source.is_file():
        raise APEResponseError("Need existing combined Q HDF5 spectrum")
    if dst.suffix != ".npz" or dst.exists() or not dst.parent.is_dir():
        raise APEResponseError("Need a new .npz path in an existing folder")
    native = LPSDOutput(str(source)).freq
    beta_q = make_ape_q_peak_norm(
        native, finesse_s=finesse_s, tau_s=tau_s, tau1_s=tau1_s,
        rho_dm_gev4=rho_dm_gev4, unit_conversion=unit_conversion,
        reference_coupling_gev_inv=reference_coupling_gev_inv)
    metadata = {
        "configuration": "A_Q",
        "mu_definition": "mu=(g_a_gamma_gamma/reference_coupling_gev_inv)^2",
        "reference_coupling_gev_inv": reference_coupling_gev_inv,
        "response_units": "rad^2/Hz per unit mu",
        "response_source": "existing ape_post_lpsd.physics.qwp_signal.qwp_amplitude_per_coupling",
        "calibration_version": "ape_rms_quadrature_calibration_v1",
        "measured_input_calibration_provenance_only": APE_MEASURED_CALIBRATION,
        "ape_qwp_parameters": {
            "finesse_s": finesse_s, "tau_s": tau_s,
            "tau1_s": tau1_s, "rho_dm_gev4": rho_dm_gev4,
            "unit_conversion": unit_conversion,
        },
        "spectral_shape_convention": "original-vendor single-bin peak_shape=[1], fractional linewidth=1e-6",
        "quadrature_combination": "S_Q=S_X+S_Y; assumes orthogonal, equally calibrated angle quadratures",
        "response_status": "PROVISIONAL_QWP_WORKING_PARAMETERS_NOT_MEASURED",
        "physical_upper_limits_authorized": False,
    }
    try:
        with dst.open("xb") as stream:
            np.savez(stream, frequency_hz=native, peak_norm_q=beta_q,
                     metadata_json=json.dumps(metadata, sort_keys=True))
    except BaseException:
        dst.unlink(missing_ok=True)
        raise
    return metadata


def load_peak_norm_q(path, native_frequency_hz):
    """Load the combined Q response with no interpolation or extra scaling."""
    f = np.asarray(native_frequency_hz, dtype=np.float64)
    try:
        with np.load(path, allow_pickle=False) as data:
            if set(data.files) != {"frequency_hz", "peak_norm_q", "metadata_json"}:
                raise APEResponseError("Combined Q response requires exactly three NPZ fields")
            recorded = np.asarray(data["frequency_hz"], dtype=np.float64)
            q = np.asarray(data["peak_norm_q"], dtype=np.float64)
            metadata = json.loads(str(data["metadata_json"].item()))
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise APEResponseError(f"Cannot load combined Q response: {exc}") from exc
    if not (recorded.shape == f.shape and np.array_equal(recorded, f)):
        raise APEResponseError("APE Q response frequency grid mismatches original C-LPSD")
    if q.shape != f.shape or not np.isfinite(q).all() or np.any(q <= 0):
        raise APEResponseError("APE Q response must be finite and positive")
    if not isinstance(metadata, dict) or metadata.get("response_units") != "rad^2/Hz per unit mu":
        raise APEResponseError("Incorrect APE Q response units")
    return q, metadata
