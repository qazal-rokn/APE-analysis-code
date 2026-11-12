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
