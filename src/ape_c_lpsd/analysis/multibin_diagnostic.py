"""Opt-in synthetic-only multi-bin APE template diagnostic.

Input NPZ arrays are externally supplied, not synthesized here:
  frequency_support_hz (K,), normalized_halo_weights (K,),
  observed_frequency_windows_hz (C,N), window_power_response (C,2,N,K).
The window matrix is *power response*, not a probability matrix.
No C output is modified and no additional K or other normalization is applied.
"""
import hashlib
from pathlib import Path

import numpy as np

from ape_c_lpsd.analysis.ape_likelihood import (
    prepare_ape_likelihood, log_likelihood_at_coupling,
)


def _sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def evaluate_multibin_windows(*, frequency_hz, observed_ps, fitted_backgrounds,
                              candidate_start_indices, config, external_vendor_root):
    """Evaluate unchanged vendor log likelihood at explicit nonnegative g values.

    Returns only a diagnostic likelihood scan, *not* significance, MLE, or UL.
    All spectra are native-C PS in rad^2. No interpolation is performed.
    """
    source = Path(config["multibin_template_file"]).expanduser().resolve(strict=True)
    if source.suffix != ".npz":
        raise ValueError("Multi-bin input must be a .npz file")
    for description in ("halo_source_description", "window_response_source_description"):
        if not isinstance(config.get(description), str) or not config[description].strip():
            raise ValueError("Explicit physical/source provenance required: " + description)
    couplings = np.asarray(config["coupling_grid_gev_inv"], dtype=float)
    if (couplings.ndim != 1 or len(couplings) < 2 or
            not np.isfinite(couplings).all() or couplings[0] != 0 or
            np.any(np.diff(couplings) <= 0)):
        raise ValueError("Provide strictly increasing finite coupling grid starting at zero")
    if config.get("n_effective_trials") is not None:
        raise ValueError("No global trials inference authorized in multi-bin diagnostic")
    starts = np.asarray(candidate_start_indices, dtype=int)
    frequency = np.asarray(frequency_hz, dtype=float)
    observed = np.asarray(observed_ps, dtype=float)
    if (starts.ndim != 1 or not starts.size or np.any(np.diff(starts) <= 0)
            or np.any(starts < 0) or observed.shape != (2, len(frequency))):
        raise ValueError("Invalid candidate selection or X/Y C spectra")

    expected_keys = {
        "frequency_support_hz", "normalized_halo_weights",
        "observed_frequency_windows_hz", "window_power_response",
    }
    with np.load(source, allow_pickle=False) as data:
        if set(data.files) != expected_keys:
            raise ValueError("Require exactly the four declared physical template arrays")
        support = np.asarray(data["frequency_support_hz"], dtype=float)
        weights = np.asarray(data["normalized_halo_weights"], dtype=float)
        window_freqs = np.asarray(data["observed_frequency_windows_hz"], dtype=float)
        W = np.asarray(data["window_power_response"], dtype=float)
    if (support.ndim != 1 or not len(support) or
            not np.isfinite(support).all() or np.any(support <= 0)
            or np.any(np.diff(support) <= 0) or weights.shape != support.shape
            or not np.isfinite(weights).all() or np.any(weights < 0)
            or not np.isclose(weights.sum(), 1, atol=1e-12, rtol=1e-12)):
        raise ValueError("Invalid external intrinsic-field support or halo probability masses")
    if (window_freqs.ndim != 2 or window_freqs.shape[0] != len(starts)
            or window_freqs.shape[1] < 2 or
            W.shape != (len(starts), 2, window_freqs.shape[1], len(support))
            or not np.isfinite(W).all() or np.any(W < 0)):
        raise ValueError("Window power response shape/values are invalid")
    n = window_freqs.shape[1]
    if np.any(starts + n > len(frequency)):
        raise ValueError("Candidate window exceeds native C frequency grid")
    if not np.array_equal(window_freqs, np.stack([frequency[k:k+n] for k in starts])):
        raise ValueError("Window frequencies must match actual native C output exactly")
    physics = config["physics"]
    results = []
    for i, start in enumerate(starts):
        stop = int(start + n)
        prepared = prepare_ape_likelihood(
            observed_ps=observed[:, start:stop],
            background_ps=np.stack([
                np.asarray(fit["background_psd"])[start:stop]
                for fit in fitted_backgrounds
            ]),
            skewnormal_parameters=np.stack([
                fit["skewnormal_loc_scale_alpha"] for fit in fitted_backgrounds
            ]),
            frequency_support_hz=support,
            normalized_halo_weights=weights,
            window_power_response=W[i],
            quadrature_power_fraction=physics["quadrature_power_fraction"],
            reference_coupling_gev_inv=physics["reference_coupling_gev_inv"],
            finesse_s=physics["finesse_s"], tau_s=physics["tau_s"],
            tau1_s=physics["tau1_s"], rho_dm_gev4=physics["rho_dm_gev4"],
            unit_conversion=physics["unit_conversion"],
        )
        power = prepared["signal_template_ps"]
        if not np.isfinite(power).all() or np.any(power < 0):
            raise ValueError("Invalid physical signal power")
        if not np.allclose(
                prepared["peak_norm"] * prepared["peak_shape"][None, :],
                power, atol=0, rtol=5e-14):
            raise RuntimeError("Vendor factorization does not reproduce signal PS")
        log_values = [
            log_likelihood_at_coupling(
                float(g), prepared=prepared, external_root=external_vendor_root)
            for g in couplings
        ]
        if not np.isfinite(log_values).all():
            raise ValueError("Vendor likelihood returned nonfinite values")
        results.append({
            "candidate_start_index": int(start),
            "frequencies_hz": window_freqs[i].tolist(),
            "signal_template_ps_rad2": power.tolist(),
            "log_likelihood_by_coupling": [float(x) for x in log_values],
        })
    return {
        "status": "provisional_multibin_forward_likelihood_diagnostic_only",
        "coupling_grid_gev_inv": couplings.tolist(),
        "windows": results,
        "template_source": str(source),
        "template_sha256": _sha256(source),
        "halo_source_description": config["halo_source_description"],
        "window_response_source_description": config["window_response_source_description"],
        "native_c_power_units": "rad^2",
        "additional_normalization_factor": None,
        "candidate_significance": None,
        "maximum_likelihood_fit": None,
        "expected_upper_limit": None,
        "physical_coverage_validated": False,
        "overlapping_bins_independence_validated": False,
    }
