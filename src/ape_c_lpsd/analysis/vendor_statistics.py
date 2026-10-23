"""Thin read-only adapter to UNMODIFIED external C-LPSD statistical functions.

This provides the background spline + skew-normal fit and a generic PSD-signal
likelihood evaluation. It does not import LIGO response templates or define
the axion signal, significance, or upper limits.
"""
from __future__ import annotations
import hashlib
import sys
import warnings
from pathlib import Path

import numpy as np

FIT_SETTINGS = {
    "kernel_size": 30, "plot_mean": False, "verbose": False,
    "plot_path": None, "plot_prefix": "ape_background",
    "k_min": 4, "k_pruning": 1, "k_plateau": 3,
    "buffer": 50, "bin_scaling": 0.15, "cfd_fraction": 0.1,
    "disable_tqdm": True,
}


def load_external(external_root):
    """Locate source without editing, copying, or replacing external functions."""
    root = Path(external_root).expanduser().resolve()
    fit_file = root / "scalar_dark_matter_lpsd/lvk/fit_background.py"
    stats_file = root / "scalar_dark_matter_lpsd/lvk/stats.py"
    if not fit_file.is_file() or not stats_file.is_file():
        raise FileNotFoundError("Original external fit_background.py and stats.py required")
    for name, expected in (("scalar_dark_matter_lpsd.lvk.stats", stats_file),
                           ("scalar_dark_matter_lpsd.lvk.fit_background", fit_file)):
        cached = sys.modules.get(name)
        if cached is not None and Path(cached.__file__).resolve() != expected:
            raise RuntimeError("Wrong vendor module already imported: " + name)
    sys.path.insert(0, str(root))
    from scalar_dark_matter_lpsd.lvk import fit_background, stats
    if Path(stats.__file__).resolve() != stats_file or Path(fit_background.__file__).resolve() != fit_file:
        raise RuntimeError("Vendor import provenance mismatch")
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (fit_file, stats_file)}
    return fit_background, stats, hashes


def fit_window(frequency_hz, psd, external_root, *, settings=None):
    """Fit ONE local window; no alteration of vendor spline/skew-normal methods."""
    f = np.asarray(frequency_hz, dtype=float)
    p = np.asarray(psd, dtype=float)
    if (f.ndim != 1 or p.shape != f.shape or f.size < 1100
            or not np.isfinite(f).all() or not np.isfinite(p).all()
            or np.any(np.diff(f) <= 0) or np.any(f <= 0) or np.any(p <= 0)):
        raise ValueError("Need at least 1100 sorted positive frequencies and PSD values")
    external_fit, _, source_hashes = load_external(external_root)
    config = dict(FIT_SETTINGS if settings is None else settings)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        _, knots, pars, chi2 = external_fit.process_iteration(
            (0, np.log(f), np.log(p), config))
    if knots is None or pars is None or len(pars) != 3:
        raise RuntimeError("Unmodified background fitter returned no valid result")
    knots = np.asarray(knots, dtype=float)
    pars = np.asarray(pars, dtype=float)
    if (len(knots) % 2 or not np.isfinite(knots).all()
            or not np.isfinite(pars).all() or pars[1] <= 0
            or not np.isfinite(chi2)):
        raise ValueError("Invalid external fitting result")
    bg = np.exp(external_fit.bkg_model(np.log(f), knots))
    if bg.shape != p.shape or not np.isfinite(bg).all() or np.any(bg <= 0):
        raise ValueError("Invalid fitted background PSD")
    return {
        "background_psd": bg,
        "log_residual": np.log(p / bg),
        "knots": knots,
        "skewnormal_loc_scale_alpha": pars,
        "vendor_reduced_chi2": float(chi2),
        "vendor_quality_pass": bool(0.5 <= chi2 <= 10.0),
        "warnings": [str(x.message) for x in caught],
        "vendor_source_sha256": source_hashes,
        "settings": config,
    }


def log_likelihood(params, Y, bkg, peak_norm, peak_shape, model_args, external_root):
    """Delegate original vendor likelihood. Y/bkg LOG power (M,N)."""
    Y = np.asarray(Y, dtype=float)
    bkg = np.asarray(bkg, dtype=float)
    peak_norm = np.asarray(peak_norm, dtype=float)
    peak_shape = np.asarray(peak_shape, dtype=float)
    model_args = np.asarray(model_args, dtype=float)
    if Y.ndim != 2 or 0 in Y.shape or bkg.shape != Y.shape:
        raise ValueError("Y and bkg require equal nonempty (M,N)")
    m, n = Y.shape
    if (peak_shape.shape != (n,) or np.any(peak_shape < 0)
            or not np.isclose(peak_shape.sum(), 1., rtol=1e-12, atol=1e-12)):
        raise ValueError("peak_shape must be normalized nonnegative (N,)")
    if model_args.shape != (m, 3) or np.any(model_args[:, 1] <= 0):
        raise ValueError("model_args must be (M,3) with positive scale")
    try:
        np.broadcast_shapes(peak_norm.shape, Y.shape)
    except ValueError as exc:
        raise ValueError("peak_norm must broadcast to (M,N)") from exc
    if not all(np.isfinite(x).all() for x in (Y, bkg, peak_norm, peak_shape, model_args)):
        raise ValueError("Nonfinite inputs")
    if np.any(peak_norm < 0) or not np.isscalar(params) or not np.isfinite(params):
        raise ValueError("Invalid normalization or signal parameter")
    _, stats, _ = load_external(external_root)
    return float(stats.log_likelihood(
        float(params), Y, bkg, peak_norm, peak_shape, model_args))


def factor_power_template(signal_power):
    """Exactly factor reference signal PS (M,N) into vendor shape and norm."""
    power = np.asarray(signal_power, dtype=float)
    if power.ndim != 2 or 0 in power.shape or not np.isfinite(power).all() or np.any(power < 0):
        raise ValueError("Expected finite nonnegative signal PS (M,N)")
    sums = power.sum(axis=0)
    if sums.sum() <= 0:
        raise ValueError("Zero signal template")
    shape = sums / sums.sum()
    norm = np.zeros_like(power)
    np.divide(power, shape[None, :], out=norm, where=shape[None, :] > 0)
    return norm, shape
