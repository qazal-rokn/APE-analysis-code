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
    sys.path.insert(0, str(root))
    from scalar_dark_matter_lpsd.lvk import fit_background, stats
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


def log_likelihood(mu, observed_psd, background_psd, signal_psd_template,
                   skewnormal_parameters, external_root, *, weights=None):
    """Delegate to UNMODIFIED vendor stats.log_likelihood.

    signal_psd_template is REQUIRED. Physical APE axion response must be
    supplied and validated separately. No default scalar-DM/LIGO template.
    """
    data = np.asarray(observed_psd, dtype=float)
    bg = np.asarray(background_psd, dtype=float)
    template = np.asarray(signal_psd_template, dtype=float)
    pars = np.asarray(skewnormal_parameters, dtype=float)
    if data.ndim != 1 or len(data) == 0 or bg.shape != data.shape or template.shape != data.shape:
        raise ValueError("Observed/background/template must be matching 1D vectors")
    if not all(np.isfinite(v).all() for v in (data, bg, template)) or (
            np.any(data <= 0) or np.any(bg <= 0) or np.any(template < 0)):
        raise ValueError("Invalid PSD or template values")
    if pars.shape != (3,) or not np.isfinite(pars).all() or pars[1] <= 0:
        raise ValueError("Expected finite (loc, scale, alpha)")
    if not np.isfinite(mu) or mu < 0:
        raise ValueError("Signal amplitude mu must be nonnegative")
    w = np.ones(1) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != (1,) or not np.isfinite(w).all() or w[0] <= 0:
        raise ValueError("Expected one positive scalar weight")
    _, stats, _ = load_external(external_root)
    return float(stats.log_likelihood(
        float(mu), np.log(data)[None, :], np.log(bg)[None, :],
        template[None, :], w, pars[None, :]))
