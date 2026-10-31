"""UNTESTED: original-vendor multi-bin likelihood bridge for APE.

This module does NOT modify or reimplement vendor statistics. It accepts
already-prepared vendor-format multi-bin inputs. No production pipeline wiring,
scientific validation, or parity testing is claimed.
"""
from __future__ import annotations

import hashlib
import importlib
from pathlib import Path

import numpy as np

from ape_c_lpsd.analysis.vendor_statistics import load_external
from ape_c_lpsd.analysis.ape_vendor_fit import VENDOR_MLE_SHA256


def fit_vendor_multibin_untested(*, prepared, external_root, min_log10mu,
                                max_log10mu, n_seed_lkl):
    """Call the original vendor calc_max_lkl with a multi-bin peak shape.

    Inputs follow vendor make_args: Y and bkg are LOG power (M,N),
    peak_norm is (M,N), peak_shape is (N,) and sums to one,
    model_args is (M,3) in (loc, scale, alpha) order.

    Caller must construct the window, background and APE response separately.
    Vendor's signed-mu search, Fisher sigma, and error behavior are untouched.
    """
    root = Path(external_root).expanduser().resolve()
    vendor_path = root / "scalar_dark_matter_lpsd/lvk/max_lkl_executable.py"
    if not vendor_path.is_file():
        raise FileNotFoundError(vendor_path)
    digest = hashlib.sha256(vendor_path.read_bytes()).hexdigest()
    if digest != VENDOR_MLE_SHA256:
        raise RuntimeError("Original vendor optimizer SHA256 mismatch")

    Y = np.asarray(prepared["Y"], dtype=float)
    bkg = np.asarray(prepared["bkg"], dtype=float)
    norm = np.asarray(prepared["peak_norm"], dtype=float)
    shape = np.asarray(prepared["peak_shape"], dtype=float)
    params = np.asarray(prepared["model_args"], dtype=float)
    if Y.ndim != 2 or Y.shape[0] == 0 or Y.shape[1] < 2:
        raise ValueError("Require nonempty vendor multi-bin Y (M,N), N >= 2")
    m, n_bins = Y.shape
    if (bkg.shape != Y.shape or norm.shape != Y.shape or
            shape.shape != (n_bins,) or params.shape != (m, 3)):
        raise ValueError("Vendor likelihood input shapes do not match")
    if not all(np.isfinite(v).all() for v in (Y, bkg, norm, shape, params)):
        raise ValueError("Nonfinite vendor likelihood input")
    if (np.any(norm < 0) or np.any(shape < 0) or
            not np.isclose(shape.sum(), 1.0, rtol=1e-12, atol=1e-12)
            or np.any(params[:, 1] <= 0)):
        raise ValueError("Invalid vendor normalization or skew-normal scale")
    a, b, seeds = float(min_log10mu), float(max_log10mu), int(n_seed_lkl)
    if not np.isfinite([a, b]).all() or a >= b or seeds < 10 or seeds != n_seed_lkl:
        raise ValueError("Invalid vendor optimizer settings")

    load_external(root)
    module = importlib.import_module("scalar_dark_matter_lpsd.lvk.max_lkl_executable")
    if Path(module.__file__).resolve() != vendor_path:
        raise RuntimeError("Wrong vendor optimizer module imported")

    # The entire vendor likelihood and optimizer run unchanged.
    l0, lmax, muhat, sigma = module.calc_max_lkl(
        Y, bkg, params, norm, shape,
        min_log10mu=a, max_log10mu=b, n_seed_lkl=seeds)
    return {
        "status": "UNTESTED_NOT_VALIDATED_NOT_FOR_SCIENTIFIC_RESULTS",
        "log_likelihood_null": float(l0),
        "log_likelihood_max": float(lmax),
        "mu_hat": float(muhat),
        "sigma_mu": float(sigma),
        "vendor_mle_sha256": digest,
        "n_bins": n_bins,
        "n_records": m,
    }
