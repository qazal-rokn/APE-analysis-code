"""Use original vendor max-likelihood optimization and original expected UL.

This module is an interface, not an APE physical sensitivity/coverage claim.
Pass already-prepared calibrated APE PS likelihood inputs, not raw detector data.
The vendor optimizer is called unchanged, including its signed-mu treatment.
"""
from __future__ import annotations
import hashlib
import importlib
from pathlib import Path
import numpy as np

from ape_c_lpsd.analysis.vendor_statistics import load_external
from ape_c_lpsd.analysis.ape_expected_upper_limit import expected_upper_limit_ape

VENDOR_MLE_SHA256 = '031ad21dde68b08abb33efa18ac535ff9e7ea00dbc03cf13b6ec2f479328967b'


def fit_vendor_mle_and_expected_ul(*, prepared, external_root, min_log10mu,
                                   max_log10mu, n_seed_lkl, alpha_cl=0.95):
    """Return vendor (L0, Lmax, muhat, sigma), one-sided q0 and expected UL.

    This pilot supports the original vendor's DEFAULT single-bin peak shape.
    Multi-bin validity-range behavior requires separate review; do not silently
    treat the one-bin optimizer validation as a multi-bin halo validation.
    """
    root = Path(external_root).expanduser().resolve()
    vendor = root / 'scalar_dark_matter_lpsd/lvk/max_lkl_executable.py'
    if not vendor.is_file():
        raise FileNotFoundError(vendor)
    digest = hashlib.sha256(vendor.read_bytes()).hexdigest()
    if digest != VENDOR_MLE_SHA256:
        raise RuntimeError('Original vendor max_lkl_executable.py SHA256 mismatch')

    Y = np.asarray(prepared['Y'], dtype=float)
    bkg = np.asarray(prepared['bkg'], dtype=float)
    norm = np.asarray(prepared['peak_norm'], dtype=float)
    shape = np.asarray(prepared['peak_shape'], dtype=float)
    params = np.asarray(prepared['model_args'], dtype=float)
    g_ref = float(prepared['reference_coupling_gev_inv'])
    a, b, n = float(min_log10mu), float(max_log10mu), int(n_seed_lkl)
    if (Y.ndim != 2 or Y.shape[1] != 1 or bkg.shape != Y.shape or
            shape.shape != (1,) or not np.array_equal(shape, [1.0]) or
            norm.shape != Y.shape or params.shape != (len(Y), 3) or
            Y.size == 0 or not all(np.isfinite(v).all() for v in (Y, bkg, norm, shape, params)) or
            np.any(norm < 0) or np.any(params[:, 1] <= 0) or not np.isfinite(g_ref) or g_ref <= 0):
        raise ValueError('Require finite single-bin vendor likelihood arrays and g_ref')
    if (not np.isfinite(a) or not np.isfinite(b) or
            not a < b or n < 10 or n != n_seed_lkl):
        raise ValueError('Invalid explicitly supplied vendor optimizer scan settings')

    # Existing APE vendor loader checks source provenance and dependencies.
    load_external(root)
    module = importlib.import_module('scalar_dark_matter_lpsd.lvk.max_lkl_executable')
    if Path(module.__file__).resolve() != vendor:
        raise RuntimeError('Wrong vendor optimizer module imported')

    l0, lmax, muhat, sigma = module.calc_max_lkl(
        Y, bkg, params, norm, shape,
        min_log10mu=a, max_log10mu=b, n_seed_lkl=n)
    if not np.isfinite([l0, lmax, muhat, sigma]).all() or sigma <= 0:
        raise ValueError('Vendor returned invalid fitted maximum/uncertainty')
    if lmax < l0 - 1e-8:
        raise ValueError('Vendor optimum is below null likelihood')

    # Identical one-sided q0 definition used by the original vendor search.
    q0 = max(0.0, 2.0 * (lmax - l0)) if muhat > 0 else 0.0
    limit = expected_upper_limit_ape(
        sigma_mu=float(sigma), g_ref_gev_inv=g_ref,
        external_root=root, alpha_cl=alpha_cl)
    return dict(
        status='provisional_vendor_single_bin_software_fit_not_ape_coverage_validated',
        log_likelihood_null=float(l0), log_likelihood_max=float(lmax),
        mu_hat=float(muhat), sigma_mu=float(sigma), q0=float(q0),
        local_sqrt_q0=float(np.sqrt(q0)),
        expected_upper_limit=limit,
        vendor_mle_sha256=digest,
        optimizer=dict(min_log10mu=a, max_log10mu=b, n_seed_lkl=n),
    )
