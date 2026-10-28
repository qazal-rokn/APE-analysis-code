"""Thin APE adapter for UNMODIFIED external 2026 vendor expected limits.

The vendor's get_upper_limits_approx accepts fitted sigma_mu and returns
sqrt(mu_UL) and a propagated uncertainty. This is an asymptotic EXPECTED
upper limit, not an observed profile likelihood limit or coverage proof.
"""
from __future__ import annotations
import hashlib
from pathlib import Path
import numpy as np

VENDOR_UL_SHA256 = "8f800fed30a0327763c87f983df7f507b44a879aa6ccfa8f57558146e1d26745"


def expected_upper_limit_ape(*, sigma_mu, g_ref_gev_inv, external_root,
                             alpha_cl=0.95, tol=1e-10):
    sigma = float(sigma_mu)
    g_ref = float(g_ref_gev_inv)
    alpha = float(alpha_cl)
    tol = float(tol)
    if not all(np.isfinite(x) for x in (sigma, g_ref, alpha, tol)):
        raise ValueError("Expected finite statistical inputs")
    if sigma <= 0 or g_ref <= 0 or not 0 < alpha < 1 or tol <= 0:
        raise ValueError("Invalid sigma_mu, g_ref, confidence level or tolerance")
    root = Path(external_root).expanduser().resolve()
    target = root / "scalar_dark_matter_lpsd/lvk/calc_upper_limit.py"
    if not target.is_file():
        raise FileNotFoundError(target)
    checksum = hashlib.sha256(target.read_bytes()).hexdigest()
    if checksum != VENDOR_UL_SHA256:
        raise RuntimeError("Unexpected vendor upper limit source version")

    # Reuse the existing guarded vendor loader and verify module provenance.
    from ape_c_lpsd.analysis.vendor_statistics import load_external
    load_external(root)
    import importlib
    source = importlib.import_module("scalar_dark_matter_lpsd.lvk.calc_upper_limit")
    if Path(source.__file__).resolve() != target:
        raise RuntimeError("Wrong vendor upper limit module imported")

    sqrt_mu_ul, sigma_sqrt_mu = source.get_upper_limits_approx(
        sigma, alpha_CL=alpha, tol=tol)
    limit = np.asarray(sqrt_mu_ul, dtype=float).reshape(-1)
    uncertainty = np.asarray(sigma_sqrt_mu, dtype=float).reshape(-1)
    if limit.size != 1 or uncertainty.size != 1:
        raise RuntimeError("Unexpected vendor upper limit output shape")
    if not np.isfinite(limit[0]) or limit[0] <= 0:
        raise RuntimeError("Invalid vendor upper limit result")
    if not np.isfinite(uncertainty[0]) or uncertainty[0] < 0:
        raise RuntimeError("Invalid vendor upper limit uncertainty")
    return {
        "status": "provisional_asymptotic_expected_limit_not_ape_coverage_validated",
        "confidence_level": alpha,
        "sigma_mu": sigma,
        "sqrt_mu_ul": float(limit[0]),
        "sigma_sqrt_mu": float(uncertainty[0]),
        "g_ref_gev_inv": g_ref,
        "g_ul_gev_inv": float(g_ref * limit[0]),
        "sigma_g_ul_gev_inv": float(g_ref * uncertainty[0]),
        "vendor_upper_limit_sha256": checksum,
        "method": "unchanged_vendor_get_upper_limits_approx",
    }
