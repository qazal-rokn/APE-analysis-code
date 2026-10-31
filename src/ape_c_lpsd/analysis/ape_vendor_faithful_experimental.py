"""EXPERIMENTAL vendor-faithful APE post-C-LPSD analysis.

Runs the ORIGINAL vendor background fit on native APE spectra converted to
vendor LPSDOutput containers. The vendor source is never modified. This
stage is not a scientifically validated end-to-end pipeline.

The APE power response must be injected downstream at the vendor's peak_norm
boundary; no LIGO transfer function is applied here.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ape_c_lpsd.analysis.vendor_statistics import load_external


STATUS = "EXPERIMENTAL_UNTESTED_VENDOR_FAITHFUL_NOT_PRODUCTION"


def vendor_background_from_ape_lpsd(*, ape_lpsd_file, output_json,
                                    external_root, segment_size=1000,
                                    buffer=50, pruning=1, ana_fmin=10.,
                                    ana_fmax=5000., n_processes=1,
                                    **vendor_options):
    """Call the vendor's COMPLETE fit_background_in_file unchanged.

    ape_lpsd_file must already be an original-vendor-readable LPSDOutput
    file. It is NOT a generic two-column APE spectrum. The exact vendor
    segmentation, buffering, fit, JSON serialization and failure handling
    are delegated to the vendor.
    """
    src = Path(ape_lpsd_file).expanduser().resolve(strict=True)
    root = Path(external_root).expanduser().resolve(strict=True)
    dest = Path(output_json).expanduser().resolve()
    if dest.exists():
        raise FileExistsError("Refusing to overwrite existing vendor fit JSON")
    if not dest.parent.is_dir():
        raise FileNotFoundError(dest.parent)
    fit, _, hashes = load_external(root)
    fit.fit_background_in_file(
        str(src), str(dest), segment_size=int(segment_size),
        buffer=int(buffer), pruning=int(pruning),
        ana_fmin=float(ana_fmin), ana_fmax=float(ana_fmax),
        n_processes=int(n_processes), **vendor_options)
    return {"status": STATUS, "vendor_background_json": str(dest),
            "vendor_source_sha256": hashes,
            "input_sha256": hashlib.sha256(src.read_bytes()).hexdigest()}


def vendor_mle_from_prepared_ape(*, prepared, external_root,
                                 min_log10mu=-40, max_log10mu=-32,
                                 n_seed_lkl=1000):
    """Run vendor optimizer verbatim on prepared vendor arrays (one or many bins).

    Preparation must use the vendor's segment selection, quality cuts and
    peak-shape conventions, replacing ONLY vendor get_beta with calibrated
    APE optical response. This function does not pretend that preparation
    has been verified.
    """
    root = Path(external_root).expanduser().resolve(strict=True)
    load_external(root)
    from scalar_dark_matter_lpsd.lvk.max_lkl_executable import calc_max_lkl
    Y = np.asarray(prepared["Y"], dtype=float)
    bkg = np.asarray(prepared["bkg"], dtype=float)
    pars = np.asarray(prepared["model_args"], dtype=float)
    norm = np.asarray(prepared["peak_norm"], dtype=float)
    shape = np.asarray(prepared["peak_shape"], dtype=float)
    if Y.ndim != 2 or Y.size == 0 or bkg.shape != Y.shape or norm.shape != Y.shape:
        raise ValueError("Require vendor arrays of shape (records, bins)")
    if pars.shape != (Y.shape[0], 3) or shape.shape != (Y.shape[1],):
        raise ValueError("Vendor skewnormal parameters or peak shape mismatch")
    if not all(np.isfinite(x).all() for x in (Y, bkg, pars, norm, shape)):
        raise ValueError("Nonfinite input")
    if np.any(norm < 0) or np.any(shape < 0) or not np.isclose(shape.sum(), 1.):
        raise ValueError("Invalid normalized vendor peak shape or APE response")
    result = calc_max_lkl(Y, bkg, pars, norm, shape,
                          min_log10mu=min_log10mu,
                          max_log10mu=max_log10mu,
                          n_seed_lkl=n_seed_lkl)
    return {"status": STATUS, "zero_lkl": float(result[0]),
            "max_lkl": float(result[1]), "mu_hat": float(result[2]),
            "sigma_mu": float(result[3])}
