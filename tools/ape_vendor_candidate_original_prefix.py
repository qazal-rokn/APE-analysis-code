"""APE single-channel candidate selection: call original vendor cluster directly.

Only the detector-independent statements from original get_candidates.main
are required. Its complete main cannot be called: it unconditionally creates
a two-detector DataManager and runs H1/L1 consistency and plots.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import norm
from scalar_dark_matter_lpsd.LPSDIO import MaxLKLIO
from scalar_dark_matter_lpsd.lvk.get_candidates import cluster


def run_vendor_candidate_prefix(max_lkl_path, output_csv, sigma_limit=5.):
    lkl = MaxLKLIO(infile=str(max_lkl_path))
    q0 = np.zeros_like(lkl.zero_lkl)
    mask = lkl.mu > 0
    q0[mask] = -2 * (lkl.zero_lkl[mask] - lkl.max_lkl[mask])
    mask = np.isnan(q0) | (q0 >= 500)
    q0 = q0[~mask]
    if not len(q0):
        raise ValueError("No valid vendor q0 bins")
    significance_sigma = -norm.ppf(norm.cdf(-sigma_limit) / len(q0))
    candidate_mask = q0 >= significance_sigma**2
    candidate_idx, candidate_lambda = cluster(candidate_mask, lkl.mu[~mask])
    frequencies = lkl.frequency[~mask][candidate_idx]
    np.savetxt(output_csv, frequencies, fmt="%.12e")
    return {"candidate_csv": str(output_csv), "n_clusters": int(len(frequencies)),
            "original_vendor_function": "get_candidates.cluster"}
