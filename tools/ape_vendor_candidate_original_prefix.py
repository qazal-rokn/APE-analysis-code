"""APE candidate selection, copied verbatim from original get_candidates.main.

Only excluded: LIGO DataManager, H1/L1 consistency tests, and plots.
Original vendor GPLv3-or-later, copyright 2026 Alexandre Goettel.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from scipy.stats import norm
from scalar_dark_matter_lpsd.LPSDIO import MaxLKLIO
from scalar_dark_matter_lpsd.lvk.get_candidates import cluster


def run_vendor_candidate_prefix(max_lkl_path, output_csv, sigma_limit=5.):
    # Original get_candidates.main: exact statistical statements.
    lkl = MaxLKLIO(infile=str(max_lkl_path))
    q0 = np.zeros_like(lkl.zero_lkl)
    mask = lkl.mu > 0
    q0[mask] = -2 * (lkl.zero_lkl[mask] - lkl.max_lkl[mask])

    # Cleaning
    mask = np.isnan(q0) | (q0 >= 500)
    q0 = q0[~mask]

    # Find 5sigma candidates (correct for look-elsewhere-effect)
    significance_sigma = -norm.ppf(norm.cdf(-sigma_limit) / len(q0))
    print(f"Corrected significance: {significance_sigma:.2f} sigma")
    candidate_mask = q0 > significance_sigma ** 2
    print(f"Number of candidates: {sum(candidate_mask)}")

    # Clustering
    candidate_idx, candidate_lambda = cluster(candidate_mask, lkl.mu[~mask])
    print(f"Number of clusters: {len(candidate_idx)}")

    # Original candidate frequency and amplitude assignments, relocated
    # immediately after clustering because intervening code is H1/L1-only.
    candidate_freqs = lkl.frequency[~mask][candidate_idx]
    candidate_lambdas = candidate_lambda

    # Original vendor CSV formatting; original writes after LIGO vetoes.
    dest = Path(output_csv)
    if dest.exists():
        raise FileExistsError(dest)
    candidate_frequencies = []
    for freq in candidate_freqs:
        candidate_frequencies.append(freq)
    with open(dest, "w") as f:
        f.write("\n".join([f"{x:.6e}" for x in candidate_frequencies]))
    return {"candidate_csv": str(dest), "n_clusters": int(len(candidate_idx)),
            "stage": "original vendor pre-H1/L1-veto candidate selection"}
