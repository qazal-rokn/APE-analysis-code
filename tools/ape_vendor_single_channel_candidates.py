"""Single-channel APE candidate selection, reusing vendor trial correction and cluster.

Only the H1/L1 consistency test and LIGO-specific candidate plots are omitted.
No custom significance estimator or clustering algorithm is introduced.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm

def select_candidates(max_lkl_file, output_csv, sigma_limit=5.0):
    from scalar_dark_matter_lpsd.LPSDIO import MaxLKLIO
    from scalar_dark_matter_lpsd.lvk.get_candidates import cluster
    lkl = MaxLKLIO(infile=str(max_lkl_file))
    q0 = np.zeros_like(lkl.zero_lkl)
    positive = lkl.mu > 0
    q0[positive] = -2 * (lkl.zero_lkl[positive] - lkl.max_lkl[positive])
    # Identical vendor cleaning and look-elsewhere correction.
    bad = np.isnan(q0) | (q0 >= 500)
    q0 = q0[~bad]
    if len(q0) == 0:
        raise ValueError("No valid likelihood bins remain after vendor cleaning")
    significance_sigma = -norm.ppf(norm.cdf(-sigma_limit) / len(q0))
    candidate_mask = q0 > significance_sigma ** 2
    # Preserve vendor cluster() as written, including its boundary behavior.
    indices, strengths = cluster(candidate_mask, lkl.mu[~bad])
    frequency = lkl.frequency[~bad][indices]
    target = Path(output_csv).expanduser().resolve()
    if target.exists():
        raise FileExistsError(target)
    if not target.parent.is_dir():
        raise NotADirectoryError(target.parent)
    np.savetxt(target, frequency, fmt="%.12e")
    return {
        "status": "EXPERIMENTAL_UNTESTED_NOT_VALIDATED",
        "method": "original_vendor_q0_cleaning_trial_factor_and_cluster",
        "input_bins": int(len(q0)),
        "corrected_significance_sigma": float(significance_sigma),
        "candidate_bins": int(np.sum(candidate_mask)),
        "candidate_clusters": int(len(indices)),
        "two_interferometer_consistency": "excluded_not_applicable_to_single_APE_channel",
        "candidate_csv": str(target),
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-lkl-file", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--sigma-limit", type=float, default=5)
    args = parser.parse_args()
    print(json.dumps(select_candidates(args.max_lkl_file, args.output_csv, args.sigma_limit), indent=2))
