"""Original-vendor single-bin MLE candidates, with APE-safe local clustering.

The original vendor optimizer and expected-UL calculation are called unchanged
through their previously tested adapters. Candidate q0 and contiguous grouping
are delegated to existing APE select_candidates (which fixes the vendor's
terminal-bin omission). No calibrated global significance or real-data use.
"""
from __future__ import annotations
import numpy as np
from ape_c_lpsd.analysis.ape_vendor_fit import fit_vendor_mle_and_expected_ul
from ape_c_lpsd.analysis.ape_candidates import select_candidates


def fit_and_report_single_bin_candidates(*, candidate_frequencies_hz,
                                         prepared_windows, external_root,
                                         min_log10mu, max_log10mu, n_seed_lkl,
                                         local_sigma_threshold, alpha_cl=0.95):
    """Fit independent selected windows, then rank local q0 candidates.

    Frequency is the candidate's positive, strictly increasing C-grid value.
    Each prepared window is the previously validated single-bin vendor input
    dictionary; no spectral construction or background fitting occurs here.
    The returned upper limit is *asymptotic expected*, NOT observed.
    """
    f = np.asarray(candidate_frequencies_hz, dtype=float)
    windows = list(prepared_windows)
    if (f.ndim != 1 or len(f) == 0 or len(windows) != len(f)
            or not np.isfinite(f).all() or np.any(f <= 0)
            or np.any(np.diff(f) <= 0)):
        raise ValueError('Require sorted positive frequencies and aligned windows')
    threshold = float(local_sigma_threshold)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError('Require positive finite local sigma threshold')

    fitted = [fit_vendor_mle_and_expected_ul(
        prepared=window, external_root=external_root,
        min_log10mu=min_log10mu, max_log10mu=max_log10mu,
        n_seed_lkl=n_seed_lkl, alpha_cl=alpha_cl)
        for window in windows]

    null = [r['log_likelihood_null'] for r in fitted]
    best = [r['log_likelihood_max'] for r in fitted]
    muhat = [r['mu_hat'] for r in fitted]
    selected = select_candidates(f, null, best, muhat,
                                 local_sigma_threshold=threshold,
                                 n_effective_trials=None)
    # No trial factor or cross-quadrature veto is assumed.
    return {
        'status': 'provisional_single_bin_local_candidates_no_global_significance_or_ape_coverage',
        'method': 'unchanged_vendor_MLE_q0_and_expected_UL_with_APE_terminal_safe_clustering',
        'frequencies_hz': f.tolist(),
        'fits': fitted,
        'selection': selected,
        'global_significance_calibrated': False,
        'ape_physical_coverage_validated': False,
    }
