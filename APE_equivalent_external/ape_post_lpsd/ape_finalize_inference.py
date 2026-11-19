"""Minimal APE-only candidate screening and median-expected upper-limit bridge.

Original vendor code is imported and called unchanged:
- LPSDIO.MaxLKLIO
- lvk.get_candidates.cluster (with a terminal False sentinel to close last cluster)
- lvk.calc_upper_limit.get_upper_limits_approx

The original LIGO main() scripts are NOT invoked: they require H1/L1 and
LIGO strain transfer functions. No new likelihood, statistics or optical model.

Results are SOFTWARE-DIAGNOSTIC ONLY until the APE optical response, native-C
kernel, halo and statistical coverage have been independently validated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

from scalar_dark_matter_lpsd.LPSDIO import MaxLKLIO
from scalar_dark_matter_lpsd.lvk.get_candidates import cluster
from scalar_dark_matter_lpsd.lvk.calc_upper_limit import get_upper_limits_approx


class APEInferenceError(ValueError):
    """Unsafe inference interface input, without changing vendor algorithms."""


def _hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def _read_validated_inputs(maxlkl_path: Path, provenance_path: Path):
    if not maxlkl_path.is_file() or not provenance_path.is_file():
        raise APEInferenceError('Step-4 MaxLKLIO output and sidecar JSON required')
    with provenance_path.open(encoding='utf-8') as stream:
        provenance = json.load(stream)
    if not isinstance(provenance, dict) or provenance.get('channel') not in ('X', 'Y', 'Q'):
        raise APEInferenceError('Step-4 sidecar must declare APE channel X or Y')
    if provenance.get('joint_xy_likelihood_used') is not False:
        raise APEInferenceError('Joint X/Y likelihood not approved')
    if provenance.get('signal_model') != 'single-bin additive PSD: exp(bkg)+mu*peak_norm; peak_shape=1':
        raise APEInferenceError('Unexpected Step-4 signal model')
    metadata = provenance.get('response_metadata')
    if not isinstance(metadata, dict) or metadata.get('response_units') != 'rad^2/Hz per unit mu':
        raise APEInferenceError('Missing APE PSD response units')
    mu_definition = metadata.get('mu_definition', '')
    g_ref = None
    if mu_definition == 'mu=(g_a_gamma_gamma/reference_coupling_gev_inv)^2':
        g_ref = metadata.get('reference_coupling_gev_inv')
        if not isinstance(g_ref, (float, int)) or not np.isfinite(g_ref) or g_ref <= 0:
            raise APEInferenceError('APE coupling parameter requires finite positive g_ref')
        g_ref = float(g_ref)
    data = MaxLKLIO(infile=str(maxlkl_path)).data
    if (data.ndim != 2 or data.shape[1] != 5 or len(data) == 0
            or not np.isfinite(data[:, 0]).all() or np.any(data[:, 0] <= 0)
            or np.any(np.diff(data[:, 0]) <= 0)):
        raise APEInferenceError('Original MaxLKLIO rows must be nonempty, ordered, unique frequencies')
    return data, provenance, g_ref


def original_candidate_screen(data: np.ndarray, sigma_limit: float = 5.0):
    """Reuse vendor q0/Bonferroni screening and unchanged cluster function.

    No H1/L1 consistency veto is possible for APE quadratures without a
    validated joint-signal/noise model. This is *screening*, not candidates
    cleared by the original full LIGO pipeline or calibrated global p-values.
    """
    if not np.isfinite(sigma_limit) or sigma_limit <= 0:
        raise APEInferenceError('Positive finite screening sigma threshold required')
    frequency, zero, best, mu, _ = data.T
    q0 = np.zeros(len(data), dtype=float)
    positive = mu > 0
    q0[positive] = -2.0 * (zero[positive] - best[positive])
    # Same raw cleaning rule as original get_candidates.main().
    rejected = np.isnan(q0) | (q0 >= 500)
    frequency, q0, mu = frequency[~rejected], q0[~rejected], mu[~rejected]
    if len(q0) == 0:
        raise APEInferenceError('No usable bins for original candidate screening')
    correction = float(-norm.ppf(norm.cdf(-sigma_limit) / len(q0)))
    active = (q0 > correction**2)
    # The vendor cluster() does not append a cluster that ends on the final
    # element. Add an inactive sentinel rather than modifying vendor source.
    indices, sqrt_mu = cluster(np.r_[active, False], np.r_[mu, 0.0])
    result = np.column_stack([frequency[indices], q0[indices],
                              mu[indices], sqrt_mu]) if len(indices) else np.empty((0, 4))
    return result, {
        'screening_sigma_requested': float(sigma_limit),
        'vendor_bonferroni_sigma_threshold': correction,
        'bins_tested_after_original_q0_cleaning': len(q0),
        'bins_above_threshold': int(np.count_nonzero(active)),
        'clusters_retained': len(indices),
        'cross_quadrature_veto_applied': False,
        'look_elsewhere_calibration': 'ORIGINAL VENDOR HEURISTIC; APE EFFECTIVE TRIALS NOT VALIDATED',
    }


def original_expected_ul(data: np.ndarray, *, g_ref: float | None = None,
                         alpha_cl: float = .95, tol: float = 1e-10):
    """Delegate original expected q_mu asymptotic solver *without alteration*.

    The original solver returns sqrt(mu_UL) and its approximate uncertainty;
    it does NOT calculate a data-dependent observed upper limit. Physical
    g conversion is allowed only when supplied mu=(g/g_ref)^2 metadata.
    """
    if not (np.isfinite(alpha_cl) and 0 < alpha_cl < 1 and
            np.isfinite(tol) and tol > 0):
        raise APEInferenceError('Invalid expected limit confidence level/tolerance')
    out = np.full((len(data), 5), np.nan, dtype=np.float64)
    out[:, 0] = data[:, 0]
    # Columns: frequency, vendor sqrt(mu_UL), vendor sqrt(mu) uncertainty,
    #          scaled g_UL GeV^-1, scaled vendor uncertainty GeV^-1.
    for i, sigma in enumerate(data[:, 4]):
        if not np.isfinite(sigma) or sigma <= 0:
            continue  # Preserve invalid fits as NaN; do not silently smooth.
        ul, sigma_sqrt_mu = get_upper_limits_approx(float(sigma), alpha_CL=alpha_cl, tol=tol)
        ul, sigma_sqrt_mu = float(np.asarray(ul).item()), float(np.asarray(sigma_sqrt_mu).item())
        if not np.isfinite(ul) or ul <= 0:
            raise APEInferenceError('Original expected-limit function returned invalid result')
        out[i, 1:3] = [ul, sigma_sqrt_mu]
        if g_ref is not None:
            out[i, 3:5] = [g_ref * ul, g_ref * sigma_sqrt_mu]
    return out


def finalize_single_quadrature(maxlkl_path: str | Path, step4_manifest_json: str | Path,
                               output_prefix: str | Path, *, sigma_limit: float = 5.,
                               alpha_cl: float = .95, tol: float = 1e-10):
    """Produce APE-only provisional screening and vendor expected-limit files."""
    maxlkl = Path(maxlkl_path).expanduser().resolve()
    manifest = Path(step4_manifest_json).expanduser().resolve()
    prefix = Path(output_prefix).expanduser().resolve()
    targets = [Path(str(prefix) + suffix) for suffix in
               ('_screened_candidates.csv', '_expected_ul.csv', '_analysis_status.json')]
    if any(path.exists() for path in targets):
        raise APEInferenceError('Refusing to overwrite candidate/upper-limit output')
    if not prefix.parent.is_dir():
        raise APEInferenceError('Output parent directory must already exist')
    data, step4, g_ref = _read_validated_inputs(maxlkl, manifest)
    candidates, screen_details = original_candidate_screen(data, sigma_limit)
    limits = original_expected_ul(data, g_ref=g_ref, alpha_cl=alpha_cl, tol=tol)
    status = {
        'analysis': 'APE separate-quadrature, single-bin original-vendor-statistics SOFTWARE DIAGNOSTIC',
        'channel': step4['channel'],
        'maxlkl_file_sha256': _hash(maxlkl),
        'step4_manifest_sha256': _hash(manifest),
        'response_metadata': step4['response_metadata'],
        'reference_coupling_gev_inv': g_ref,
        'candidate_screening': screen_details,
        'upper_limit_type': 'median EXPECTED asymptotic upper limit, original get_upper_limits_approx',
        'upper_limit_confidence_level_nominal': float(alpha_cl),
        'observed_upper_limit_calculated': False,
        'physics_validated': False,
        'statistical_coverage_validated_for_APE': False,
        'LIGO_detector_consistency_checks_applied': False,
        'required_remaining_validation': 'QWP response and phase; native C D PSD kernel; halo; X/Y and bin covariance; coverage',
    }
    try:
        np.savetxt(targets[0], candidates, delimiter=',', comments='# ',
                   header='frequency_Hz,q0_vendor,signed_mu_hat,sqrt_mu_cluster', fmt='%.17g')
        np.savetxt(targets[1], limits, delimiter=',', comments='# ',
                   header='frequency_Hz,vendor_sqrt_mu_UL,vendor_sqrt_mu_uncertainty,g_reference_scaled_UL_GeV-1,g_reference_scaled_uncertainty_GeV-1', fmt='%.17g')
        with targets[2].open('x', encoding='utf-8') as stream:
            json.dump(status, stream, indent=2)
            stream.write('\n')
    except BaseException:
        for path in targets:
            path.unlink(missing_ok=True)
        raise
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maxlkl-npy', required=True)
    parser.add_argument('--step4-manifest-json', required=True)
    parser.add_argument('--output-prefix', required=True)
    parser.add_argument('--sigma-limit', type=float, default=5.)
    parser.add_argument('--alpha-cl', type=float, default=.95)
    parser.add_argument('--tol', type=float, default=1e-10)
    args = parser.parse_args()
    status = finalize_single_quadrature(
        args.maxlkl_npy, args.step4_manifest_json, args.output_prefix,
        sigma_limit=args.sigma_limit, alpha_cl=args.alpha_cl, tol=args.tol)
    print(json.dumps(status, indent=2))


if __name__ == '__main__':
    main()
