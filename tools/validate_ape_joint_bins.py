"""One bounded, reproducible APE joint-bin statistical benchmark.

This is NOT an axion confidence-limit pipeline. The simulated observation is
the *complex Gaussian Fourier coefficient* of windowed stationary white noise,
with an optional deterministic coherent mean field whose per-bin mean power
matches the supplied APE template. The shared field phase is one per draw
across all bins; it is NOT an independently randomized bin-power model.

Noise-only log powers have a log-exponential marginal; a skew-normal is fitted
once as a diagnostic approximation to the original vendor's marginal model.
Both the correlated-bin experiment and independent-bin control use that SAME
marginal, injection, settings, sample counts, and likelihood function.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2, skewnorm

from ape_c_lpsd.analysis.vendor_statistics import (
    factor_power_template, load_external,
)
from physics.qwp_signal import qwp_signal_power_template


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def kaiser_correlation(frequencies_hz, psll_db, *, integration_nodes=2049):
    """Same-length-window circular Gaussian-noise coefficient covariance.

    Centered analytical Kaiser window. For neighbouring actual C bins with
    slightly different segment lengths this is an explicit approximation.
    Frequencies must be from one strictly increasing native C candidate window.
    """
    f = np.asarray(frequencies_hz, dtype=float)
    if (f.ndim != 1 or not 2 <= len(f) <= 32
            or not np.isfinite(f).all() or np.any(f <= 0)
            or np.any(np.diff(f) <= 0)
            or not 25 <= psll_db <= 250):
        raise ValueError("Invalid bin grid or Kaiser PSLL")
    d = np.diff(f)
    if np.max(d) / np.min(d) > 1.1:
        raise ValueError("Window is not locally uniform enough for common-T benchmark")
    # Exact original genwin.c Kaiser alpha interpolation (not user normalization).
    x = psll_db / 100.0
    alpha = ((0.0889732*x - 0.493285)*x + 4.71469)*x - 0.0821377
    t = np.linspace(-0.5, 0.5, integration_nodes)
    w2 = np.i0(np.pi*alpha*np.sqrt(np.maximum(0, 1-(2*t)**2)))**2
    w2 /= np.sum(w2)
    common_t = 1.0 / np.median(d)
    offsets = f[:, None] - f[None, :]
    gamma = np.exp(2j*np.pi*offsets[..., None]*common_t*t).dot(w2)
    gamma = (gamma + gamma.conj().T)/2
    if not np.allclose(np.diag(gamma), 1, rtol=0, atol=2e-12):
        raise RuntimeError("Invalid Fourier coefficient unit variance")
    eig = np.linalg.eigvalsh(gamma)
    if eig[0] <= 0:
        raise ValueError("Non-positive Gaussian noise covariance; grid not supported")
    return gamma, common_t, alpha


def draw_noise(rng, count, covariance, *, xy_corr=0):
    """Proper complex Gaussian coefficients, shape (draw, X/Y, bin)."""
    if not 0 <= xy_corr < 1:
        raise ValueError("Use explicitly declared 0 <= X/Y coefficient correlation < 1")
    n = covariance.shape[0]
    chol = np.linalg.cholesky(covariance)
    white = (rng.standard_normal((count, 2, n))
             + 1j*rng.standard_normal((count, 2, n)))/np.sqrt(2)
    z = white @ chol.T
    z[:, 1, :] = xy_corr*z[:, 0, :] + np.sqrt(1-xy_corr**2)*z[:, 1, :]
    return z


def loglikelihood_grid(power, background, template, params, mu_grid, *, batch=64):
    """Vectorized *numerical parity* of unchanged vendor skew-normal product.

    Not a new inference likelihood: evaluates the same one-bin marginal
    skew-normal expression and the same sum over channels and spectral bins.
    """
    p = np.asarray(power, dtype=float)
    bg = np.asarray(background, dtype=float)
    shape = np.asarray(template, dtype=float)
    mus = np.asarray(mu_grid, dtype=float)
    if (p.ndim != 3 or p.shape[1:] != bg.shape or bg.shape != shape.shape
            or mus.ndim != 1 or not np.isfinite(mus).all()
            or np.any(mus < 0) or
            not np.isfinite(p).all() or np.any(p <= 0)
            or not np.isfinite(bg).all() or np.any(bg <= 0)
            or not np.isfinite(shape).all() or np.any(shape < 0)):
        raise ValueError("Invalid likelihood arrays")
    alpha, loc, scale = params
    out = np.empty((len(p), len(mus)))
    model = np.log(bg[None, None, :, :] + mus[None, :, None, None]*shape[None, None, :, :])
    for start in range(0, len(p), batch):
        stop = min(start + batch, len(p))
        residual = np.log(p[start:stop, None, :, :]) - model
        out[start:stop] = skewnorm.logpdf(
            residual, a=alpha, loc=loc, scale=scale).sum(axis=(2, 3))
    if not np.isfinite(out).all():
        raise ValueError("Nonfinite vendor-equivalent log-likelihood")
    return out


def make_summary(ll_null, ll_signal, mu_grid, injected_mu):
    """Empirical null q0 and nominal likelihood-ratio interval coverage."""
    threshold = float(chi2.ppf(.90, 1))  # Chernoff 50:50 mixture; nominal local 5%.
    ci_threshold = float(chi2.ppf(.95, 1))  # Wilks nominal two-sided 95%.
    def fit(ll):
        best = ll.argmax(axis=1)
        return mu_grid[best], 2*np.maximum(0, ll[np.arange(len(ll)), best]-ll[:, 0])
    null_mu, null_q0 = fit(ll_null)
    signal_mu, signal_q0 = fit(ll_signal)
    inj_index = np.flatnonzero(np.isclose(mu_grid, injected_mu, atol=1e-12, rtol=0))
    if len(inj_index) != 1:
        raise ValueError("Injected mu must occur exactly once on the prespecified grid")
    q_true = 2*(np.max(ll_signal, axis=1)-ll_signal[:, inj_index[0]])
    return {
        "null": {
            "n": int(len(ll_null)),
            "q0_95pct_empirical": float(np.quantile(null_q0, .95)),
            "false_positive_rate_nominal_5pct": float(np.mean(null_q0 > threshold)),
            "nominal_5pct_threshold": threshold,
            "null_best_mu_boundary_fraction": float(np.mean(null_mu == mu_grid[-1])),
        },
        "injection": {
            "n": int(len(ll_signal)),
            "mu_true": float(injected_mu),
            "median_mu_hat": float(np.median(signal_mu)),
            "mean_mu_hat": float(np.mean(signal_mu)),
            "bias_mean_mu": float(np.mean(signal_mu)-injected_mu),
            "nominal_95pct_profile_coverage": float(np.mean(q_true <= ci_threshold)),
            "nominal_95pct_threshold": ci_threshold,
            "mu_hat_upper_grid_boundary_fraction": float(np.mean(signal_mu == mu_grid[-1])),
        },
    }


def run(args):
    file = Path(args.template).expanduser().resolve(strict=True)
    if file.suffix != ".npz":
        raise ValueError("Existing multi-bin physical template .npz required")
    config_file = Path(args.config).expanduser().resolve(strict=True)
    config = json.loads(config_file.read_text())
    if not args.synthetic_only:
        raise ValueError("Require explicit --synthetic-only; never use this for real data")
    if (not 0 < args.snr_at_mu1 < 100 or args.mu_injected <= 0
            or args.n_null < 128 or args.n_signal < 128
            or args.n_grid < 20 or args.mu_max <= args.mu_injected):
        raise ValueError("Invalid predeclared ensemble or signal-to-background settings")

    with np.load(file, allow_pickle=False) as data:
        expected = {"frequency_support_hz", "normalized_halo_weights",
                    "observed_frequency_windows_hz", "window_power_response"}
        if set(data.files) != expected:
            raise ValueError("Template file lacks approved multi-bin data contract")
        support, weights, observed_f, response = (
            np.asarray(data[k], dtype=float) for k in
            ("frequency_support_hz", "normalized_halo_weights",
             "observed_frequency_windows_hz", "window_power_response"))
    if (observed_f.ndim != 2 or observed_f.shape[0] != 1
            or response.ndim != 4 or response.shape[0] != 1
            or response.shape[1] != 2 or response.shape[2] != observed_f.shape[1]):
        raise ValueError("One prespecified X/Y candidate window required")
    f = observed_f[0]
    if len(f) != args.window_bins:
        raise ValueError("Expected exact predeclared candidate-bin count")
    physics = config["physics"]
    template = qwp_signal_power_template(
        signal_frequency_hz=support, normalized_halo_weights=weights,
        window_power_response=response[0],
        quadrature_power_fraction=physics["quadrature_power_fraction"],
        reference_coupling_gev_inv=physics["reference_coupling_gev_inv"],
        finesse_s=physics["finesse_s"], tau_s=physics["tau_s"],
        tau1_s=physics["tau1_s"], rho_dm_gev4=physics["rho_dm_gev4"],
        unit_conversion=physics["unit_conversion"])
    if (not np.isfinite(template).all() or np.any(template < 0)
            or not np.any(template > 0)):
        raise ValueError("Invalid physical benchmark template")
    peak_norm, peak_shape = factor_power_template(template)
    np.testing.assert_allclose(peak_norm*peak_shape, template, rtol=1e-13, atol=0)
    psll = float(config["c_settings"]["psll_db"])
    gamma, common_t, alpha_window = kaiser_correlation(f, psll)
    independent = np.eye(len(f), dtype=complex)

    # Synthetic noise level is an explicit BENCHMARK SNR, not a signal calibration.
    baseline = float(np.median(template[template > 0]) / args.snr_at_mu1)
    B = np.full_like(template, baseline)
    grid = np.linspace(0, args.mu_max, args.n_grid)
    if not np.any(np.isclose(grid, args.mu_injected, atol=1e-12, rtol=0)):
        raise ValueError("Injection mu must be an exact prespecified grid point")

    rng = np.random.default_rng(args.seed)
    # Fit the null log(Exp(1)) marginal ONCE, shared by both test conditions.
    null_fit_sample = np.log(rng.exponential(size=30000))
    shape, loc, scale = map(float, skewnorm.fit(null_fit_sample))
    params = (shape, loc, scale)

    _, original_stats, vendor_hashes = load_external(args.vendor)
    results = {}
    for label, cov in (("correlated_kaiser", gamma), ("independent_bins_control", independent)):
        total = args.n_null + args.n_signal
        z = draw_noise(rng, total, cov, xy_corr=args.xy_corr)
        signal_phase = rng.uniform(0, 2*np.pi, size=(args.n_signal, 1, 1))
        amplitude = np.sqrt(args.mu_injected*template)[None, :, :]
        z[args.n_null:] += amplitude/np.sqrt(B)[None, :, :]*np.exp(1j*signal_phase)
        power = np.abs(z)**2 * B[None, :, :]

        # Original vendor parity check on actual (nonmocked) imported stats.py.
        model_args = np.tile(np.array([loc, scale, shape]), (2, 1))
        ll_check = loglikelihood_grid(power[:1], B, template, params, grid[:2])
        for k, mu in enumerate(grid[:2]):
            original = float(original_stats.log_likelihood(
                float(mu), np.log(power[0]), np.log(B),
                peak_norm, peak_shape, model_args))
            if not np.isclose(original, ll_check[0, k], rtol=0, atol=1e-9):
                raise RuntimeError("Original vendor log-likelihood numerical parity FAILED")

        all_ll = loglikelihood_grid(power, B, template, params, grid)
        summary = make_summary(all_ll[:args.n_null],
                               all_ll[args.n_null:], grid, args.mu_injected)
        measured = np.corrcoef(np.abs(z[:args.n_null, 0, :])**2, rowvar=False)
        summary["empirical_null_adjacent_power_correlation"] = float(np.mean(
            np.diag(measured, k=1)))
        results[label] = summary

    result = {
        "status": "joint_bin_statistical_benchmark_complete_not_ape_coverage_validation",
        "sources": {"template": str(file), "template_sha256": digest(file),
                    "config_sha256": digest(config_file), "vendor_hashes": vendor_hashes},
        "assumptions": {
            "noise": "stationary proper complex Gaussian white-noise Fourier coefficients",
            "bin_window": "same-length centered Kaiser approximation of actual close C bins",
            "xy_complex_noise_correlation": args.xy_corr,
            "signal": "fixed coherent complex template envelope, one common random phase per realization; NOT a stochastic SHM/APE signal",
            "null_logpower_marginal": "fitted skew-normal approximation to log exponential",
            "signal_to_background_median_at_mu1": args.snr_at_mu1,
            "background_power_rad2": baseline,
            "native_c_reexecuted": False,
            "additional_K_normalization": False,
            "calibrated_ape_95pct_coverage": False,
        },
        "prespecified": {
            "seed": args.seed, "n_null_per_condition": args.n_null,
            "n_signal_per_condition": args.n_signal,
            "mu_injected": args.mu_injected, "mu_grid_max": args.mu_max,
            "n_mu_grid": args.n_grid, "confidence_criterion": "nominal LR 3.841",
            "false_positive_criterion": "nominal local 5pct q0 2.706",
        },
        "geometry": {
            "frequency_hz": f.tolist(),
            "common_window_seconds_approximate": common_t,
            "kaiser_alpha": alpha_window,
            "theoretical_adjacent_power_correlation": float(abs(gamma[0, 1])**2),
            "background_units": "rad^2",
            "QWP_reference_signal_template_units": "rad^2",
        },
        "results": results,
        "interpretation": "Benchmark only: assess differences between correlated and independent bins; X/Y physics, physical axion halo, real-data coverage remain unvalidated.",
    }
    output = Path(args.outdir).expanduser().resolve()
    if output.exists():
        raise FileExistsError("Output directory must be brand-new: " + str(output))
    output.mkdir(parents=True, exist_ok=False)
    report = output/"joint_bin_benchmark_report.json"
    with report.open("x") as stream:
        json.dump(result, stream, indent=2)
    print("ORIGINAL_VENDOR_LIKELIHOOD_PARITY_PASS")
    print("CORRELATED_BIN_AND_INDEPENDENT_CONTROL_COMPLETE")
    for name, entry in results.items():
        print(name,
              "adjacent_power_corr", entry["empirical_null_adjacent_power_correlation"],
              "null_FPR_at_nominal_5pct", entry["null"]["false_positive_rate_nominal_5pct"],
              "mu_hat_mean", entry["injection"]["mean_mu_hat"],
              "nominal_profile_coverage", entry["injection"]["nominal_95pct_profile_coverage"])
    print("APE_JOINT_BIN_BENCHMARK_ONLY_NO_PHYSICAL_LIMIT")
    print("REPORT", report)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--template", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--vendor", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--synthetic-only", action="store_true")
    p.add_argument("--n-null", type=int, default=2048)
    p.add_argument("--n-signal", type=int, default=2048)
    p.add_argument("--seed", type=int, default=20261028)
    p.add_argument("--mu-injected", type=float, default=1.)
    p.add_argument("--mu-max", type=float, default=4.)
    p.add_argument("--n-grid", type=int, default=161)
    p.add_argument("--window-bins", type=int, default=7)
    p.add_argument("--xy-corr", type=float, default=0.)
    p.add_argument("--snr-at-mu1", type=float, default=0.3)
    run(p.parse_args())


if __name__ == "__main__":
    main()
