"""Postprocess existing native C-LPSD X/Y spectra with the unchanged vendor.

Synthetic benchmark only. Does not execute C, generate stochastic dark-matter
signals, modify vendor sources, estimate global significance, or infer observed
limits. Requires an explicit frozen deterministic signal-template JSON.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from ape_c_lpsd.analysis.vendor_statistics import (
    factor_power_template, fit_window, load_external,
)
from ape_c_lpsd.analysis.ape_candidates import select_candidates
from ape_c_lpsd.analysis.ape_expected_upper_limit import expected_upper_limit_ape


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def execute(*, rho_x, rho_y, model_path, config_path, output_dir,
            reference_report=None):
    """Run one frozen 12-bin/17-hypothesis *synthetic* postprocessing analysis.

    The deterministic template is already precomputed for the central rest
    frequency. Neighbor hypotheses must supply corresponding templates explicitly
    in a future, scientifically approved template input; never shift or
    manufacture halo templates silently.
    """
    paths = [Path(x).expanduser().resolve(strict=True)
             for x in (rho_x, rho_y, model_path, config_path)]
    dest = Path(output_dir).expanduser().resolve()
    if dest.exists():
        raise FileExistsError(dest)
    if not dest.name.startswith("ape-postprocess-"):
        raise ValueError("Output directory must start with ape-postprocess-")
    model = json.loads(paths[2].read_text())
    cfg = json.loads(paths[3].read_text())

    if (model.get("deterministic_DM_signal_power") is not True
            or model.get("random_DM_phase") is not False
            or model.get("random_DM_amplitude") is not False):
        raise ValueError("Require explicitly deterministic dark-matter power model")
    if cfg.get("external_vendor_root") != "/source":
        raise ValueError("Require approved original vendor path /source")
    if cfg.get("n_effective_trials") is not None:
        raise ValueError("Uncalibrated global trial-factor correction prohibited")

    freqs = np.asarray(model["native_frequency_bins_hz"], dtype=float)
    signal = np.asarray(
        model["deterministic_power_template_per_channel_rad2_at_gref"],
        dtype=float,
    )
    if (freqs.ndim != 1 or len(freqs) != 12 or signal.shape != (2, len(freqs))
            or not np.isfinite(signal).all() or np.any(signal < 0)
            or not np.isfinite(freqs).all()
            or np.any(np.diff(freqs) <= 0)):
        raise ValueError("Invalid frozen deterministic 12-bin model")

    arrays = [np.loadtxt(p, ndmin=2) for p in paths[:2]]
    if any(x.ndim != 2 or x.shape[1] != 2 for x in arrays):
        raise ValueError("Native C spectrum must contain two columns")
    np.testing.assert_array_equal(arrays[0][:, 0], arrays[1][:, 0])
    grid = arrays[0][:, 0]
    idx = np.searchsorted(grid, freqs)
    if np.any(idx >= len(grid)):
        raise ValueError("Model outside native C grid")
    np.testing.assert_array_equal(grid[idx], freqs)
    if not np.all(np.diff(idx) == 1):
        raise ValueError("Frozen template bins must be consecutive")

    observed = np.stack([a[:, 1] for a in arrays])
    if not np.isfinite(observed).all() or np.any(observed <= 0):
        raise ValueError("Native C power must be finite and positive")
    fits = [fit_window(grid, p, cfg["external_vendor_root"])
            for p in observed]
    if not all(f["vendor_quality_pass"] for f in fits):
        raise RuntimeError("Original vendor background quality check failed")
    background = np.stack([f["background_psd"][idx] for f in fits])
    args = np.stack([f["skewnormal_loc_scale_alpha"] for f in fits])
    peak_norm, peak_shape = factor_power_template(signal)
    _, stats, vendor_hashes = load_external(cfg["external_vendor_root"])
    Y = np.log(observed[:, idx])
    B = np.log(background)

    def logl(mu):
        return float(stats.log_likelihood(
            float(mu), Y, B, peak_norm, peak_shape, args))

    positive = signal > 0
    lower = -0.90 * float(np.min(background[positive] / signal[positive]))
    optimum = minimize_scalar(lambda mu: -logl(mu),
                              method="bounded", bounds=(lower, 10.0),
                              options={"xatol": 1e-8})
    if not optimum.success or optimum.x <= lower + 1e-4 or optimum.x >= 9.999:
        raise RuntimeError("Signed likelihood optimum not reliable")
    muhat = float(optimum.x)
    lmax = logl(muhat)
    q0 = float(max(0., 2. * (lmax - logl(0)))) if muhat > 0 else 0.
    sigma = float(stats.sigma_at_point(
        logl, muhat, initial_dx=0.02, tolerance=1e-4, max_iterations=40))
    if not np.isfinite(sigma) or sigma <= 0:
        raise RuntimeError("Original vendor sigma invalid")

    gref = float(cfg["physics"]["reference_coupling_gev_inv"])
    expected = expected_upper_limit_ape(
        sigma_mu=sigma, g_ref_gev_inv=gref,
        external_root=cfg["external_vendor_root"], alpha_cl=0.95)
    # The source model is central-only. We do not claim 17-frequency coverage
    # or candidate-clustering results from this single-hypothesis execution.
    center = 5000.0
    if not any(np.isclose(freqs, center, rtol=0, atol=1e-7)):
        raise ValueError("Frozen central rest frequency not represented")

    result = {
        "status": "synthetic_c_ps_single_hypothesis_postprocessing",
        "rest_frequency_hz": center,
        "mu_hat_signed": muhat,
        "sigma_mu_original_vendor": sigma,
        "local_q0": q0,
        "local_sqrt_q0": float(np.sqrt(q0)),
        "expected_g_ul_gev_inv": expected["g_ul_gev_inv"],
        "confidence_level": 0.95,
        "expected_not_observed": True,
        "candidate_scan_or_global_significance": False,
        "correlated_native_c_coverage_validated": False,
        "physical_ape_response_validated": False,
        "no_random_dark_matter_signal": True,
        "native_C_rerun": False,
        "original_vendor_modified": False,
        "provenance": {
            "rho_x_sha256": _hash(paths[0]),
            "rho_y_sha256": _hash(paths[1]),
            "deterministic_model_sha256": _hash(paths[2]),
            "configuration_sha256": _hash(paths[3]),
            "vendor_source_hashes": vendor_hashes,
            "vendor_ul_sha256": expected["vendor_upper_limit_sha256"],
        },
    }

    if reference_report is not None:
        reference = json.loads(
            Path(reference_report).expanduser().resolve(strict=True).read_text())
        target = next(x for x in reference["expected_limits"]
                      if abs(x["frequency_hz"] - center) < 1e-7)
        ref_g = float(target["expected_g_ul_gev_inv"])
        difference = abs(result["expected_g_ul_gev_inv"] / ref_g - 1.)
        result["frozen_reference_relative_difference"] = difference
        # The vendor fit/sigma step is intentionally not replaced with a
        # hardcoded number. Require close parity or stop the integration.
        if difference > 0.02:
            raise RuntimeError(
                f"Frozen 5 kHz expected-limit mismatch: {difference:.3%}")

    dest.mkdir(parents=True, exist_ok=False)
    with (dest / "summary.json").open("x") as handle:
        json.dump(result, handle, indent=2)
    with (dest / "summary.csv").open("x", newline="") as handle:
        fields = ["rest_frequency_hz", "mu_hat_signed",
                  "sigma_mu_original_vendor", "local_sqrt_q0",
                  "expected_g_ul_gev_inv"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow({key: result[key] for key in fields})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho-x", required=True)
    parser.add_argument("--rho-y", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--reference-report")
    args = parser.parse_args()
    result = execute(
        rho_x=args.rho_x, rho_y=args.rho_y, model_path=args.model,
        config_path=args.config, output_dir=args.output_dir,
        reference_report=args.reference_report)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
