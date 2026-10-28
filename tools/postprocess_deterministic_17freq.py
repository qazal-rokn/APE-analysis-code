"""Deterministic 17-frequency APE synthetic native-C-PS postprocessor.

Reuses existing native-C output and original vendor fitting/likelihood/expected
limit. No C execution, simulated dark-matter field, global significance or
observed upper limits. The halo distribution and kernel are frozen to the
approved 2026-10-28 5-kHz synthetic benchmark, not a physical APE calibration.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.optimize import minimize_scalar
from scipy.special import roots_legendre

from ape_c_lpsd.analysis.vendor_statistics import (
    factor_power_template, fit_window, load_external,
)
from ape_c_lpsd.analysis.ape_candidates import select_candidates
from ape_c_lpsd.analysis.ape_expected_upper_limit import expected_upper_limit_ape
from physics.qwp_signal import qwp_amplitude_per_coupling
from tools.run_paper_axion_like_injection import window_power_kernel


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def execute(*, rho_x, rho_y, model_path, config_path, output_dir,
            reference_report=None):
    paths = [Path(p).expanduser().resolve(strict=True)
             for p in (rho_x, rho_y, model_path, config_path)]
    out = Path(output_dir).expanduser().resolve()
    if out.exists():
        raise FileExistsError(out)
    if not out.name.startswith("ape-postprocess-"):
        raise ValueError("Require new ape-postprocess-* directory")

    model = json.loads(paths[2].read_text())
    config = json.loads(paths[3].read_text())
    if (model.get("deterministic_DM_signal_power") is not True
            or model.get("random_DM_phase") is not False
            or model.get("random_DM_amplitude") is not False):
        raise ValueError("Deterministic signal power required; random DM forbidden")
    if config["external_vendor_root"] != "/source":
        raise ValueError("Frozen original vendor path required")
    if config.get("n_effective_trials") is not None:
        raise ValueError("Do not infer uncalibrated global significance")

    f_ref = np.asarray(model["native_frequency_bins_hz"], float)
    T_ref = np.asarray(
        model["deterministic_power_template_per_channel_rad2_at_gref"], float)
    if (f_ref.shape != (12,) or T_ref.shape != (2, 12)
            or not np.isfinite(T_ref).all() or np.any(T_ref < 0)):
        raise ValueError("Unexpected frozen 12-bin reference template")

    native = [np.loadtxt(p, ndmin=2) for p in paths[:2]]
    if any(a.ndim != 2 or a.shape[1] != 2 for a in native):
        raise ValueError("Native spectra require frequency/PS columns")
    np.testing.assert_array_equal(native[0][:, 0], native[1][:, 0])
    freq = native[0][:, 0]
    if len(freq) != 1301 or not np.all(np.diff(freq) > 0):
        raise ValueError("Unexpected native C grid")
    P = np.stack([a[:, 1] for a in native])
    if not np.isfinite(P).all() or np.any(P <= 0):
        raise ValueError("Native C PS must be positive and finite")

    vroot = config["external_vendor_root"]
    fits = [fit_window(freq, channel, vroot) for channel in P]
    if not all(x["vendor_quality_pass"] for x in fits):
        raise RuntimeError("Original vendor background quality gate failed")
    B = np.stack([x["background_psd"] for x in fits])
    pars = np.stack([x["skewnormal_loc_scale_alpha"] for x in fits])
    _, stats, hashes = load_external(vroot)

    idx_ref = np.searchsorted(freq, f_ref)
    if np.any(idx_ref >= len(freq)):
        raise ValueError("Reference frequency beyond C grid")
    np.testing.assert_array_equal(freq[idx_ref], f_ref)
    center_idx = int(np.argmin(abs(freq - 5000.0)))
    if abs(freq[center_idx] - 5000.0) > 1e-7:
        raise ValueError("Frozen reference center missing")
    offsets = idx_ref - center_idx
    if not np.array_equal(np.diff(idx_ref), np.ones(len(idx_ref)-1)):
        raise ValueError("Frozen model bins not adjacent")
    q_native = 1.0 / np.median(freq[1:] / freq[:-1] - 1.0)
    if not np.isclose(q_native, 1e6, rtol=1e-3):
        raise ValueError("Original native Q differs from frozen model")

    # Fixed deterministic truncated 3D Gaussian halo quadrature.
    theta = 229000.0
    vesc = 528000.0
    vlab = 229000.0
    sigma_v = theta * np.sqrt(1.5)
    vn, vw = roots_legendre(128)
    an, aw = roots_legendre(96)
    vg = 0.5 * vesc * (vn + 1)
    weight_v = (0.5 * vesc * vw * vg**2
                * np.exp(-vg**2 / (2 * sigma_v**2)))
    v2 = (vg[:, None]**2 + vlab**2
          - 2 * vg[:, None] * vlab * an[None, :])
    epsilon = (v2 / (2 * constants.c**2)).reshape(-1)
    masses = (weight_v[:, None] * 0.5 * aw[None, :]).reshape(-1)
    masses /= masses.sum()
    phys = config["physics"]
    fraction = np.asarray(phys["quadrature_power_fraction"], float)
    gref = float(phys["reference_coupling_gev_inv"])
    if (fraction.shape != (2,) or np.any(fraction < 0)
            or not np.isclose(fraction.sum(), 1.0)
            or gref <= 0):
        raise ValueError("Invalid frozen X/Y fractions or coupling")

    def template(fout, frest):
        support = frest * (1 + epsilon)
        amp = qwp_amplitude_per_coupling(
            support, finesse_s=phys["finesse_s"], tau_s=phys["tau_s"],
            tau1_s=phys["tau1_s"], rho_dm_gev4=phys["rho_dm_gev4"],
            unit_conversion=phys["unit_conversion"])
        W = window_power_kernel(fout, support, q=q_native, psll=238.13)
        return ((gref**2 / 2) * fraction[:, None]
                * (W @ (masses * amp**2))[None, :])

    np.testing.assert_allclose(template(f_ref, 5000.0), T_ref,
                               rtol=2e-10, atol=0)

    scan = []
    for ci in range(center_idx - 8, center_idx + 9):
        idx = ci + offsets
        if idx.min() < 0 or idx.max() >= len(freq):
            raise ValueError("Window beyond native C grid")
        rest = float(freq[ci])
        power = template(freq[idx], rest)
        peak_norm, peak_shape = factor_power_template(power)
        y, b = np.log(P[:, idx]), np.log(B[:, idx])

        def logl(mu):
            return float(stats.log_likelihood(
                float(mu), y, b, peak_norm, peak_shape, pars))

        positive = power > 0
        lower = -0.9 * float(np.min(B[:, idx][positive] / power[positive]))
        fit = minimize_scalar(lambda mu: -logl(mu),
                              bounds=(lower, 10.0), method="bounded",
                              options={"xatol": 1e-8})
        if not fit.success or fit.x <= lower + 1e-4 or fit.x >= 9.999:
            raise RuntimeError(f"Unreliable likelihood maximum at {rest}")
        muhat = float(fit.x)
        lmax = logl(muhat)
        q0 = max(0., 2 * (lmax - logl(0.))) if muhat > 0 else 0.
        sigma = float(stats.sigma_at_point(
            logl, muhat, initial_dx=0.02, tolerance=1e-4,
            max_iterations=40))
        if not np.isfinite(sigma) or sigma <= 0:
            raise RuntimeError(f"Vendor sigma invalid at {rest}")
        ul = expected_upper_limit_ape(
            sigma_mu=sigma, g_ref_gev_inv=gref, external_root=vroot,
            alpha_cl=0.95)
        scan.append(dict(
            frequency_hz=rest, mu_hat_signed=muhat,
            sigma_mu_vendor=sigma, local_sqrt_q0=float(np.sqrt(q0)),
            expected_g_ul_gev_inv=float(ul["g_ul_gev_inv"]),
            log_likelihood_null=float(logl(0.)),
            log_likelihood_max=float(lmax)))

    selected = select_candidates(
        [x["frequency_hz"] for x in scan],
        [x["log_likelihood_null"] for x in scan],
        [x["log_likelihood_max"] for x in scan],
        [x["mu_hat_signed"] for x in scan],
        local_sigma_threshold=3.0, n_effective_trials=None)

    result = {
        "status": "synthetic_native_C_deterministic_17_frequency_postprocessing",
        "frequency_count": 17, "frequencies": scan,
        "local_candidate_selection": selected,
        "expected_not_observed": True,
        "global_significance_validated": False,
        "correlated_noise_coverage_validated": False,
        "physical_ape_calibration_validated": False,
        "negative_MLE_observed_limits_resolved": False,
        "dark_matter_signal_stochasticity": False,
        "native_C_rerun": False, "original_vendor_modified": False,
        "provenance": {
            "rho_x_sha256": digest(paths[0]),
            "rho_y_sha256": digest(paths[1]),
            "model_sha256": digest(paths[2]),
            "config_sha256": digest(paths[3]),
            "vendor_hashes": hashes,
        },
    }
    if reference_report is not None:
        reference = json.loads(
            Path(reference_report).expanduser().resolve(strict=True).read_text())
        refrows = reference["expected_limits"]
        if len(refrows) != 17:
            raise RuntimeError("Frozen reference must have 17 hypotheses")
        errors = []
        for got, ref in zip(scan, refrows):
            if abs(got["frequency_hz"] - ref["frequency_hz"]) > 1e-7:
                raise RuntimeError("Frozen frequency alignment failure")
            errors.append(abs(got["expected_g_ul_gev_inv"] /
                              ref["expected_g_ul_gev_inv"] - 1))
        result["reference_max_relative_error"] = float(max(errors))
        if max(errors) > 0.02:
            raise RuntimeError("Frozen 17-frequency expected UL parity failed")

        # Only enforce the noise-only cluster criterion if the input
        # PS is exactly the frozen noise-only X/Y data by SHA256.
        frozen_x = "482b3106bd895bad9467b8dd2c67f5677ee5d81803d12997bcd9e6835c40d108"
        frozen_y = "a8d802edbefbd2c45042cdabbaa2645f54534904ceb5e41b7fa8eca312de58d7"
        if (result["provenance"]["rho_x_sha256"] == frozen_x
                and result["provenance"]["rho_y_sha256"] == frozen_y
                and selected["candidates"]):
            raise RuntimeError("Noise-only frozen benchmark generated a candidate")

    out.mkdir(parents=True, exist_ok=False)
    with (out / "summary.json").open("x") as file:
        json.dump(result, file, indent=2)
    with (out / "frequencies.csv").open("x", newline="") as file:
        columns = ["frequency_hz", "mu_hat_signed", "sigma_mu_vendor",
                   "local_sqrt_q0", "expected_g_ul_gev_inv"]
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows({c: row[c] for c in columns} for row in scan)
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
    report = execute(
        rho_x=args.rho_x, rho_y=args.rho_y, model_path=args.model,
        config_path=args.config, output_dir=args.output_dir,
        reference_report=args.reference_report)
    print("FREQUENCIES", report["frequency_count"])
    print("LOCAL_CANDIDATES", len(report["local_candidate_selection"]["candidates"]))
    print("EXPECTED_G_UL_5000", report["frequencies"][8]["expected_g_ul_gev_inv"])
    print("MAX_RELATIVE_REFERENCE_ERROR", report.get("reference_max_relative_error"))
    print("NATIVE_C_RERUN", report["native_C_rerun"])
    print("RANDOM_DM_SIGNAL", report["dark_matter_signal_stochasticity"])
    print("REPORT", str(Path(args.output_dir) / "summary.json"))


if __name__ == "__main__":
    main()
