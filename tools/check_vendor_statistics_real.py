"""Read-only integration check using real unchanged vendor fit and likelihood.

Uses the previously archived APE X 100-Hz 1100-bin validation window.
No spectrum recomputation, no vendor edits, no output file writes.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ape_c_lpsd.analysis.vendor_statistics import fit_window, log_likelihood, load_external

HOME = Path.home()
EXTERNAL = HOME / "Desktop/APE_Research/Scalar-Dark-Matter-LPSD-1.0.0"
SAVED = HOME / "APE-LPSD-validation/results/step09"
WINDOW = SAVED / "x_100Hz_values.csv"
REPORT = SAVED / "production_background_report.json"


def main():
    if not WINDOW.is_file() or not REPORT.is_file():
        raise FileNotFoundError("Missing original read-only Step 9 validation outputs")
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    record = next(r for r in report["records"]
                  if r["channel"] == "X" and float(r["center_requested_Hz"]) == 100.0)
    if not record.get("fit_quality_pass", False):
        raise RuntimeError("Reference X/100-Hz fit did not pass its quality cut")
    existing = np.loadtxt(WINDOW, delimiter=",", skiprows=1, ndmin=2)
    if existing.shape != (1100, 5) or not np.isfinite(existing).all():
        raise RuntimeError("Unexpected original validation window")
    f, psd, navg, saved_background, saved_residual = existing.T
    if np.any(psd <= 0) or np.any(saved_background <= 0) or np.any(navg < 1):
        raise RuntimeError("Invalid saved PSD / fit values")
    if not np.all(np.diff(f) > 0):
        raise RuntimeError("Saved frequency grid is not increasing")
    sha_before = hashlib.sha256(WINDOW.read_bytes()).hexdigest()

    settings = dict(record.get("released_code_configuration", report["released_code_configuration"]))
    fit = fit_window(f, psd, EXTERNAL, settings=settings)
    _, vendor_stats, vendor_hashes = load_external(EXTERNAL)
    core = slice(50, -50)
    # Original saved values were written to twelve significant digits and
    # optimizers can vary across installed SciPy versions. Report differences.
    ratio = fit["background_psd"][core] / saved_background[core]
    relative_background_difference = float(np.max(np.abs(ratio - 1.0)))
    residual_diff = float(np.max(np.abs(fit["log_residual"][core] - saved_residual[core])))
    pars_diff = float(np.max(np.abs(
        fit["skewnormal_loc_scale_alpha"]
        - np.asarray(record["skewnormal_loc_scale_alpha"], dtype=float))))
    chi_diff = abs(fit["vendor_reduced_chi2"]
                   - float(record["fit_quality_chi2_reduced"]))

    # Check the original likelihood implementation directly and via adapter.
    # This beta=B is the established algorithmic PSD-excess test, NOT
    # the APE axion signal model or a significance estimate.
    idx = 550
    observed = np.array([psd[idx]], dtype=float)
    background = np.array([fit["background_psd"][idx]], dtype=float)
    template = background.copy()
    parameters = fit["skewnormal_loc_scale_alpha"]
    mu = 0.25
    via_adapter = log_likelihood(
        mu, np.log(observed)[None, :], np.log(background)[None, :],
        1.0, np.ones(1), parameters[None, :], EXTERNAL)
    direct = float(vendor_stats.log_likelihood(
        mu, np.log(observed)[None, :], np.log(background)[None, :],
        1.0, template, parameters[None, :]))
    likelihood_difference = abs(via_adapter - direct)
    sha_after = hashlib.sha256(WINDOW.read_bytes()).hexdigest()
    result = {
        "window": "X quadrature, 100 Hz, 1100 saved bins (1000 core)",
        "external_source_sha256": vendor_hashes,
        "reference_quality_pass": True,
        "new_fit_quality_pass": fit["vendor_quality_pass"],
        "vendor_chi2": fit["vendor_reduced_chi2"],
        "reference_chi2": float(record["fit_quality_chi2_reduced"]),
        "absolute_chi2_difference": chi_diff,
        "max_core_relative_background_difference": relative_background_difference,
        "max_core_log_residual_difference": residual_diff,
        "max_skewnormal_parameter_difference": pars_diff,
        "likelihood_adapter_minus_direct_abs": likelihood_difference,
        "saved_window_unchanged": sha_before == sha_after,
        "interpretation": "Algorithmic integration only; not an APE axion significance test",
    }
    print(json.dumps(result, indent=2, allow_nan=False))
    if sha_before != sha_after:
        raise AssertionError("Read-only validation window unexpectedly changed")
    if not np.isfinite([chi_diff, residual_diff, pars_diff,
                        relative_background_difference,
                        likelihood_difference]).all():
        raise AssertionError("Nonfinite statistical integration result")
    if likelihood_difference > 1e-10:
        raise AssertionError("Adapter likelihood differs from original function")
    print("PASS: unchanged external functions used; likelihood agrees directly.")
    print("Fit-versus-archive differences are reported, not silently declared exact.")


if __name__ == "__main__":
    main()
