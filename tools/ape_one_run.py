"""Unified entry point for the existing, *scientifically scoped* APE workflows.

This is a thin dispatcher, not a new implementation of spectral estimation,
background fitting, likelihood or limits. The underlying verified modules and
their scientific gates are unchanged.

Modes:
* synthetic-full: original integrated synthetic C PS -> vendor -> candidates.
* frozen-17freq: frozen deterministic 17-frequency native-C PS postprocessor.
* real-calibration: automatic real DXD -> calibrated HDF5 (Mac ARM64).
* real-calibration-reuse: verify completed 0.0284 HDF5 without reprocessing.
There is intentionally no purported production real-data limits mode yet.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import h5py

EXPECTED = {
    "cavity_transmission_w_per_v": 0.0284,
    "responsivity_a_per_w": 0.7,
    "transimpedance_gain_v_per_a": 1e6,
    "second_harmonic_pd_response": 0.7,
}
DATASETS = {"i0_w": "W", "eta0": "dimensionless",
            "rho_x": "rad", "rho_y": "rad"}


def _sha256(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def verify_completed_real_calibration(stage_result, *, verify_full_hash=False):
    """Validate immutable existing result; no scans over 1.2B samples."""
    source = Path(stage_result).expanduser().resolve(strict=True)
    report = json.loads(source.read_text())
    if report.get("status") != "calibration_complete_C_backend_pending_scientific_gate":
        raise ValueError("Not a successfully completed real calibration stage")
    if report.get("native_c_executed") is not False:
        raise ValueError("Unexpected C execution claim")
    if report.get("physical_g_limits_authorized") is not False:
        raise ValueError("Cannot treat this product as an approved physical limit")
    target = Path(report["calibrated_hdf5"]).resolve(strict=True)
    if target.parent != source.parent:
        raise ValueError("Stage-result HDF5 must reside in same immutable run directory")
    expected_sha = report["calibrated_hdf5_sha256"]
    if not isinstance(expected_sha, str) or len(expected_sha) != 64:
        raise ValueError("Invalid provenance hash")
    with h5py.File(target, "r") as file:
        n = int(file.attrs["sample_count"])
        fs = float(file.attrs["sample_rate_hz"])
        if n != int(report["sample_count"]) or fs != float(report["sample_rate_hz"]):
            raise ValueError("HDF5 and stage record disagree on rate/count")
        if str(file.attrs["manifest_sha256"]) != report["frozen_manifest_sha256"]:
            raise ValueError("Manifest provenance mismatch")
        if bool(file.attrs.get("x_y_combined", True)):
            raise ValueError("Combined X/Y not supported")
        calibration = json.loads(file.attrs["calibration_parameters_json"])
        for k, v in EXPECTED.items():
            if (k not in calibration or
                    not math.isclose(float(calibration[k]), v, rel_tol=1e-12)):
                raise ValueError(f"Calibration mismatch: {k}")
            if not math.isclose(float(report["calibration"][k]), v, rel_tol=1e-12):
                raise ValueError(f"Report calibration mismatch: {k}")
        for key, unit in DATASETS.items():
            dset = file[key]
            if dset.shape != (n,) or dset.dtype.name != "float64":
                raise ValueError(f"Bad calibrated dataset: {key}")
            if dset.attrs.get("unit") != unit:
                raise ValueError(f"Bad unit for {key}")
    if verify_full_hash and _sha256(target) != expected_sha:
        raise ValueError("Full HDF5 SHA-256 mismatch")
    return {
        "status": "verified_existing_real_calibration_only",
        "calibrated_hdf5": str(target),
        "samples_per_channel": n,
        "sample_rate_hz": fs,
        "hdf5_sha256_recorded": expected_sha,
        "full_content_sha256_recomputed": bool(verify_full_hash),
        "native_c_executed": False,
        "physical_limits_authorized": False,
        "next_unconnected_stage": "real_data_C_PS_and_inference_adapter",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)

    s = modes.add_parser("synthetic-full", help="Existing synthetic C + vendor controller")
    s.add_argument("--config", type=Path, required=True)

    f = modes.add_parser("frozen-17freq", help="Frozen deterministic 17-frequency analysis")
    f.add_argument("--rho-x", type=Path, required=True)
    f.add_argument("--rho-y", type=Path, required=True)
    f.add_argument("--model", type=Path, required=True)
    f.add_argument("--config", type=Path, required=True)
    f.add_argument("--output-dir", type=Path, required=True)
    f.add_argument("--reference-report", type=Path)

    r = modes.add_parser("real-calibration", help="Automatic DXD calibration on Mac")
    r.add_argument("--run-config", type=Path, required=True)
    r.add_argument("--execute-calibration", action="store_true")

    v = modes.add_parser("real-calibration-reuse", help="Read-only completed HDF5 preflight")
    v.add_argument("--stage-result", type=Path, required=True)
    v.add_argument("--verify-full-hash", action="store_true")

    args = parser.parse_args()
    if args.mode == "synthetic-full":
        from tools.run_ape_pipeline import execute
        config = json.loads(args.config.expanduser().resolve(strict=True).read_text())
        result = execute(config)
        print(json.dumps({
            "status": result["status"],
            "mode": args.mode,
            "scientific_scope": "synthetic_only",
        }, indent=2))
    elif args.mode == "frozen-17freq":
        from tools.postprocess_deterministic_17freq import execute
        result = execute(
            rho_x=args.rho_x, rho_y=args.rho_y, model_path=args.model,
            config_path=args.config, output_dir=args.output_dir,
            reference_report=args.reference_report)
        print(json.dumps({
            "mode": args.mode,
            "frequency_count": result["frequency_count"],
            "local_candidates": len(result["local_candidate_selection"]["candidates"]),
            "native_C_rerun": result["native_C_rerun"],
            "physical_limits_authorized": False,
            "summary": str(args.output_dir / "summary.json"),
        }, indent=2))
    elif args.mode == "real-calibration":
        from tools.run_ape_automatic import execute
        execute(args.run_config.expanduser().resolve(strict=True),
                execute_calibration=args.execute_calibration)
    elif args.mode == "real-calibration-reuse":
        print(json.dumps(verify_completed_real_calibration(
            args.stage_result, verify_full_hash=args.verify_full_hash), indent=2))
    else:
        raise AssertionError("Unreachable")


if __name__ == "__main__":
    main()
