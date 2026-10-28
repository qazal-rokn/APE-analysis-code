"""Derive a new APE HDF5 calibration from an immutable prior HDF5.

Exact algebraic re-expression when ONLY the transmission W/V calibration
changes and q, G, H2 and all RMS conventions remain fixed. No raw DXD
replay, interpolation, filtering, X/Y mixing, native C run or physical
coupling-limit claim. Planning is read-only; writing requires --execute.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import uuid
from pathlib import Path

import h5py
import numpy as np

DATASETS = {
    "i0_w": "W",
    "eta0": "dimensionless",
    "rho_x": "rad",
    "rho_y": "rad",
}
SAFETY_MARGIN_BYTES = 512 * 1024**2


def _equal(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=0)


def _inspect(source, old_c, new_c, h2):
    with h5py.File(source, "r") as h5:
        if str(h5.attrs.get("product", "")) != "ape_frozen_calibrated_segment":
            raise ValueError("Unexpected source calibration product")
        if str(h5.attrs.get("calibration_method", "")) != "ape_rms_quadrature_calibration_v1":
            raise ValueError("Requires original RMS quadrature calibration v1")
        if bool(h5.attrs.get("x_y_combined", True)):
            raise ValueError("X/Y combination in source is forbidden")
        for attr in ("filtering_applied", "resampling_applied",
                     "interpolation_applied", "detrending_applied"):
            if bool(h5.attrs.get(attr, True)):
                raise ValueError(f"Source has {attr}=true")
        parameters = json.loads(h5.attrs["calibration_parameters_json"])
        checked = {
            "cavity_transmission_w_per_v": old_c,
            "responsivity_a_per_w": 0.7,
            "transimpedance_gain_v_per_a": 1e6,
            "second_harmonic_pd_response": h2,
        }
        for name, expected in checked.items():
            if name not in parameters or not _equal(parameters[name], expected):
                raise ValueError(f"Source calibration {name} does not match approved input")
        n = int(h5.attrs["sample_count"])
        fs = float(h5.attrs["sample_rate_hz"])
        if n <= 0 or not math.isfinite(fs) or fs <= 0:
            raise ValueError("Invalid sample count/rate")
        for name, unit in DATASETS.items():
            if name not in h5:
                raise ValueError(f"Missing {name}")
            ds = h5[name]
            if (not isinstance(ds, h5py.Dataset) or ds.shape != (n,)
                    or ds.dtype != np.dtype("float64")
                    or ds.attrs.get("unit") != unit):
                raise ValueError(f"Invalid dataset shape, type or unit: {name}")
        manifest = str(h5.attrs.get("manifest_sha256", ""))
        if len(manifest) != 64 or any(ch not in "0123456789abcdefABCDEF" for ch in manifest):
            raise ValueError("Missing valid frozen-manifest SHA256")
        if not h5.attrs.get("source_parts_json"):
            raise ValueError("Missing original source-part provenance")
        return {"samples": n, "sample_rate_hz": fs,
                "source_manifest_sha256": manifest,
                "source_calibration_parameters": parameters,
                "source_code_commit": str(h5.attrs.get("code_commit", "unknown"))}


def run(args):
    source = Path(args.source).expanduser().resolve(strict=True)
    target = Path(args.output).expanduser().resolve()
    if not source.is_file() or source == target:
        raise ValueError("Source/target must be separate files")
    if target.exists():
        raise FileExistsError(f"Target already exists: {target}")
    if not target.parent.is_dir():
        raise FileNotFoundError(f"Target directory must already exist: {target.parent}")
    provenance = target.with_suffix(target.suffix + ".provenance.json")
    if provenance.exists():
        raise FileExistsError(provenance)
    old_c, new_c, h2 = (float(args.expected_old_c_trans),
                        float(args.new_c_trans), float(args.h2))
    if (not all(math.isfinite(x) and x > 0 for x in (old_c, new_c, h2))
            or _equal(old_c, new_c)):
        raise ValueError("Require explicit different positive old/new C_trans and H2")
    meta = _inspect(source, old_c, new_c, h2)
    k = new_c / old_c
    factor = 1.0 / math.sqrt(k)
    required = 4 * meta["samples"] * np.dtype("float64").itemsize
    free = shutil.disk_usage(target.parent).free
    plan = {
        "status": "plan_only_unless_execute_is_set",
        "source": str(source), "target": str(target),
        "source_samples_per_dataset": meta["samples"],
        "sample_rate_hz": meta["sample_rate_hz"],
        "calibration": {
            "old_c_trans_w_per_v": old_c,
            "new_c_trans_w_per_v": new_c,
            "unchanged_q_a_per_w": 0.7,
            "unchanged_G_v_per_a": 1e6,
            "unchanged_H2": h2,
            "H2_measurement_provenance": "researcher working acceptance; independent measurement not supplied",
        },
        "exact_algebra": {
            "i0_w_new_over_old": k,
            "eta0_new_over_old": factor,
            "rho_x_new_over_old": factor,
            "rho_y_new_over_old": factor,
        },
        "minimum_output_bytes": required,
        "available_bytes_in_output_filesystem": free,
        "safety_margin_bytes": SAFETY_MARGIN_BYTES,
        "native_C_rerun": False,
        "raw_DXD_replay": False,
        "original_input_modified": False,
        "physical_g_limits_authorized": False,
    }
    print(json.dumps(plan, indent=2))
    if not args.execute:
        print("PLAN_COMPLETE_NO_FILES_WRITTEN")
        return
    if free < required + SAFETY_MARGIN_BYTES:
        raise RuntimeError(
            "Insufficient space for a new, uncompressed calibrated HDF5; "
            "source remains untouched")
    block = int(args.chunk_samples)
    if block <= 0 or block > 2_000_000:
        raise ValueError("Chunk size must be 1..2000000")

    temp = target.with_name(target.name + ".partial-" + uuid.uuid4().hex)
    hashes = {}
    try:
        with h5py.File(source, "r") as old, h5py.File(temp, "x") as new:
            for name, value in old.attrs.items():
                new.attrs[name] = value
            parameters = dict(meta["source_calibration_parameters"])
            parameters["cavity_transmission_w_per_v"] = new_c
            new.attrs["calibration_parameters_json"] = json.dumps(
                parameters, sort_keys=True)
            new.attrs["calibration_method"] = (
                "ape_rms_quadrature_calibration_v1_derived_c_trans_update")
            new.attrs["parent_calibration_method"] = (
                "ape_rms_quadrature_calibration_v1")
            new.attrs["parent_calibration_code_commit"] = meta["source_code_commit"]
            new.attrs["parent_hdf5_path"] = str(source)
            new.attrs["parent_hdf5_size_bytes"] = source.stat().st_size
            new.attrs["derived_i0_ratio"] = k
            new.attrs["derived_eta0_and_rho_ratio"] = factor
            new.attrs["h2_provenance_status"] = (
                "researcher working acceptance; not independently confirmed")
            new.attrs["physical_ape_coupling_approved"] = False

            for name, unit in DATASETS.items():
                prev = old[name]
                ds = new.create_dataset(
                    name, shape=prev.shape, dtype="f8",
                    chunks=(min(block, meta["samples"]),))
                for key, val in prev.attrs.items():
                    ds.attrs[key] = val
                ds.attrs["unit"] = unit
                old_digest = hashlib.sha256()
                new_digest = hashlib.sha256()
                gain = k if name == "i0_w" else factor
                for start in range(0, meta["samples"], block):
                    stop = min(start + block, meta["samples"])
                    original = np.asarray(prev[start:stop], dtype=np.float64)
                    if not np.isfinite(original).all():
                        raise ValueError(f"Nonfinite source data in {name}")
                    if name in ("i0_w", "eta0") and np.any(original <= 0):
                        raise ValueError(f"Nonpositive source data in {name}")
                    converted = original * gain
                    if not np.isfinite(converted).all():
                        raise ValueError(f"Nonfinite output in {name}")
                    ds[start:stop] = converted
                    old_digest.update(original.astype("<f8", copy=False).tobytes())
                    new_digest.update(converted.astype("<f8", copy=False).tobytes())
                hashes[name] = {
                    "parent_values_sha256_little_endian_f64": old_digest.hexdigest(),
                    "derived_values_sha256_little_endian_f64": new_digest.hexdigest(),
                }
                new.attrs[name + "_parent_values_sha256"] = old_digest.hexdigest()
                new.attrs[name + "_derived_values_sha256"] = new_digest.hexdigest()
                print("COMPLETED", name, flush=True)
            new.flush()
        # The completed HDF5 is made visible only after all four arrays finish.
        if target.exists():
            raise FileExistsError(target)
        os.link(temp, target)
        temp.unlink()
    finally:
        temp.unlink(missing_ok=True)

    report = dict(plan)
    report["status"] = "derived_calibration_created_not_physical_limit_approval"
    report["data_hashes"] = hashes
    report["source_manifest_sha256"] = meta["source_manifest_sha256"]
    with provenance.open("x") as out:
        json.dump(report, out, indent=2)
    print("NEW_CALIBRATED_HDF5_CREATED", target)
    print("PROVENANCE", provenance)
    print("NO_C_RERUN_NO_VENDOR_MODIFICATIONS")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-old-c-trans", type=float, required=True)
    parser.add_argument("--new-c-trans", type=float, required=True)
    parser.add_argument("--h2", type=float, required=True)
    parser.add_argument("--chunk-samples", type=int, default=500000)
    parser.add_argument("--execute", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
