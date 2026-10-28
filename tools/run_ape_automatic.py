"""Automatic APE frozen-segment calibration stage.

ONE run configuration -> immutable DXD slices -> calibrated HDF5 -> handoff.
No manual intermediate file construction or silent scientific changes.
Dry-run by default. Execution requires --execute on macOS for DWDataReader.
Native C and vendor likelihood are never executed by this stage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import yaml


REQUIRED_CALIBRATION = (
    "cavity_transmission_w_per_v",
    "responsivity_a_per_w",
    "transimpedance_gain_v_per_a",
    "second_harmonic_pd_response",
)
UNITS = {"i0_w": "W", "eta0": "dimensionless",
         "rho_x": "rad", "rho_y": "rad"}
MIN_FREE_MARGIN = 512 * 1024**2


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve(value: str, base: Path) -> Path:
    p = Path(value).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def _read_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as file:
        obj = yaml.safe_load(file)
    if not isinstance(obj, dict):
        raise ValueError(f"Expected mapping: {path}")
    return obj


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as file:
        obj = json.load(file)
    if not isinstance(obj, dict):
        raise ValueError(f"Expected object: {path}")
    return obj


def _positive_number(data: dict, key: str) -> float:
    x = data.get(key)
    if isinstance(x, bool):
        raise ValueError(f"Invalid calibration parameter {key}")
    try:
        f = float(x)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Missing calibration parameter {key}") from exc
    if not (math.isfinite(f) and f > 0):
        raise ValueError(f"Invalid calibration parameter {key}")
    return f


def prepare(config_path: Path, *, require_mac: bool = False):
    configuration = _read_yaml(config_path)
    expected = {
        "source_root", "frozen_segment_manifest", "dxd_ingestion_config",
        "output_root", "run_id", "calibration",
    }
    if set(configuration) != expected:
        raise ValueError(
            "Run YAML must contain exactly: " + ", ".join(sorted(expected))
        )
    base = config_path.parent
    root = _resolve(configuration["source_root"], base)
    manifest_path = _resolve(configuration["frozen_segment_manifest"], base)
    dxd_path = _resolve(configuration["dxd_ingestion_config"], base)
    output_root = _resolve(configuration["output_root"], base)
    if not root.is_dir() or not output_root.is_dir():
        raise FileNotFoundError("source_root and output_root must already exist")
    if not manifest_path.is_file() or not dxd_path.is_file():
        raise FileNotFoundError("Frozen manifest or DXD ingestion YAML missing")
    run_id = configuration["run_id"]
    if (not isinstance(run_id, str) or not run_id.startswith("ape-run-")
            or not run_id.replace("-", "").replace("_", "").isalnum()
            or len(run_id) > 96):
        raise ValueError("run_id must be a safe, unique ape-run-* name")
    out = (output_root / run_id).resolve()
    if out.exists() or out == root or out == output_root:
        raise FileExistsError(f"Non-overwrite policy: {out}")
    if root in out.parents or out in root.parents:
        raise ValueError("Outputs must be outside raw source root")
    cal = configuration["calibration"]
    if not isinstance(cal, dict) or set(cal) != set(REQUIRED_CALIBRATION):
        raise ValueError("Four explicit calibration parameters required")
    approved = {key: _positive_number(cal, key)
                for key in REQUIRED_CALIBRATION}

    manifest = _read_json(manifest_path)
    if manifest.get("kind") != "frozen_validated_segment":
        raise ValueError("Expected frozen_validated_segment")
    if manifest.get("unit") != "V":
        raise ValueError("Frozen source must be raw volts")
    parts = manifest.get("parts")
    if not isinstance(parts, list) or not parts:
        raise ValueError("Frozen source parts missing")
    sample_count = int(manifest["sample_count"])
    fs = float(manifest["sample_rate_hz"])
    if not (sample_count > 0 and math.isfinite(fs) and fs > 0):
        raise ValueError("Invalid source sample count/rate")
    total = 0
    files = []
    for part in parts:
        rel = Path(part["path"])
        if rel.is_absolute():
            source = rel.resolve()
        else:
            source = (root / rel).resolve()
            if source != root and root not in source.parents:
                raise ValueError("Relative raw source escapes declared root")
        if not source.is_file():
            raise FileNotFoundError(f"Raw DXD source absent: {source}")
        first, last = int(part["raw_start"]), int(part["raw_stop"])
        if first < 0 or last <= first:
            raise ValueError("Invalid raw sample range")
        if "sample_count" in part and int(part["sample_count"]) != last - first:
            raise ValueError("Frozen source-part sample count mismatch")
        total += last - first
        files.append({"path": str(source), "raw_start": first,
                      "raw_stop": last})
    if total != sample_count:
        raise ValueError("Frozen segment parts do not sum to sample_count")
    digest = _sha256(manifest_path)
    provenance_path = manifest_path.with_suffix(".sha256.json")
    if not provenance_path.is_file() or _read_json(
        provenance_path).get("sha256") != digest:
        raise ValueError("Frozen-segment manifest provenance hash mismatch")

    dxd = _read_yaml(dxd_path)
    for field in ("channels", "cleaning", "reference_project"):
        if field not in dxd:
            raise ValueError(f"DXD config missing {field}")
    ref = _resolve(dxd["reference_project"], dxd_path.parent)
    wrapper = ref / "vendor/dewesoft/5.0.8/DWDataReaderHeader.py"
    library = ref / "vendor/dewesoft/5.0.8/DWDataReaderLib64.dylib"
    if not wrapper.is_file() or not library.is_file():
        raise FileNotFoundError("Original macOS Dewesoft library/wrapper missing")
    if require_mac and sys.platform != "darwin":
        raise RuntimeError(
            "DXD calibration runs on macOS with the original .dylib, "
            "not inside Linux Docker"
        )

    # The historical input YAML may contain C_trans=0.020. Never edit it.
    # Build a run-specific snapshot with the researcher's confirmed values.
    snapshot = dict(dxd)
    snapshot["reference_project"] = str(ref)
    snapshot["calibration"] = approved
    snapshot["output_directory"] = str(out)
    nominal_bytes = sample_count * len(UNITS) * 8
    free_bytes = shutil.disk_usage(output_root).free
    plan = {
        "status": "plan_only",
        "run_id": run_id,
        "platform": platform.system(),
        "source_root": str(root),
        "frozen_manifest": str(manifest_path),
        "frozen_manifest_sha256": digest,
        "dxd_config": str(dxd_path),
        "dxd_config_sha256": _sha256(dxd_path),
        "raw_parts": files,
        "sample_count": sample_count,
        "sample_rate_hz": fs,
        "calibration": approved,
        "calibration_method": "ape_rms_quadrature_calibration_v1",
        "output_hdf5": str(out / "calibrated_ape.h5"),
        "estimated_uncompressed_bytes": nominal_bytes,
        "free_bytes": free_bytes,
        "scientific_approval": {
            "H2_measurement_independently_verified": False,
            "physical_XY_axion_response_validated": False,
            "real_data_coupling_limits_authorized": False,
        },
        "native_c_executed": False,
        "original_vendor_modified": False,
        "source_data_modified": False,
        "next_stage": "native_C_PS_handoff_requires_separate_approval",
    }
    if nominal_bytes + MIN_FREE_MARGIN > free_bytes:
        raise RuntimeError("Insufficient disk space for new calibrated output")
    return plan, snapshot, out, root, manifest_path, dxd_path


def execute(config_path: Path, *, execute_calibration: bool):
    plan, snapshot, out, root, manifest, dxd_path = prepare(
        config_path, require_mac=execute_calibration
    )
    print(json.dumps(plan, indent=2), flush=True)
    if not execute_calibration:
        print("APE_AUTOMATIC_CALIBRATION_PLAN_OK_NO_DATA_WRITTEN")
        return
    # Hard safety guard: output is a new, dedicated run folder.
    out.mkdir(parents=False, exist_ok=False)
    (out / "config_snapshot.yaml").write_text(
        yaml.safe_dump(snapshot, sort_keys=True), encoding="utf-8")
    (out / "plan.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    tool = Path(__file__).resolve().parent / "calibrate_frozen_segment_to_hdf5.py"
    output = out / "calibrated_ape.h5"
    command = [
        sys.executable, str(tool), "--config", str(out / "config_snapshot.yaml"),
        "--manifest", str(manifest), "--output", str(output),
    ]
    print("CALIBRATING_RAW_DXD_AUTOMATICALLY", flush=True)
    try:
        subprocess.run(command, cwd=root, check=True)
        with h5py.File(output, "r") as h5:
            n, fs = int(h5.attrs["sample_count"]), float(h5.attrs["sample_rate_hz"])
            if n != plan["sample_count"] or fs != plan["sample_rate_hz"]:
                raise RuntimeError("Calibrated HDF5 rate/count mismatch")
            if str(h5.attrs["manifest_sha256"]) != plan["frozen_manifest_sha256"]:
                raise RuntimeError("Output does not match frozen source manifest")
            recorded = json.loads(h5.attrs["calibration_parameters_json"])
            if recorded != plan["calibration"]:
                raise RuntimeError("Calibrated constants do not match run contract")
            if h5.attrs["x_y_combined"]:
                raise RuntimeError("Unexpected combined quadratures")
            for name, unit in UNITS.items():
                ds = h5[name]
                if (ds.shape != (n,) or ds.dtype.name != "float64"
                        or ds.attrs.get("unit") != unit):
                    raise RuntimeError(f"Invalid calibrated output {name}")
        outcome = {
            "status": "calibration_complete_C_backend_pending_scientific_gate",
            "run_id": plan["run_id"],
            "calibrated_hdf5": str(output),
            "calibrated_hdf5_sha256": _sha256(output),
            "frozen_manifest_sha256": plan["frozen_manifest_sha256"],
            "calibration": plan["calibration"],
            "sample_count": plan["sample_count"],
            "sample_rate_hz": plan["sample_rate_hz"],
            "next_stage": "transport_or_mount_HDF5_to_verified_native_C_backend",
            "native_c_executed": False,
            "physical_g_limits_authorized": False,
        }
        (out / "stage_result.json").write_text(
            json.dumps(outcome, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(outcome, indent=2), flush=True)
        print("APE_AUTOMATIC_CALIBRATION_COMPLETE")
    except BaseException:
        # Preserve an incomplete run as forensic evidence. Never overwrite
        # it or silently resume: each attempt requires a fresh run_id.
        print(f"APE_AUTOMATIC_CALIBRATION_FAILED_OUTPUT_RETAINED: {out}",
              file=sys.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-config", type=Path, required=True)
    parser.add_argument("--execute-calibration", action="store_true",
                        help="Explicit opt-in on macOS; dry-run otherwise")
    args = parser.parse_args()
    execute(args.run_config.expanduser().resolve(strict=True),
            execute_calibration=args.execute_calibration)


if __name__ == "__main__":
    main()
