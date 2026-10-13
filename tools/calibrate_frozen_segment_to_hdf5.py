#!/usr/bin/env python3
"""Calibrate one frozen APE segment and write the canonical HDF5 product."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ape_c_lpsd.calibrated_hdf5 import (
    create_calibrated_hdf5,
    validate_completed_hdf5,
    write_calibrated_chunk,
)
from ape_c_lpsd.calibration import calibrate_chunk
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _load_frozen_manifest(path: Path) -> dict:
    with path.open() as handle:
        manifest = json.load(handle)
    if manifest.get("kind") != "frozen_validated_segment":
        raise ValueError("manifest is not a frozen_validated_segment")
    if manifest.get("unit") != "V":
        raise ValueError("frozen input unit must be V")
    parts = manifest.get("parts")
    if not isinstance(parts, list) or not parts:
        raise ValueError("manifest must contain non-empty parts")
    return manifest


def _verify_manifest_provenance(path: Path, digest: str) -> Path:
    provenance = path.with_suffix(".sha256.json")
    if not provenance.is_file():
        raise FileNotFoundError(f"frozen manifest provenance file not found: {provenance}")
    with provenance.open() as handle:
        payload = json.load(handle)
    expected = payload.get("sha256")
    if expected != digest:
        raise ValueError(
            "frozen manifest SHA-256 does not match its provenance file"
        )
    return provenance


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calibrate a frozen APE segment into canonical HDF5."
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument("--manifest", required=True, help="frozen segment JSON")
    parser.add_argument("--output", required=True, help="new calibrated HDF5 path")
    parser.add_argument(
        "--chunk-size-samples",
        type=int,
        default=1_000_000,
        help="raw samples per chunk (default: 1000000)",
    )
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = _load_frozen_manifest(manifest_path)
    manifest_sha256 = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    provenance_path = _verify_manifest_provenance(manifest_path, manifest_sha256)

    sample_rate_hz = float(manifest["sample_rate_hz"])
    sample_count = int(manifest["sample_count"])
    if sample_rate_hz <= 0 or sample_count <= 0:
        raise ValueError("frozen manifest has invalid rate or sample count")

    expected_from_parts = sum(
        int(part["raw_stop"]) - int(part["raw_start"])
        for part in manifest["parts"]
    )
    if expected_from_parts != sample_count:
        raise ValueError(
            f"manifest sample_count mismatch: {sample_count} vs {expected_from_parts}"
        )

    output = Path(args.output).expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(
        output.name + f".partial-{uuid.uuid4().hex}"
    )

    calibration_parameters = {
        "cavity_transmission_w_per_v":
            float(config.calibration.cavity_transmission_w_per_v),
        "responsivity_a_per_w":
            float(config.calibration.responsivity_a_per_w),
        "transimpedance_gain_v_per_a":
            float(config.calibration.transimpedance_gain_v_per_a),
        "second_harmonic_pd_response":
            float(config.calibration.second_harmonic_pd_response),
    }

    h5 = None
    output_offset = 0
    try:
        h5 = create_calibrated_hdf5(
            temporary,
            sample_count=sample_count,
            sample_rate_hz=sample_rate_hz,
            manifest_sha256=manifest_sha256,
            calibration_method="ape_rms_quadrature_calibration_v1",
            calibration_parameters=calibration_parameters,
            source_parts=manifest["parts"],
            code_commit=_git_commit(),
            chunk_samples=args.chunk_size_samples,
        )
        h5.attrs["frozen_manifest_path"] = str(manifest_path)
        h5.attrs["frozen_manifest_provenance_path"] = str(provenance_path)

        for part_number, part in enumerate(manifest["parts"], start=1):
            source = Path(part["path"]).expanduser()
            raw_start = int(part["raw_start"])
            raw_stop = int(part["raw_stop"])

            print(
                f"Calibrating part {part_number}/{len(manifest['parts'])}: "
                f"{source.name} [{raw_start}, {raw_stop})",
                file=sys.stderr,
                flush=True,
            )

            for chunk in iter_ape_dxd_chunks(
                source,
                channel_names=config.channels.as_mapping(),
                wrapper_path=dewesoft.wrapper_path,
                library_path=dewesoft.library_path,
                chunk_size_samples=args.chunk_size_samples,
                start_index=raw_start,
                stop_index=raw_stop,
            ):
                if chunk.sample_rate_hz != sample_rate_hz:
                    raise RuntimeError(
                        f"sample-rate mismatch in {source}: "
                        f"{chunk.sample_rate_hz} vs {sample_rate_hz}"
                    )

                calibrated = calibrate_chunk(chunk, config.calibration)
                if not np.all(calibrated.valid_mask):
                    bad = int(np.count_nonzero(~calibrated.valid_mask))
                    raise RuntimeError(
                        f"calibration produced {bad} invalid samples in "
                        f"{source.name} [{chunk.start_index}, {chunk.stop_index})"
                    )

                output_offset = write_calibrated_chunk(
                    h5,
                    output_start=output_offset,
                    i0_w=calibrated.i0_w,
                    eta0=calibrated.eta0,
                    rho_x=calibrated.rho_x,
                    rho_y=calibrated.rho_y,
                )

        if output_offset != sample_count:
            raise RuntimeError(
                f"calibrated output length mismatch: {output_offset} vs {sample_count}"
            )

        validate_completed_hdf5(h5, sample_count)
        h5.flush()
        h5.close()
        h5 = None

        os.replace(temporary, output)
    except Exception:
        if h5 is not None:
            h5.close()
        if temporary.exists():
            temporary.unlink()
        raise

    size_gib = output.stat().st_size / (1024 ** 3)
    summary = {
        "product": "ape_frozen_calibrated_segment",
        "output": str(output),
        "sample_count": sample_count,
        "sample_rate_hz": sample_rate_hz,
        "duration_s": sample_count / sample_rate_hz,
        "duration_h": sample_count / sample_rate_hz / 3600.0,
        "manifest_sha256": manifest_sha256,
        "calibration_method": "ape_rms_quadrature_calibration_v1",
        "datasets": ["i0_w", "eta0", "rho_x", "rho_y"],
        "dtype": "float64",
        "file_size_gib": size_gib,
        "all_calibrated_samples_valid": True,
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
