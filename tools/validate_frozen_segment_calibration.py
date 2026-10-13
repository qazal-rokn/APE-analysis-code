#!/usr/bin/env python3
"""Validate APE calibration on a short slice from a frozen segment manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ape_c_lpsd.calibration import calibrate_chunk
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _load_manifest(path: Path) -> dict:
    with path.open() as handle:
        manifest = json.load(handle)
    if manifest.get("kind") != "frozen_validated_segment":
        raise ValueError("manifest is not a frozen_validated_segment")
    if float(manifest.get("sample_rate_hz", 0.0)) <= 0:
        raise ValueError("manifest sample_rate_hz must be positive")
    parts = manifest.get("parts")
    if not isinstance(parts, list) or not parts:
        raise ValueError("manifest must contain non-empty parts")
    return manifest


def _slice_parts(manifest: dict, duration_s: float) -> tuple[int, int, list[dict]]:
    fs = float(manifest["sample_rate_hz"])
    total = int(manifest["sample_count"])
    requested = int(round(duration_s * fs))
    if requested <= 0 or requested > total:
        raise ValueError("requested validation duration is outside frozen segment")

    segment_start = (total - requested) // 2
    segment_stop = segment_start + requested
    selected = []
    cursor = 0

    for part in manifest["parts"]:
        count = int(part["raw_stop"]) - int(part["raw_start"])
        part_segment_start = cursor
        part_segment_stop = cursor + count
        overlap_start = max(segment_start, part_segment_start)
        overlap_stop = min(segment_stop, part_segment_stop)
        if overlap_stop > overlap_start:
            selected.append({
                "path": part["path"],
                "raw_start": int(part["raw_start"]) + overlap_start - part_segment_start,
                "raw_stop": int(part["raw_start"]) + overlap_stop - part_segment_start,
            })
        cursor = part_segment_stop

    covered = sum(p["raw_stop"] - p["raw_start"] for p in selected)
    if covered != requested:
        raise RuntimeError(f"slice coverage mismatch: expected {requested}, got {covered}")
    return segment_start, segment_stop, selected


def _stats(values: np.ndarray) -> dict:
    finite = np.isfinite(values)
    if not np.any(finite):
        return {"count": int(values.size), "finite_fraction": 0.0}
    x = values[finite]
    return {
        "count": int(values.size),
        "finite_fraction": float(x.size / values.size),
        "mean": float(np.mean(x)),
        "std": float(np.std(x)),
        "median": float(np.median(x)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--duration-s", type=float, default=60.0)
    parser.add_argument("--chunk-size-samples", type=int, default=1_000_000)
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)
    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = _load_manifest(manifest_path)

    fs = float(manifest["sample_rate_hz"])
    manifest_sha256 = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    slice_start, slice_stop, selected_parts = _slice_parts(manifest, args.duration_s)

    i0_parts, eta_parts, rho_x_parts, rho_y_parts = [], [], [], []
    x1_parts, y1_parts, x1_rec_parts, y1_rec_parts = [], [], [], []

    q = config.calibration.responsivity_a_per_w
    gain = config.calibration.transimpedance_gain_v_per_a

    for selected in selected_parts:
        source = Path(selected["path"]).expanduser()
        for chunk in iter_ape_dxd_chunks(
            source,
            channel_names=config.channels.as_mapping(),
            wrapper_path=dewesoft.wrapper_path,
            library_path=dewesoft.library_path,
            chunk_size_samples=args.chunk_size_samples,
            start_index=int(selected["raw_start"]),
            stop_index=int(selected["raw_stop"]),
        ):
            calibrated = calibrate_chunk(chunk, config.calibration)
            denominator = np.sqrt(2.0) * q * gain * calibrated.i0_w * calibrated.eta0

            i0_parts.append(calibrated.i0_w)
            eta_parts.append(calibrated.eta0)
            rho_x_parts.append(calibrated.rho_x)
            rho_y_parts.append(calibrated.rho_y)
            x1_parts.append(np.asarray(chunk.x_1f, dtype=np.float64))
            y1_parts.append(np.asarray(chunk.y_1f, dtype=np.float64))
            x1_rec_parts.append(denominator * calibrated.rho_x)
            y1_rec_parts.append(denominator * calibrated.rho_y)

    i0 = np.concatenate(i0_parts)
    eta0 = np.concatenate(eta_parts)
    rho_x = np.concatenate(rho_x_parts)
    rho_y = np.concatenate(rho_y_parts)
    x1 = np.concatenate(x1_parts)
    y1 = np.concatenate(y1_parts)
    x1_rec = np.concatenate(x1_rec_parts)
    y1_rec = np.concatenate(y1_rec_parts)

    expected = slice_stop - slice_start
    if i0.size != expected:
        raise RuntimeError(f"sample count mismatch: expected {expected}, got {i0.size}")

    finite_all = np.isfinite(i0) & np.isfinite(eta0) & np.isfinite(rho_x) & np.isfinite(rho_y)
    valid_x = np.isfinite(x1_rec) & np.isfinite(x1)
    valid_y = np.isfinite(y1_rec) & np.isfinite(y1)
    xerr = x1_rec[valid_x] - x1[valid_x]
    yerr = y1_rec[valid_y] - y1[valid_y]

    summary = {
        "method": "ape_frozen_segment_calibration_validation_v1",
        "writes_output_files": False,
        "manifest": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "duration_s": args.duration_s,
        "sample_rate_hz": fs,
        "sample_count": int(expected),
        "segment_relative_start_index": int(slice_start),
        "segment_relative_stop_index_exclusive": int(slice_stop),
        "selected_parts": selected_parts,
        "sanity": {
            "all_outputs_finite_fraction": float(np.mean(finite_all)),
            "i0_positive_fraction": float(np.mean(i0[np.isfinite(i0)] > 0.0)) if np.any(np.isfinite(i0)) else 0.0,
            "eta0_positive_fraction": float(np.mean(eta0[np.isfinite(eta0)] > 0.0)) if np.any(np.isfinite(eta0)) else 0.0,
        },
        "statistics": {
            "i0_w": _stats(i0),
            "eta0": _stats(eta0),
            "rho_x": _stats(rho_x),
            "rho_y": _stats(rho_y),
        },
        "inverse_reconstruction": {
            "x1f_max_abs_error_v": float(np.max(np.abs(xerr))) if xerr.size else None,
            "x1f_rms_error_v": float(np.sqrt(np.mean(xerr**2))) if xerr.size else None,
            "y1f_max_abs_error_v": float(np.max(np.abs(yerr))) if yerr.size else None,
            "y1f_rms_error_v": float(np.sqrt(np.mean(yerr**2))) if yerr.size else None,
        },
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
