"""Chunked HDF5 output for the frozen calibrated APE segment."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError as exc:  # pragma: no cover - exercised by runtime environment
    raise ImportError(
        "h5py is required for calibrated HDF5 output. Install it in the project venv."
    ) from exc


_DATASETS = {
    "i0_w": "W",
    "eta0": "dimensionless",
    "rho_x": "rad",
    "rho_y": "rad",
}


def create_calibrated_hdf5(
    path: str | Path,
    *,
    sample_count: int,
    sample_rate_hz: float,
    manifest_sha256: str,
    calibration_method: str,
    calibration_parameters: dict[str, float],
    source_parts: list[dict],
    code_commit: str,
    chunk_samples: int = 1_000_000,
):
    """Create an empty calibrated HDF5 product and return the open file.

    The caller owns the returned file handle and must close it.
    """
    path = Path(path)
    if sample_count <= 0:
        raise ValueError("sample_count must be positive")
    if sample_rate_hz <= 0:
        raise ValueError("sample_rate_hz must be positive")
    if chunk_samples <= 0:
        raise ValueError("chunk_samples must be positive")

    h5 = h5py.File(path, "x")
    h5.attrs["product"] = "ape_frozen_calibrated_segment"
    h5.attrs["format_version"] = "1"
    h5.attrs["sample_count"] = int(sample_count)
    h5.attrs["sample_rate_hz"] = float(sample_rate_hz)
    h5.attrs["manifest_sha256"] = str(manifest_sha256)
    h5.attrs["calibration_method"] = str(calibration_method)
    h5.attrs["calibration_parameters_json"] = json.dumps(
        calibration_parameters, sort_keys=True
    )
    h5.attrs["source_parts_json"] = json.dumps(source_parts, sort_keys=True)
    h5.attrs["code_commit"] = str(code_commit)
    h5.attrs["x_y_combined"] = False
    h5.attrs["filtering_applied"] = False
    h5.attrs["detrending_applied"] = False
    h5.attrs["resampling_applied"] = False
    h5.attrs["interpolation_applied"] = False

    chunks = (min(int(chunk_samples), int(sample_count)),)
    for name, unit in _DATASETS.items():
        ds = h5.create_dataset(
            name,
            shape=(sample_count,),
            dtype=np.float64,
            chunks=chunks,
        )
        ds.attrs["unit"] = unit

    return h5


def write_calibrated_chunk(
    h5,
    *,
    output_start: int,
    i0_w: np.ndarray,
    eta0: np.ndarray,
    rho_x: np.ndarray,
    rho_y: np.ndarray,
) -> int:
    """Write one aligned calibrated chunk and return the next output offset."""
    arrays = {
        "i0_w": np.asarray(i0_w, dtype=np.float64),
        "eta0": np.asarray(eta0, dtype=np.float64),
        "rho_x": np.asarray(rho_x, dtype=np.float64),
        "rho_y": np.asarray(rho_y, dtype=np.float64),
    }
    shapes = {name: values.shape for name, values in arrays.items()}
    if len(set(shapes.values())) != 1:
        raise ValueError(f"calibrated arrays are misaligned: {shapes}")

    shape = next(iter(shapes.values()))
    if len(shape) != 1 or shape[0] <= 0:
        raise ValueError("calibrated chunk arrays must be non-empty and one-dimensional")

    n = int(shape[0])
    output_stop = int(output_start) + n
    if output_start < 0 or output_stop > int(h5.attrs["sample_count"]):
        raise ValueError("calibrated chunk lies outside HDF5 output bounds")

    for name, values in arrays.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} contains non-finite calibrated samples")
        h5[name][output_start:output_stop] = values

    return output_stop


def validate_completed_hdf5(h5, expected_samples: int) -> None:
    """Validate shape and finite endpoints before the product is finalized."""
    if int(h5.attrs["sample_count"]) != int(expected_samples):
        raise ValueError("HDF5 sample_count metadata mismatch")
    for name in _DATASETS:
        ds = h5[name]
        if ds.shape != (expected_samples,):
            raise ValueError(f"{name} has wrong final shape {ds.shape}")
        if expected_samples:
            endpoints = np.asarray([ds[0], ds[-1]], dtype=np.float64)
            if not np.all(np.isfinite(endpoints)):
                raise ValueError(f"{name} has non-finite endpoint samples")
