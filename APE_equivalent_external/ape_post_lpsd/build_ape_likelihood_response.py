"""APE-only bridge from existing g_a-gamma-gamma/QWP PSD model to vendor response NPZ.

No change to original C, vendor likelihood, background fitter or APE physics.
Inputs must be supplied explicitly: an existing native C D/PSD kernel and
physical response configuration. Output remains PROVISIONAL.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from .ape_axion_signal_psd import save_provisional_response


class ResponseBuildError(ValueError):
    """An input violates the previously established APE response contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_from_existing_kernel(*, spectrum_hdf5, kernel_npz, physics_json,
                               output_npz, require_full_grid=True):
    """Only connect existing measured C kernel and approved QWP inputs.

    No rebinning, replacement halo, PSD/PS conversion guess, or vendor edits.
    Full native frequency coverage is mandatory unless explicitly disabled
    for a diagnostic-only partial-kernel test.
    """
    spectrum = Path(spectrum_hdf5).expanduser().resolve()
    kernel_path = Path(kernel_npz).expanduser().resolve()
    physics_path = Path(physics_json).expanduser().resolve()
    output = Path(output_npz).expanduser().resolve()
    for path in (spectrum, kernel_path, physics_path):
        if not path.is_file():
            raise ResponseBuildError(f"Missing input: {path}")
    if output.suffix != ".npz" or output.exists():
        raise ResponseBuildError("Output must be a new NPZ file; no overwrites")

    try:
        with h5py.File(spectrum, "r") as source:
            native = np.asarray(source["frequency"][:], dtype=np.float64)
        with np.load(kernel_path, allow_pickle=False) as source:
            required = {"native_frequency_hz", "output_bin_indices",
                        "halo_frequency_hz",
                        "native_c_psd_kernel_per_variance_hz", "metadata_json"}
            if set(source.files) != required:
                raise ResponseBuildError(f"Expected exact kernel fields: {sorted(required)}")
            kernel_grid = np.asarray(source["native_frequency_hz"], dtype=np.float64)
            indices = np.asarray(source["output_bin_indices"])
            modes = np.asarray(source["halo_frequency_hz"], dtype=np.float64)
            kernel = np.asarray(source["native_c_psd_kernel_per_variance_hz"],
                                dtype=np.float64)
            raw = source["metadata_json"]
            if raw.ndim != 0:
                raise ResponseBuildError("Kernel metadata must be a scalar JSON string")
            kernel_meta = json.loads(str(raw.item()))
        physics = json.loads(physics_path.read_text(encoding="utf-8"))
    except (OSError, KeyError, ValueError, TypeError, json.JSONDecodeError) as exc:
        if isinstance(exc, ResponseBuildError):
            raise
        raise ResponseBuildError(f"Cannot read response inputs: {exc}") from exc

    if not np.array_equal(native, kernel_grid):
        raise ResponseBuildError("Kernel and observed C frequency grids differ")
    if not isinstance(kernel_meta, dict) or (
        kernel_meta.get("kernel_units")
        != "Hz^-1 per unit time-domain rad^2 variance"
    ) or kernel_meta.get("source_format") != "native C method-2 D = PSD rad^2/Hz":
        raise ResponseBuildError("Missing native C D/PSD normalization provenance")
    if not isinstance(kernel_meta.get("c_configuration_source"), str) or not kernel_meta["c_configuration_source"].strip():
        raise ResponseBuildError("C configuration provenance is missing")
    if indices.ndim != 1 or indices.dtype.kind not in "iu" or not len(indices):
        raise ResponseBuildError("Kernel indices must be nonempty integers")
    if require_full_grid and not np.array_equal(indices, np.arange(native.size)):
        raise ResponseBuildError("Full likelihood requires kernel response at every native bin")

    fields = {"normalized_halo_weights", "quadrature_power_fraction",
              "reference_coupling_gev_inv", "finesse_s", "tau_s", "tau1_s",
              "rho_dm_gev4", "unit_conversion", "halo_source",
              "calibration_version"}
    if not isinstance(physics, dict) or set(physics) != fields:
        raise ResponseBuildError(f"Physics JSON requires exactly {sorted(fields)}")
    for key in ("halo_source", "calibration_version"):
        if not isinstance(physics[key], str) or not physics[key].strip():
            raise ResponseBuildError(f"Missing physical provenance: {key}")

    info = save_provisional_response(
        output,
        native_frequency_hz=native,
        output_bin_indices=indices,
        halo_frequency_hz=modes,
        normalized_halo_weights=physics["normalized_halo_weights"],
        native_c_psd_kernel_per_variance_hz=kernel,
        quadrature_power_fraction=physics["quadrature_power_fraction"],
        reference_coupling_gev_inv=physics["reference_coupling_gev_inv"],
        finesse_s=physics["finesse_s"],
        tau_s=physics["tau_s"],
        tau1_s=physics["tau1_s"],
        rho_dm_gev4=physics["rho_dm_gev4"],
        unit_conversion=physics["unit_conversion"],
        halo_source=physics["halo_source"],
        native_c_kernel_source=f"{kernel_path} sha256={_sha256(kernel_path)}",
        calibration_version=physics["calibration_version"],
    )
    return {
        "response_npz": str(output),
        "native_bins": int(native.size),
        "modeled_bins": int(indices.size),
        "halo_modes": int(modes.size),
        "full_grid_required": bool(require_full_grid),
        "physics_status": info["optical_response_status"],
        "spectrum_sha256": _sha256(spectrum),
        "kernel_sha256": _sha256(kernel_path),
        "physics_sha256": _sha256(physics_path),
        "response_sha256": _sha256(output),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum-hdf5", required=True)
    parser.add_argument("--kernel-npz", required=True)
    parser.add_argument("--physics-json", required=True)
    parser.add_argument("--output-npz", required=True)
    parser.add_argument("--allow-partial-grid", action="store_true",
                        help="Diagnostic only: accept an incomplete kernel")
    options = parser.parse_args()
    report = build_from_existing_kernel(
        spectrum_hdf5=options.spectrum_hdf5,
        kernel_npz=options.kernel_npz,
        physics_json=options.physics_json,
        output_npz=options.output_npz,
        require_full_grid=not options.allow_partial_grid,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
