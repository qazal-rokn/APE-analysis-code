"""Prepare ONE vendor-readable APE spectrum by summing separate X/Y C-LPSD power.

EXPERIMENTAL, UNTESTED. Does not modify original vendor functions or inputs.
Input must be linear POWER per C-LPSD bin, not PSD per Hz. No hidden
density-to-power conversion is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def _load(path):
    path = Path(path).expanduser().resolve(strict=True)
    arr = np.loadtxt(path, comments="#", ndmin=2)
    if arr.ndim != 2 or arr.shape[1] != 2 or len(arr) < 2:
        raise ValueError(f"Expected two-column frequency,PS: {path}")
    f, ps = arr[:, 0], arr[:, 1]
    if (not np.isfinite(arr).all() or np.any(f <= 0) or
            np.any(np.diff(f) <= 0) or np.any(ps <= 0)):
        raise ValueError(f"Invalid native power spectrum: {path}")
    return f, ps, path


def prepare(*, rho_x_ps, rho_y_ps, output_h5, fs_hz,
            t0, t1, provenance_json=None):
    fx, x, px = _load(rho_x_ps)
    fy, y, py = _load(rho_y_ps)
    if not np.array_equal(fx, fy):
        raise ValueError("X and Y C-LPSD grids differ: no interpolation permitted")
    if not np.isfinite(fs_hz) or fs_hz <= 0:
        raise ValueError("Invalid original sample rate")
    if int(t1) <= int(t0):
        raise ValueError("Invalid observation bounds")
    total = x + y
    if not np.isfinite(total).all() or np.any(total <= 0):
        raise ValueError("Invalid combined linear power")
    target = Path(output_h5).expanduser().resolve()
    if target.exists() or not target.parent.is_dir():
        raise FileExistsError("Output exists or parent directory missing")
    from scalar_dark_matter_lpsd.LPSDIO import LPSDData
    native = LPSDData(
        np.log(total), freq=fx, fmin=float(fx[0]), fmax=float(fx[-1]),
        fs=float(fs_hz), Jdes=int(len(fx)))
    native.save_to_hdf5(str(target))
    def digest(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = {
        "status": "EXPERIMENTAL_UNTESTED_NOT_VALIDATED",
        "input_unit": "native_C_LPSD_linear_power_per_bin",
        "operation": "PS_total[j] = PS_X[j] + PS_Y[j]",
        "quadratures_combined_before_background_and_likelihood": True,
        "frequency_grid_exactly_equal": True,
        "x_file": str(px), "y_file": str(py),
        "x_sha256": digest(px), "y_sha256": digest(py),
        "output_h5": str(target), "output_sha256": digest(target),
        "t0": int(t0), "t1": int(t1),
        "note": "Vendor group filenames must end _t0_t1_APE.h5; "
                "APE response replaces vendor detector beta downstream. "
                "Sum noise distribution and response normalization not validated.",
    }
    if provenance_json is not None:
        path = Path(provenance_json).expanduser().resolve()
        if path.exists():
            raise FileExistsError(path)
        path.write_text(json.dumps(manifest, indent=2))
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rho-x-ps", required=True)
    p.add_argument("--rho-y-ps", required=True)
    p.add_argument("--output-h5", required=True)
    p.add_argument("--fs-hz", type=float, required=True)
    p.add_argument("--t0", type=int, required=True)
    p.add_argument("--t1", type=int, required=True)
    p.add_argument("--provenance-json")
    args = vars(p.parse_args())
    print(json.dumps(prepare(**{k.replace("-", "_"): v for k,v in args.items()}), indent=2))


if __name__ == "__main__":
    main()
