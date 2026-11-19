"""Optional APE two-QWP calibration preparation, using the EXISTING APE optics code.

This does not infer unmeasured X/Y projections or C-LPSD kernel values.
Outputs exactly the numerical input expected by configuration_A.build_response.
The reused QWP formula is PROVISIONAL, and this output is not a validated
coupling limit or observed signal.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling
from configuration_A.build_response import ACalibrationError, _numeric_1d


def derive_provisional_qwp_calibration(source_npz, optics_json, output_npz):
    src = Path(source_npz).expanduser().resolve()
    cfg = Path(optics_json).expanduser().resolve()
    dst = Path(output_npz).expanduser().resolve()
    if not src.is_file() or not cfg.is_file() or dst.suffix != ".npz" or dst.exists():
        raise ACalibrationError("Require existing input files and a new .npz output")
    with cfg.open(encoding="utf-8") as stream:
        parameters = json.load(stream)
    expected = {"finesse_s", "tau_s", "tau1_s", "rho_dm_gev4", "unit_conversion"}
    if not isinstance(parameters, dict) or set(parameters) != expected:
        raise ACalibrationError("Optics JSON must define exactly the five existing QWP model parameters")
    with np.load(src, allow_pickle=False) as z:
        required = {
            "frequency_hz", "quadrature_amplitude_projection_x",
            "quadrature_amplitude_projection_y",
            "c_psd_hz_inv_per_unit_variance_x",
            "c_psd_hz_inv_per_unit_variance_y",
        }
        if set(z.files) != required:
            raise ACalibrationError("Kernel/projection NPZ does not match QWP interface")
        frequency = _numeric_1d(z["frequency_hz"], "frequency_hz", positive=True)
        if np.any(np.diff(frequency) <= 0):
            raise ACalibrationError("frequency grid must be increasing")
        projections = []
        kernels = []
        for channel in ("x", "y"):
            p = _numeric_1d(z[f"quadrature_amplitude_projection_{channel}"], channel)
            k = _numeric_1d(z[f"c_psd_hz_inv_per_unit_variance_{channel}"], channel)
            if p.shape != frequency.shape or k.shape != frequency.shape:
                raise ACalibrationError("Projections/kernels must match native frequency grid")
            if np.any(np.abs(p) > 1) or np.any(k < 0):
                raise ACalibrationError("Invalid quadrature projection or C kernel")
            projections.append(p)
            kernels.append(k)
    # The function and formula below are the existing APE code, not a
    # new optics derivation introduced for Configuration A.
    response = qwp_amplitude_per_coupling(frequency, **parameters)
    if response.shape != frequency.shape or not np.isfinite(response).all():
        raise ACalibrationError("Existing QWP response is invalid")
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        with dst.open("xb") as stream:
            np.savez(stream, frequency_hz=frequency,
                     peak_rotation_rad_per_gev_inv_x=response*projections[0],
                     peak_rotation_rad_per_gev_inv_y=response*projections[1],
                     c_psd_hz_inv_per_unit_variance_x=kernels[0],
                     c_psd_hz_inv_per_unit_variance_y=kernels[1])
    except BaseException:
        dst.unlink(missing_ok=True)
        raise
    return dst


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kernel-projection-npz", required=True)
    p.add_argument("--optics-json", required=True)
    p.add_argument("--output-calibration-npz", required=True)
    a = p.parse_args()
    print(derive_provisional_qwp_calibration(
        a.kernel_projection_npz, a.optics_json, a.output_calibration_npz
    ))


if __name__ == "__main__":
    main()
