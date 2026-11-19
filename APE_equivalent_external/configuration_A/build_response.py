"""Configuration A APE-only peak_norm producer for ORIGINAL single-bin likelihood.

No dark-matter halo template, multi-bin treatment, new statistics, or changes
to the original vendor library. The caller must supply calibrated *peak*
polarization-angle response per coupling and the native C LPSD D/PSD output
per unit time-domain angular variance. No model or numerical constants are
guessed. A response file is NOT scientific certification of these inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


class ACalibrationError(ValueError):
    """Invalid/incomplete Configuration A calibration."""


def _numeric_1d(value, name, *, positive=False, allow_zero=True):
    a = np.asarray(value, dtype=np.float64)
    if a.ndim != 1 or not a.size or not np.isfinite(a).all():
        raise ACalibrationError(f"{name} must be a nonempty finite 1D array")
    if positive and np.any(a <= 0):
        raise ACalibrationError(f"{name} must be strictly positive")
    if not allow_zero and np.any(a == 0):
        raise ACalibrationError(f"{name} contains zero")
    return a


def make_singlebin_peak_norm(
    frequency_hz, *,
    peak_rotation_rad_per_gev_inv_x, peak_rotation_rad_per_gev_inv_y,
    c_psd_hz_inv_per_unit_variance_x, c_psd_hz_inv_per_unit_variance_y,
    reference_coupling_gev_inv,
):
    """Compute X/Y additive PSD [rad^2/Hz] per mu=(g/gref)^2.

    R_k(f) is measured/calibrated PEAK angular amplitude (rad / GeV^-1),
    which includes all APE optical/cavity/PEM readout and X/Y projections.
    K_k(f) is original native C method-2 output D, measured in Hz^-1 for
    unit angular *time-domain variance*, for the chosen single test frequency.
    Therefore: beta_k(f) = 0.5 * (gref * abs(R_k(f)))^2 * K_k(f).
    No extra linewidth factor is applied: the empirical K must represent
    the same single-bin signal-power convention chosen for Config A.
    """
    f = _numeric_1d(frequency_hz, "frequency_hz", positive=True)
    if np.any(np.diff(f) <= 0):
        raise ACalibrationError("frequency_hz must be strictly increasing")
    gref = float(reference_coupling_gev_inv)
    if not np.isfinite(gref) or gref <= 0:
        raise ACalibrationError("Reference coupling must be finite and positive")

    result = []
    for channel, r, k in (
        ("X", peak_rotation_rad_per_gev_inv_x, c_psd_hz_inv_per_unit_variance_x),
        ("Y", peak_rotation_rad_per_gev_inv_y, c_psd_hz_inv_per_unit_variance_y),
    ):
        response = _numeric_1d(r, f"peak rotation {channel}")
        kernel = _numeric_1d(k, f"C PSD/variance kernel {channel}")
        if response.shape != f.shape or kernel.shape != f.shape:
            raise ACalibrationError(f"{channel} response and kernel must match native grid")
        if np.any(kernel < 0):
            raise ACalibrationError(f"{channel} kernel must be nonnegative")
        beta = 0.5 * np.square(gref * response) * kernel
        if not np.isfinite(beta).all() or np.any(beta < 0):
            raise ACalibrationError(f"{channel} beta is invalid")
        result.append(beta)
    if not any(np.any(beta > 0) for beta in result):
        raise ACalibrationError("Both channel peak normalizations vanish")
    return tuple(result)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_response_file(calibration_npz, provenance_json, output_npz):
    """Build existing Step-3/4 compatible four-field response NPZ, without overwrites.

    The calibration NPZ must have EXACTLY these arrays:
      frequency_hz,
      peak_rotation_rad_per_gev_inv_x/y,
      c_psd_hz_inv_per_unit_variance_x/y.
    The JSON supplies measurement provenance and gref. Input data are not
    deemed scientifically certified merely because they pass numeric checks.
    """
    src = Path(calibration_npz).expanduser().resolve()
    info = Path(provenance_json).expanduser().resolve()
    dst = Path(output_npz).expanduser().resolve()
    if dst.suffix != ".npz" or dst.exists() or not src.is_file() or not info.is_file():
        raise ACalibrationError("Require existing inputs and a NEW .npz output")
    with info.open(encoding="utf-8") as stream:
        meta = json.load(stream)
    required_strings = ("response_source", "calibration_version",
                        "native_c_kernel_source", "single_bin_convention",
                        "peak_amplitude_calibration_source",
                        "quadrature_projection_source")
    if not isinstance(meta, dict) or any(
        not isinstance(meta.get(key), str) or not meta[key].strip()
        for key in required_strings
    ):
        raise ACalibrationError("Provenance JSON needs all six nonempty source fields")
    if meta.get("response_status") not in (
        "PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS",
        "CALIBRATION_REVIEWED_AWAITING_INJECTION_COVERAGE",
    ):
        raise ACalibrationError("Response status must explicitly remain unvalidated")
    gref = meta.get("reference_coupling_gev_inv")
    if not isinstance(gref, (float, int)) or not np.isfinite(gref) or gref <= 0:
        raise ACalibrationError("reference_coupling_gev_inv is required and positive")
    fields = {
        "frequency_hz",
        "peak_rotation_rad_per_gev_inv_x", "peak_rotation_rad_per_gev_inv_y",
        "c_psd_hz_inv_per_unit_variance_x", "c_psd_hz_inv_per_unit_variance_y",
    }
    with np.load(src, allow_pickle=False) as z:
        if set(z.files) != fields:
            raise ACalibrationError("Calibration NPZ field names must match exact contract")
        f = np.asarray(z["frequency_hz"], dtype=np.float64)
        x, y = make_singlebin_peak_norm(
            f, reference_coupling_gev_inv=gref,
            **{key: z[key] for key in fields if key != "frequency_hz"}
        )
    metadata = {
        "mu_definition": "mu=(g_a_gamma_gamma/reference_coupling_gev_inv)^2",
        "reference_coupling_gev_inv": float(gref),
        "response_units": "rad^2/Hz per unit mu",
        "response_source": meta["response_source"],
        "calibration_version": meta["calibration_version"],
        "spectral_shape_convention": "original-vendor single-bin peak_shape=[1], PSD per unit time variance",
        "optical_peak_amplitude_source": meta["peak_amplitude_calibration_source"],
        "quadrature_projection_source": meta["quadrature_projection_source"],
        "native_c_kernel_source": meta["native_c_kernel_source"],
        "single_bin_convention": meta["single_bin_convention"],
        "response_status": meta["response_status"],
        "input_calibration_sha256": _sha256(src),
        "input_provenance_sha256": _sha256(info),
        "physical_upper_limits_authorized": False,
        "configuration": "A_original_vendor_single_bin",
    }
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        with dst.open("xb") as stream:
            np.savez(stream, frequency_hz=f, peak_norm_x=x, peak_norm_y=y,
                     metadata_json=json.dumps(metadata, sort_keys=True))
    except BaseException:
        dst.unlink(missing_ok=True)
        raise
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration-npz", required=True)
    parser.add_argument("--provenance-json", required=True)
    parser.add_argument("--output-npz", required=True)
    args = parser.parse_args()
    print(json.dumps(build_response_file(
        args.calibration_npz, args.provenance_json, args.output_npz
    ), indent=2))


if __name__ == "__main__":
    main()
