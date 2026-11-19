"""Minimal APE-only bridge into the UNCHANGED external likelihood optimizer.

Uses APE signal-power response in place of LIGO get_beta().
The production Config A entry point analyzes only the summed Q PSD.
Single-bin spectral template only in this validation stage (peak_shape=[1]).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from scalar_dark_matter_lpsd.LPSDIO import LPSDOutput, LPSDJSONIO, MaxLKLIO
from scalar_dark_matter_lpsd.lvk.fit_background import bkg_model
from scalar_dark_matter_lpsd.lvk.max_lkl_executable import calc_max_lkl

from ape_post_lpsd.ape_signal_response import load_peak_norm, load_peak_norm_q


class APELikelihoodError(ValueError):
    """Fail closed on invalid APE-to-vendor likelihood inputs."""


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def run_single_quadrature(
    spectrum_path: str | Path,
    background_json: str | Path,
    response_npz: str | Path,
    output_npy: str | Path,
    *,
    channel: str,
    ana_fmin: float,
    ana_fmax: float,
    min_log10mu: int,
    max_log10mu: int,
    max_chi_sqr: float = 10.0,
    max_frequencies: int | None = None,
    start_index: int = 0,
    n_seed_lkl: int = 1000,
) -> dict:
    """Run ORIGINAL calc_max_lkl on one APE channel with supplied PSD response.

    No X/Y joint-likelihood independence assumption, detector calibration,
    spectral rebinning, new optimization, or physical coupling conversion.
    Output is original 5-column MaxLKLIO format: f, lnL0, lnLmax, mu, sigma.
    """
    if channel not in ("X", "Y", "Q"):
        raise APELikelihoodError("Specify one APE PSD: X, Y or phase-independent sum Q")
    spectrum = Path(spectrum_path).expanduser().resolve()
    background = Path(background_json).expanduser().resolve()
    response = Path(response_npz).expanduser().resolve()
    output = Path(output_npy).expanduser().resolve()
    report = output.with_suffix(".json")
    tmp = output.with_name(output.stem + ".partial.npy")
    if output.suffix != ".npy" or not spectrum.name.endswith(f"_{channel}.h5"):
        raise APELikelihoodError("Expected .npy output and APE_<t0>_<t1>_<channel>.h5 input")
    if not all(path.is_file() for path in (spectrum, background, response)):
        raise APELikelihoodError("Spectrum, fitted-background JSON, and APE response NPZ must exist")
    if any(path.exists() for path in (output, report, tmp)):
        raise APELikelihoodError("Refusing to overwrite an existing result or partial file")
    if not (0 < ana_fmin < ana_fmax and min_log10mu < max_log10mu
            and 0 < max_chi_sqr and n_seed_lkl >= 3):
        raise APELikelihoodError("Invalid analysis frequencies, mu scan bounds or fit parameters")
    if max_frequencies is not None and max_frequencies <= 0:
        raise APELikelihoodError("max_frequencies must be positive")
    if not isinstance(start_index, int) or start_index < 0:
        raise APELikelihoodError("start_index must be a nonnegative integer")

    data = LPSDOutput(str(spectrum))  # unchanged vendor reader
    if channel == "Q":
        channel_peak_norm, metadata = load_peak_norm_q(response, data.freq)
    else:
        peak_norm, metadata = load_peak_norm(response, data.freq)
        channel_peak_norm = peak_norm[0 if channel == "X" else 1]
    background_data = LPSDJSONIO(str(background))
    label = background_data.get_label(filepath=str(spectrum))
    table = background_data.get_df(label)
    required = {"fmin", "fmax", "x_knots", "y_knots", "alpha_skew",
                "loc_skew", "sigma_skew", "chi_sqr"}
    if not required.issubset(table.columns):
        raise APELikelihoodError("Vendor background record has missing columns")

    start = int(np.searchsorted(data.freq, ana_fmin, side="left")) + start_index
    stop = int(np.searchsorted(data.freq, ana_fmax, side="left"))
    if max_frequencies is not None:
        stop = min(stop, start + max_frequencies)
    if start >= stop:
        raise APELikelihoodError("No native C-LPSD frequencies in requested range")

    # APE-only replacement for LIGO stats.get_beta: the previously checked
    # response is in the vendor peak_norm array, with NO numerical adjustment.
    peak_shape = np.ones(1, dtype=np.float64)
    rows: list[tuple[float, float, float, float, float]] = []
    skipped = 0
    for j in range(start, stop):
        freq = float(data.freq[j])
        match = table[(table["fmin"] <= freq) & (table["fmax"] >= freq)]
        if len(match) == 0:
            skipped += 1
            continue
        fitted = match.iloc[0]  # original vendor convention
        if not (np.isfinite(fitted["chi_sqr"]) and fitted["chi_sqr"] < max_chi_sqr):
            skipped += 1
            continue
        knots = np.concatenate((fitted["x_knots"], fitted["y_knots"]))
        fitted_log_background = np.asarray(bkg_model(np.log([freq]), knots), dtype=np.float64)
        model_args = np.array([[fitted["loc_skew"], fitted["sigma_skew"], fitted["alpha_skew"]]], dtype=np.float64)
        Y = data.logPSD[j:j + 1][None, :]
        bkg = fitted_log_background[None, :]
        beta = channel_peak_norm[j:j + 1][None, :]
        if not (np.isfinite(Y).all() and np.isfinite(bkg).all()
                and np.isfinite(beta).all() and (beta >= 0).all()
                and np.isfinite(model_args).all() and model_args[0, 1] > 0):
            raise APELikelihoodError(f"Invalid original likelihood inputs at {freq} Hz")
        # Use the ORIGINAL optimizer, including signed-mu scan and Fisher code.
        # One channel only; do NOT pretend X and Y are independent LIGO IFOs.
        zero_lkl, max_lkl, mu, sigma = calc_max_lkl(
            Y, bkg, model_args, beta, peak_shape,
            min_log10mu=min_log10mu, max_log10mu=max_log10mu,
            n_seed_lkl=n_seed_lkl)
        rows.append((freq, zero_lkl, max_lkl, mu, sigma))
    if not rows:
        raise APELikelihoodError("No frequencies survived original background-quality criteria")

    output.parent.mkdir(parents=True, exist_ok=True)
    result = MaxLKLIO(data=np.asarray(rows, dtype=np.float64))
    provenance = {
        "status": "CONDITIONAL: APE response uses provisional working QWP parameters",
        "channel": channel,
        "source_spectrum": str(spectrum),
        "source_background": str(background),
        "source_response": str(response),
        "sha256": {"spectrum": _sha256(spectrum), "background": _sha256(background),
                   "response": _sha256(response)},
        "response_metadata": metadata,
        "vendor_functions": ["LPSDOutput", "LPSDJSONIO", "bkg_model", "calc_max_lkl", "MaxLKLIO"],
        "signal_model": "single-bin additive PSD: exp(bkg)+mu*peak_norm; peak_shape=1",
        "mu_is_signed_vendor_fit_parameter": True,
        "joint_xy_likelihood_used": False,
        "quadrature_psd_sum_used": channel == "Q",
        "frequency_range_hz": [ana_fmin, ana_fmax],
        "start_index_within_analysis_range": start_index,
        "native_grid_indices": [start, stop],
        "frequency_count": len(rows),
        "frequencies_skipped_by_background": skipped,
        "scan_log10mu": [min_log10mu, max_log10mu],
    }
    try:
        result.save(str(tmp))  # unchanged vendor MaxLKLIO output format
        tmp.rename(output)
        with report.open("x", encoding="utf-8") as stream:
            json.dump(provenance, stream, indent=2)
            stream.write("\n")
    except BaseException:
        tmp.unlink(missing_ok=True)
        output.unlink(missing_ok=True)
        raise
    return provenance


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spectrum", dest="spectrum_path", required=True)
    p.add_argument("--background-json", required=True)
    p.add_argument("--response-npz", required=True)
    p.add_argument("--output-npy", required=True)
    p.add_argument("--channel", required=True, choices=("X", "Y", "Q"))
    p.add_argument("--ana-fmin", required=True, type=float)
    p.add_argument("--ana-fmax", required=True, type=float)
    p.add_argument("--min-log10mu", required=True, type=int)
    p.add_argument("--max-log10mu", required=True, type=int)
    p.add_argument("--max-chi-sqr", type=float, default=10.0)
    p.add_argument("--max-frequencies", type=int, default=None)
    p.add_argument("--start-index", type=int, default=0,
                   help="Zero-based offset from first selected APE frequency for independent chunk jobs")
    args = vars(p.parse_args())
    print(json.dumps(run_single_quadrature(**{k.replace("-", "_"): v for k, v in args.items()}), indent=2))


if __name__ == "__main__":
    main()
