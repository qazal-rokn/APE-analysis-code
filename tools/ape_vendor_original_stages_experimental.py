"""Experimental APE adapter: ORIGINAL vendor stage wiring, unchanged source.

Run on a DIRECTORY containing ONLY combined APE LPSDOutput HDF5 files
named <prefix>_<t0>_<t1>_APE.h5. X/Y must be summed in LINEAR power
before these files are made. UNTESTED; no scientific-result authorization.

The only monkeypatch is scoped to vendor stats.get_beta, which represents
the detector-specific power normalization. Vendor make_args, process_q0,
background fitting, likelihood and optimizer are all called unchanged.
"""
from __future__ import annotations
import argparse
import contextlib
import json
from pathlib import Path

import numpy as np


@contextlib.contextmanager
def ape_beta(stats, response_frequency, response_beta):
    original = stats.get_beta
    frequency = np.asarray(response_frequency, float)
    beta = np.asarray(response_beta, float)
    if frequency.ndim != 1 or beta.shape != frequency.shape:
        raise ValueError("APE beta must be on the full vendor frequency grid")
    if not np.isfinite(beta).all() or np.any(beta < 0):
        raise ValueError("APE beta must be finite and nonnegative")
    def get_ape_beta(freq_Hz, unused_transfer):
        f = np.asarray(freq_Hz, float)
        positions = np.searchsorted(frequency, f)
        if np.any(positions >= len(frequency)) or not np.array_equal(frequency[positions], f):
            raise ValueError("APE response and vendor frequency grid differ")
        return beta[positions]
    try:
        stats.get_beta = get_ape_beta
        yield
    finally:
        stats.get_beta = original


def run(config):
    from ape_c_lpsd.analysis.vendor_statistics import load_external
    cfg = json.loads(Path(config).resolve(strict=True).read_text())
    root = Path(cfg["external_vendor_root"]).resolve(strict=True)
    folder = Path(cfg["combined_lpsd_directory"]).resolve(strict=True)
    prefix = cfg["data_prefix"]
    files = sorted(folder.glob(prefix + "*.h5"))
    if not files or any(not x.name.endswith("_APE.h5") for x in files):
        raise ValueError("Expected exclusively combined _t0_t1_APE.h5 vendor inputs")
    if sorted(x for x in folder.iterdir() if x.name.startswith(prefix)) != files:
        raise ValueError("Input directory contains unexpected prefixed files")
    out = Path(cfg["output_directory"]).resolve()
    if out.exists():
        raise FileExistsError("Output directory already exists")
    out.mkdir(parents=True)
    fit, stats, hashes = load_external(root)
    from scalar_dark_matter_lpsd.LPSDIO import LPSDDataGroup, LPSDJSONIO
    from scalar_dark_matter_lpsd.lvk import max_lkl_executable as mle
    from scalar_dark_matter_lpsd.lvk import calc_upper_limit as ul
    from scalar_dark_matter_lpsd.lvk import get_candidates as candidates

    # Original vendor background workflow and original JSON structure.
    background = out / "vendor_background.json"
    for path in files:
        fit.fit_background_in_file(
            str(path), str(background), **cfg.get("vendor_background_options", {}))
    group = LPSDDataGroup(str(folder), prefix, buffer_path=str(out / "vendor_group_buffer.h5"))
    if group.logPSD.shape[0] != len(files):
        raise RuntimeError("Original vendor loader did not load all combined records")
    if set(group.metadata["ifo"]) != {"APE"}:
        raise RuntimeError("Expected one combined APE spectrum per observation")
    bkg = LPSDJSONIO(str(background))
    peak = np.load(Path(cfg["vendor_peak_shape_npy"]).resolve(strict=True), allow_pickle=False)
    if peak.ndim != 1 or not peak.size or not np.isfinite(peak).all() or np.any(peak < 0) or peak.sum() != 1.0:
        raise ValueError("Require unchanged vendor normalized peak shape")
    with np.load(Path(cfg["ape_beta_grid_npz"]).resolve(strict=True), allow_pickle=False) as d:
        response_f = np.array(d["frequency_hz"], float)
        response_beta = np.array(d["beta_ape"], float)
    if not np.array_equal(response_f, group.freq):
        raise ValueError("APE beta grid must exactly match vendor frequency grid")

    # The original vendor make_args creates every likelihood input and applies
    # original NaN, background and chi-square filtering; process_q0 invokes
    # the original signed-mu MLE, Fisher sigma, and likelihood.
    with ape_beta(stats, response_f, response_beta):
        args = mle.make_args(
            group, bkg, {"APE": None}, iteration=0, n_frequencies=-1,
            ana_fmin=float(cfg.get("ana_fmin", 10)),
            ana_fmax=float(cfg.get("ana_fmax", 5000)),
            peak_shape=peak,
            max_chi_sqr=float(cfg.get("max_chi_sqr", 10)),
            min_log10mu=float(cfg.get("min_log10mu", -40)),
            max_log10mu=float(cfg.get("max_log10mu", -32)))
        results = np.asarray([mle.process_q0(item) for item in args], dtype=float)
    if results.ndim != 2 or results.shape[1] != 5 or not len(results):
        raise RuntimeError("No valid vendor likelihood results")
    from scalar_dark_matter_lpsd.LPSDIO import MaxLKLIO
    MaxLKLIO(data=results).save(str(out / "vendor_mle_result.npy"))
    # Original vendor ul.main contains hard-coded LIGO filesystem paths,
    # transfer-function loading, and plots. Its original numerical upper-limit
    # function is portable and reused directly; the main() driver is not.
    limits = np.full((len(results), 3), np.nan)
    limits[:, 0] = results[:, 0]
    for i, sigma in enumerate(results[:, 4]):
        if np.isfinite(sigma) and sigma > 0:
            limit, uncertainty = ul.get_upper_limits_approx(
                float(sigma), alpha_CL=float(cfg.get("alpha_CL", 0.95)))
            limits[i, 1] = float(np.asarray(limit).reshape(-1)[0])
            limits[i, 2] = float(np.asarray(uncertainty).reshape(-1)[0])
    np.savetxt(out / "vendor_expected_sqrt_mu.csv", limits, delimiter=",",
               header="frequency_hz,sqrt_mu_ul,sigma_sqrt_mu", comments="")
    # The original get_candidates.main includes a two-IFO consistency test
    # that is not applicable to the single combined APE channel. Do not
    # silently run or reinterpret it.
    (out / "STATUS.json").write_text(json.dumps({
        "status": "EXPERIMENTAL_UNTESTED_NOT_VALIDATED",
        "vendor_background_source_sha256": hashes,
        "combined_observation_count": len(files),
        "likelihood_bins": len(results),
        "reused": ["fit_background_in_file", "make_args", "process_q0",
                   "get_upper_limits_approx"],
        "only_detector_response_override": "stats.get_beta scoped to make_args",
        "not_implemented": ["vendor two-observatory consistency test",
                            "verified vendor candidate workflow for one APE channel",
                            "automatic native C-LPSD to combined vendor HDF5 wiring",
                            "vendor upper-limit plotting and LIGO-only conversions"],
        "no_validation_or_parity_tests_run": True
    }, indent=2))
    return str(out / "STATUS.json")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    print(run(p.parse_args().config))
