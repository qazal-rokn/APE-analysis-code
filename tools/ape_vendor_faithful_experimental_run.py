"""EXPERIMENTAL vendor-faithful post-C-LPSD runner. NOT VALIDATED.

Input: original-vendor-readable LPSDOutput files, each with genuine observation
metadata and a matching APE response. No changes to original vendor code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ape_c_lpsd.analysis.vendor_statistics import load_external


def run(config_file):
    config_file = Path(config_file).resolve(strict=True)
    cfg = json.loads(config_file.read_text())
    root = Path(cfg["external_vendor_root"]).resolve(strict=True)
    out = Path(cfg["output_dir"]).resolve()
    if out.exists():
        raise FileExistsError("Refusing to overwrite or resume experimental results")
    out.mkdir(parents=True)
    from scalar_dark_matter_lpsd.LPSDIO import LPSDDataGroup, LPSDJSONIO
    from scalar_dark_matter_lpsd.lvk.fit_background import fit_background_in_file, bkg_model
    from scalar_dark_matter_lpsd.lvk.max_lkl_executable import process_q0
    from scalar_dark_matter_lpsd.lvk.calc_upper_limit import get_upper_limits_approx
    from scalar_dark_matter_lpsd.lvk.get_candidates import cluster
    from scipy.stats import norm

    load_external(root)
    paths = [Path(x).resolve(strict=True) for x in cfg["vendor_lpsd_files"]]
    bkg_path = out / "vendor_background.json"
    bkg_opts = dict(cfg.get("vendor_background_options", {}))
    for path in paths:
        fit_background_in_file(str(path), str(bkg_path), **bkg_opts)
    data = LPSDDataGroup(str(paths[0]) if len(paths) == 1 else str(paths[0].parent),
                         cfg.get("data_prefix", ""))
    # A multi-file directory MUST contain exactly the intended inputs.
    if len(paths) > 1:
        discovered = sorted(x.resolve() for x in paths[0].parent.glob(cfg.get("data_prefix", "") + "*")
                            if x.is_file())
        if discovered != sorted(paths):
            raise ValueError("Vendor data directory includes files outside explicit manifest")
    bkg_data = LPSDJSONIO(str(bkg_path))
    bkg_info = {key: bkg_data.get_df(key) for key in bkg_data.data}
    shape = np.load(Path(cfg["vendor_peak_shape_npy"]).resolve(strict=True), allow_pickle=False)
    if shape.ndim != 1 or not shape.size or np.any(shape < 0) or not np.isfinite(shape).all() or shape.sum() != 1.:
        raise ValueError("Require original vendor unit-sum peak shape")
    # The APE response is a per-record, per-frequency reference signal POWER,
    # stored on the original C-LPSD grid. Its origin is deliberately explicit.
    response = np.load(Path(cfg["ape_reference_signal_power_npy"]).resolve(strict=True),
                       allow_pickle=False)
    if response.shape != data.logPSD.shape or not np.isfinite(response).all() or np.any(response < 0):
        raise ValueError("APE response must match vendor data (records, frequencies)")
    f = data.freq
    start = np.where(f >= cfg.get("ana_fmin", 10.))[0]
    if not len(start):
        raise ValueError("No frequency above ana_fmin")
    start = int(start[0])
    stop = np.where(f >= cfg.get("ana_fmax", 5000.))[0]
    stop = int(stop[0]) if len(stop) else len(f)
    kwargs = {"min_log10mu": cfg.get("min_log10mu", -40),
              "max_log10mu": cfg.get("max_log10mu", -32)}
    max_chi = cfg.get("max_chi_sqr", 10.)
    results = []
    for idx in range(start, stop):
        if idx + len(shape) > len(f):
            continue
        freq = f[idx:idx + len(shape)]
        Y, B, P, N = [], [], [], []
        for j, (t0, t1, ifo) in enumerate(zip(data.metadata["t0"],
                                               data.metadata["t1"],
                                               data.metadata["ifo"])):
            observed = data.logPSD[j, idx:idx + len(shape)]
            if any(np.isnan(observed)):
                continue
            label = bkg_data.get_label(t0=t0, t1=t1, ifo=ifo)
            df = bkg_info[label]
            rows = df[(df["fmin"] <= f[idx]) & (df["fmax"] >= f[idx])]
            if len(rows) == 0:
                continue
            xk, yk, alpha, loc, scale, chi = rows[
                ["x_knots", "y_knots", "alpha_skew", "loc_skew",
                 "sigma_skew", "chi_sqr"]].iloc[0]
            if chi >= max_chi:
                continue
            Y.append(observed)
            B.append(bkg_model(np.log(freq), np.concatenate([xk, yk])))
            P.append([loc, scale, alpha])
            # This is the ONLY physical replacement for vendor get_beta().
            # peak_norm * shape must equal the APE reference signal power.
            # Zero-shape bins must have zero physical signal power.
            power = response[j, idx:idx + len(shape)]
            if np.any((shape == 0) & (power != 0)):
                raise ValueError("APE response nonzero where vendor shape is zero")
            N.append(np.divide(power, shape, out=np.zeros_like(power),
                               where=shape > 0))
        if not Y:
            continue
        args = (np.asarray(Y), np.asarray(B), np.asarray(P), np.asarray(N),
                shape, f[idx], kwargs)
        results.append(process_q0(args))
    if not results:
        raise RuntimeError("No valid vendor likelihood windows")
    results = np.asarray(results, dtype=float)
    np.save(out / "vendor_format_max_lkl.npy", results)
    # Vendor get_candidates.main cannot be reused unchanged for a single APE
    # apparatus: its H1/L1 consistency check is detector-specific and contains
    # an interactive pdb breakpoint. Preserve vendor q0, trial factor and
    # cluster functions, and explicitly omit only that detector-specific veto.
    q0 = np.zeros(len(results))
    positive = results[:, 3] > 0
    q0[positive] = -2 * (results[positive, 1] - results[positive, 2])
    good = ~(np.isnan(q0) | (q0 >= 500))
    cleaned = q0[good]
    threshold = -norm.ppf(norm.cdf(-cfg.get("sigma_limit", 5.)) / len(cleaned))
    indices, amplitudes = cluster(cleaned > threshold ** 2, results[good, 3])
    candidate_freq = results[good, 0][indices]
    np.savetxt(out / "candidate_frequencies.csv", candidate_freq, delimiter=",")
    # The original vendor expected-limit function is used directly; this is
    # not an observed upper limit or a physical-exclusion validation.
    limits = np.full((len(results), 3), np.nan)
    limits[:, 0] = results[:, 0]
    for i, sigma in enumerate(results[:, 4]):
        if np.isfinite(sigma) and sigma > 0:
            lim, err = get_upper_limits_approx(
                sigma, alpha_CL=cfg.get("alpha_CL", .95))
            limits[i, 1:] = float(np.asarray(lim).reshape(-1)[0]), float(np.asarray(err).reshape(-1)[0])
    np.savetxt(out / "vendor_expected_sqrt_mu_limits.csv", limits, delimiter=",",
               header="frequency_hz,sqrt_mu_ul,sigma_sqrt_mu", comments="")
    manifest = {
        "status": "EXPERIMENTAL_UNTESTED_NOT_VALIDATED",
        "n_likelihood_frequencies": len(results),
        "n_candidates_without_detector_veto": len(candidate_freq),
        "vendor_source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (root / "scalar_dark_matter_lpsd/lvk/max_lkl_executable.py",
                      root / "scalar_dark_matter_lpsd/lvk/get_candidates.py",
                      root / "scalar_dark_matter_lpsd/lvk/fit_background.py")},
        "differences": [
            "APE reference power replaces vendor get_beta response",
            "Original H1/L1 detector consistency veto omitted for APE",
            "Input conversion to vendor LPSDOutput format is a separate prerequisite",
            "No parity or coverage tests have been executed",
        ],
    }
    (out / "status.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.config), indent=2))


if __name__ == "__main__":
    main()
