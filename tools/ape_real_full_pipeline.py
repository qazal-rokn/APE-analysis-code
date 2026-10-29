"""APE provisional real-data native PS -> original vendor fit -> expected limits.

Reuses unmodified upstream algorithms and never represents expected limits
as observed/published exclusion limits. Does not modify frozen inputs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import numpy as np


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, obj):
    path = Path(path)
    temp = path.with_name(path.name + ".partial")
    if temp.exists() or path.exists():
        raise FileExistsError(path)
    with temp.open("x") as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write("\n")
    os.replace(temp, path)


def load_native_ps(path):
    path = Path(path).expanduser().resolve(strict=True)
    a = np.load(path, allow_pickle=False, mmap_mode="r") if path.suffix == ".npy" else np.loadtxt(path, comments="#", ndmin=2)
    if a.ndim != 2 or a.shape[1] != 2 or len(a) < 1100:
        raise ValueError("Native C input must be >=1100 two-column frequency,PS rows: " + str(path))
    f, p = np.asarray(a[:, 0], float), np.asarray(a[:, 1], float)
    if not (np.isfinite(f).all() and np.isfinite(p).all() and (f > 0).all()
            and (p > 0).all() and (np.diff(f) > 0).all()):
        raise ValueError("Nonpositive or unordered native C power spectrum")
    return f, p, path


def run(cfg, cfgpath):
    if cfg.get("input_unit") != "native_PS_rad2":
        raise ValueError("Never treat PSD as native PS")
    if cfg.get("scientific_scope") != "real_data_provisional_diagnostic":
        raise ValueError("Not an authorized physical APE exclusion analysis")
    if cfg.get("peak_shape") != "vendor_single_bin":
        raise ValueError("Only the existing single-bin vendor bridge is available")

    from ape_c_lpsd.analysis.vendor_statistics import fit_window, load_external, FIT_SETTINGS
    from ape_c_lpsd.analysis.vendor_shape_ape import prepare_vendor_shape_ape_likelihood
    from ape_c_lpsd.analysis.ape_vendor_fit import fit_vendor_mle_and_expected_ul
    from ape_c_lpsd.analysis.ape_candidates import select_candidates

    fx, x, px = load_native_ps(cfg["native_ps"]["rho_x"])
    fy, y, py = load_native_ps(cfg["native_ps"]["rho_y"])
    if not np.array_equal(fx, fy):
        raise ValueError("Native C X/Y frequency grids differ")
    f = fx
    root = Path(cfg["external_vendor_root"]).expanduser().resolve(strict=True)
    _, _, vendor_hashes = load_external(root)
    out = Path(cfg["output_dir"]).expanduser().resolve()
    if out.exists() and not out.is_dir():
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    provenance = out / "provenance.json"
    source = {"input_x": str(px), "input_y": str(py), "sha_x": sha(px),
              "sha_y": sha(py), "config_sha": sha(cfgpath),
              "vendor_hashes": vendor_hashes, "n_bins": len(f)}
    if provenance.exists():
        if json.loads(provenance.read_text()) != source:
            raise RuntimeError("Output directory belongs to different source/config")
    else:
        atomic_json(provenance, source)

    options = cfg["background"]
    segment = int(options.get("segment_bins", 2000))
    if segment < 1100:
        raise ValueError("Background segments must contain >=1100 bins")
    buffer = int(options.get("buffer_bins", 50))
    if buffer < 0:
        raise ValueError("Negative fit buffer")
    fit_settings = dict(FIT_SETTINGS)
    fit_settings.update(options.get("vendor_fit_settings", {}))
    bkg = np.empty((2, len(f)), float)
    pars = np.empty((2, len(f), 3), float)
    quality = []
    cache = out / "background_parts"
    cache.mkdir(exist_ok=True)
    for chan, p in enumerate((x, y)):
        for k, start in enumerate(range(0, len(f), segment)):
            end = min(start + segment, len(f))
            lo = max(0, start - buffer)
            hi = min(len(f), end + buffer)
            if hi - lo < 1100:
                lo = max(0, hi - 1100)
            item = cache / f"{chan}_{k:06d}.npz"
            if item.exists():
                with np.load(item, allow_pickle=False) as d:
                    if int(d["start"]) != start or int(d["end"]) != end:
                        raise RuntimeError("Bad background cache bounds")
                    bb = np.array(d["background_ps"], float)
                    pp = np.array(d["model_args"], float)
                    chi = float(d["chi2"])
                    valid = bool(d["quality_pass"])
            else:
                fit = fit_window(f[lo:hi], p[lo:hi], root, settings=fit_settings)
                bb = np.array(fit["background_psd"][start - lo:end - lo])
                pp = np.array(fit["skewnormal_loc_scale_alpha"])
                chi = fit["vendor_reduced_chi2"]
                valid = fit["vendor_quality_pass"]
                if not valid:
                    raise RuntimeError(f"Vendor background quality failed channel={chan} block={k}, chi2={chi}")
                temp = item.with_suffix(".partial.npz")
                np.savez(temp, start=start, end=end, background_ps=bb,
                         model_args=pp, chi2=chi, quality_pass=valid)
                os.replace(temp, item)
            if not valid or len(bb) != end - start or pp.shape != (3,) or (bb <= 0).any() or not np.isfinite(bb).all():
                raise RuntimeError("Invalid fitted background")
            bkg[chan, start:end] = bb
            pars[chan, start:end, :] = pp
            quality.append([chan, k, chi])

    physics, opt = cfg["physics"], cfg["optimizer"]
    stride = int(cfg.get("scan_stride", 1))
    if stride < 1:
        raise ValueError("scan_stride must be positive")
    stop = int(cfg.get("max_scan_bins", len(f)))
    if stop < 1:
        raise ValueError("max_scan_bins must be positive")
    indices = np.arange(0, min(stop, len(f)), stride, dtype=int)
    results = out / "frequency_results.csv"
    rows = []
    if results.exists():
        with results.open() as reader:
            rows = list(csv.DictReader(reader))
        seen = [int(v["bin_index"]) for v in rows]
        if seen != indices[:len(seen)].tolist():
            raise RuntimeError("Resumed scan has wrong frequency indices")
    fields = ["bin_index", "frequency_hz", "mu_hat", "sigma_mu", "q0",
              "local_sqrt_q0", "expected_g_ul_gev_inv", "status"]
    for n, k in enumerate(indices[len(rows):], start=len(rows)):
        prepared = prepare_vendor_shape_ape_likelihood(
            observed_ps=np.array([[x[k]], [y[k]]]),
            background_ps=bkg[:, k:k+1], skewnormal_parameters=pars[:, k, :],
            frequency_bins_hz=np.array([f[k]]), vendor_peak_shape=np.array([1.]),
            quadrature_power_fraction=physics["quadrature_power_fraction"],
            reference_coupling_gev_inv=physics["reference_coupling_gev_inv"],
            finesse_s=physics["finesse_s"], tau_s=physics["tau_s"],
            tau1_s=physics["tau1_s"], rho_dm_gev4=physics["rho_dm_gev4"],
            unit_conversion=physics["unit_conversion"])
        fit = fit_vendor_mle_and_expected_ul(
            prepared=prepared, external_root=root,
            min_log10mu=opt["min_log10mu"], max_log10mu=opt["max_log10mu"],
            n_seed_lkl=opt["n_seed_lkl"], alpha_cl=opt.get("alpha_cl", 0.95))
        row = dict(bin_index=int(k), frequency_hz=float(f[k]),
                   mu_hat=fit["mu_hat"], sigma_mu=fit["sigma_mu"],
                   q0=fit["q0"], local_sqrt_q0=fit["local_sqrt_q0"],
                   expected_g_ul_gev_inv=fit["expected_upper_limit"]["g_ul_gev_inv"],
                   status=fit["status"])
        with results.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            if n == 0:
                writer.writeheader()
            writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())
        rows.append({key: str(val) for key, val in row.items()})
    if not rows:
        raise RuntimeError("No likelihood results")
    selected = select_candidates(
        np.array([float(r["frequency_hz"]) for r in rows]),
        np.zeros(len(rows)), np.array([float(r["q0"])/2 for r in rows]),
        np.array([float(r["mu_hat"]) for r in rows]),
        local_sigma_threshold=float(cfg["local_sigma_threshold"]))
    summary = {
        "status": "COMPLETED_PROVISIONAL_REAL_DATA_EXPECTED_LIMITS_NOT_PHYSICAL_EXCLUSION",
        "provenance": source, "total_native_bins": len(f),
        "analyzed_bins": len(rows), "scan_stride": stride,
        "full_grid": len(rows) == len(f) and stride == 1,
        "background_quality_rows": len(quality), "local_candidates": selected,
        "not_validated": ["physical_QWP_and_XY_response", "native_C_incoherent_averaging",
                          "cross_bin_XY_correlations", "global_significance",
                          "APE_limit_coverage", "negative_MLE_observed_limits"],
        "expected_limits_not_observed": True, "physical_g_limits_authorized": False}
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    frequencies = np.array([float(r["frequency_hz"]) for r in rows])
    limits = np.array([float(r["expected_g_ul_gev_inv"]) for r in rows])
    if not np.isfinite(limits).all() or not (limits > 0).all():
        raise RuntimeError("Invalid coupling values")
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.loglog(frequencies, limits, lw=0.8)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Provisional expected 95% g upper limit (GeV$^{-1}$)")
    ax.set_title("APE real-data diagnostic: NOT a validated exclusion")
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    image = out / "expected_coupling_diagnostic.png"
    if image.exists():
        raise FileExistsError(image)
    fig.savefig(image, dpi=180)
    plt.close(fig)
    atomic_json(out / "summary.json", summary)
    return {"status": summary["status"], "results": str(results),
            "summary": str(out / "summary.json"), "frequencies_scanned": len(rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    path = args.config.expanduser().resolve(strict=True)
    print(json.dumps(run(json.loads(path.read_text()), path), indent=2))


if __name__ == "__main__":
    main()
