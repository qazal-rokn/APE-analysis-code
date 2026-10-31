"""Single post-LPSD APE entry point. Existing APE ingestion/calibration/C-LPSD unchanged.

Uses original vendor LPSDOutput, LPSDData, fit_background_in_file,
LPSDDataGroup, LPSDJSONIO, make_args, process_q0, MaxLKLIO,
get_upper_limits_approx and the vendor candidate prefix adapter.
Only new numerical operation: PS_X + PS_Y.
Experimental: do not use for scientific results until validated.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import h5py


def run(config):
    from scalar_dark_matter_lpsd.LPSDIO import (
        LPSDOutput, LPSDData, LPSDDataGroup, LPSDJSONIO, MaxLKLIO)
    from scalar_dark_matter_lpsd.lvk import (
        fit_background, max_lkl_executable, calc_upper_limit, stats)
    from ape_vendor_original_stages_experimental import ape_beta
    from ape_vendor_candidate_original_prefix import run_vendor_candidate_prefix
    from ape_vendor_upper_limit_original_prefix import run_vendor_upper_limit_prefix

    cfg = json.loads(Path(config).resolve(strict=True).read_text())
    out = Path(cfg["output_directory"]).resolve()
    if out.exists():
        raise FileExistsError(out)
    observations = cfg["observations"]
    if not observations:
        raise ValueError("No X/Y C-LPSD pairs supplied")
    with np.load(cfg["ape_beta_grid_npz"], allow_pickle=False) as response:
        beta_freq = np.array(response["frequency_hz"])
        beta = np.array(response["beta_ape"])
    peak = (np.load(cfg["vendor_peak_shape_npy"], allow_pickle=False)
            if cfg.get("vendor_peak_shape_npy") else np.ones(1))
    if peak.ndim != 1 or np.sum(peak) != 1.:
        raise ValueError("Original vendor peak-shape normalization failed")

    out.mkdir(parents=True)
    combined = out / "combined"
    combined.mkdir()
    background_path = out / "vendor_background.json"
    seen = set()
    background_options = {
        "n_processes": 1, "buffer": 50, "pruning": 1,
        "segment_size": 1000, "ana_fmin": 10., "ana_fmax": 5000.,
        "k_min": 4, "k_pruning": 1, "k_plateau": 3,
        "bin_scaling": .15, "verbose": False, "plot_path": None,
        "kernel_size": 30, "cfd_fraction": .1
    }
    background_options.update(cfg.get("vendor_background_options", {}))
    for obs in observations:
        x = LPSDOutput(str(Path(obs["rho_x_lpsd"]).resolve(strict=True)))
        y = LPSDOutput(str(Path(obs["rho_y_lpsd"]).resolve(strict=True)))
        if not np.array_equal(x.freq, y.freq, equal_nan=True):
            raise ValueError("Native X/Y grids differ")
        if set(x.kwargs) != set(y.kwargs):
            raise ValueError("Native X/Y metadata keys differ")
        for key in x.kwargs:
            if x.kwargs[key] is None or y.kwargs[key] is None:
                if x.kwargs[key] is not None or y.kwargs[key] is not None:
                    raise ValueError("Native X/Y metadata differ")
            elif not np.array_equal(np.asarray(x.kwargs[key]),
                                    np.asarray(y.kwargs[key]), equal_nan=True):
                raise ValueError("Native X/Y metadata differ")
        if x.logPSD.shape != y.logPSD.shape:
            raise ValueError("Native X/Y PS shapes differ")
        t0, t1 = int(obs["t0"]), int(obs["t1"])
        if t1 <= t0:
            raise ValueError("Invalid observation interval")
        name = f"ape_{t0}_{t1}_APE.h5"
        if name in seen:
            raise ValueError("Duplicate observation")
        seen.add(name)
        # Stable linear-power summation in log domain; preserves NaN gaps.
        summed = np.logaddexp(x.logPSD, y.logPSD)
        path = combined / name
        # Vendor LPSDData.__init__ does not set self.freq when freq is
        # provided. Preserve native missing-bin positions with the exact
        # original HDF5 datasets/attributes instead of calling that writer.
        with h5py.File(path, "w") as file:
            file.create_dataset("frequency", data=x.freq)
            file.create_dataset("logPSD", data=summed)
            file["logPSD"].attrs.update({k: (-1 if v is None else v)
                                         for k, v in x.kwargs.items()})
        check = LPSDOutput(str(path))
        if not np.array_equal(check.freq, x.freq, equal_nan=True):
            raise RuntimeError("Combined file changed native frequency grid")
        if not np.array_equal(check.logPSD, summed, equal_nan=True):
            raise RuntimeError("Combined file changed summed power")
        fit_background.fit_background_in_file(
            str(path), str(background_path),
            **background_options)

    group = LPSDDataGroup(str(combined), "ape_",
                          buffer_path=str(out / "vendor_group_buffer.h5"))
    if group.logPSD.shape[0] != len(seen):
        raise RuntimeError("Original vendor group failed to load all observations")
    if set(group.metadata["ifo"]) != {"APE"}:
        raise RuntimeError("Unexpected detector metadata")
    if not np.array_equal(group.freq, beta_freq, equal_nan=True):
        raise ValueError("APE response must use exact native C-LPSD grid")
    bkg = LPSDJSONIO(str(background_path))
    if not np.all(np.diff(beta_freq) > 0):
        raise ValueError("APE beta grid must be strictly increasing")
    with ape_beta(stats, beta_freq, beta):
        args = max_lkl_executable.make_args(
            group, bkg, {"APE": None}, 0, -1,
            float(cfg.get("ana_fmin", 10)), float(cfg.get("ana_fmax", 5000)),
            peak, float(cfg.get("max_chi_sqr", 10)),
            float(cfg.get("min_log10mu", -40)),
            float(cfg.get("max_log10mu", -32)))
        rows = np.asarray(list(map(max_lkl_executable.process_q0, args)))
    if rows.ndim != 2 or rows.shape[1] != 5 or len(rows) == 0:
        raise RuntimeError("Original vendor likelihood produced no valid rows")
    mle_file = out / "vendor_max_lkl.npy"
    MaxLKLIO(data=rows).save(str(mle_file))
    candidate_result = run_vendor_candidate_prefix(
        mle_file, out / "candidates.csv", float(cfg.get("sigma_limit", 5)))

    # Reuse the original vendor's entire 1000-point UL lookup-table
    # algorithm and interpolation. Omit only LIGO-specific plotting,
    # comparisons, and strain-equivalent conversion.
    expected = run_vendor_upper_limit_prefix(
        mle_file, out / "vendor_expected_limits.tsv",
        pruning=int(cfg.get("ul_pruning", 1)),
        alpha_CL=float(cfg.get("alpha_CL", .95)),
        tol=float(cfg.get("ul_tol", 1e-10)),
        kernel_size=int(cfg.get("ul_kernel_size", 100000)))
    status = {"status": "EXPERIMENTAL_UNTESTED", "scope": "post-C-LPSD",
              "original_vendor_stages": [
                  "LPSDOutput", "LPSDData.save_to_hdf5",
                  "fit_background_in_file", "LPSDDataGroup", "LPSDJSONIO",
                  "make_args", "process_q0", "MaxLKLIO",
                  "get_candidates.main original-statistics prefix",
                  "original upper-limit lookup-table and interpolation"],
              "necessary_ape_changes": ["X+Y power sum", "APE beta response",
                  "exclude LIGO-only consistency and plots"],
              "observations": len(seen), "likelihood_rows": len(rows),
              "candidate_result": candidate_result}
    (out / "STATUS.json").write_text(json.dumps(status, indent=2))
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    print(json.dumps(run(parser.parse_args().config), indent=2))
