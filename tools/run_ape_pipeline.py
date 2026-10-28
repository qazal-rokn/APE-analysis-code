"""Integrated, explicit APE native-C PS -> vendor -> candidates controller.

Two modes: 'synthetic_halo_file' and 'calibrated_file', both are READ ONLY.
Real observations cannot be analyzed without separately approved real-data
science validation; --permit-real-data flag is intentionally absent.
C uses exact external binary and unchanged external vendor functions.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from ape_c_lpsd.spectral.configured_ps import plan_configured_ps
from ape_c_lpsd.backends.local import run_local
from ape_c_lpsd.analysis.ape_pipeline_controller import finalize_ps_plan
from ape_c_lpsd.analysis.vendor_statistics import fit_window
from ape_c_lpsd.analysis.vendor_shape_ape import load_vendor_peak_shape
from ape_c_lpsd.analysis.ape_candidate_scan import scan_candidate_windows
from ape_c_lpsd.analysis.vendor_shape_ape import prepare_vendor_shape_ape_likelihood
from ape_c_lpsd.analysis.ape_vendor_candidate_report import fit_and_report_single_bin_candidates
from tools.run_ape_synthetic_candidates import join_ps_sidecars
from ape_c_lpsd.analysis.multibin_diagnostic import evaluate_multibin_windows


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def validate_input(config):
    source=Path(config["input_hdf5"]).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if config["input_kind"]!="synthetic_halo_file":
        raise ValueError("Real-data processing requires separate scientific approval; refused")
    with h5py.File(source,"r") as h5:
        if not bool(h5.attrs.get("synthetic_only",False)):
            raise ValueError("Missing synthetic_only=true provenance")
        if h5.attrs.get("unit")!="rad":
            raise ValueError("Require calibrated X/Y in radians")
        fs=float(h5.attrs["fs_hz"])
        expected=int(round(fs*float(config["c_settings"]["duration_s"])))
        if fs!=float(config["c_settings"]["fs_hz"]):
            raise ValueError("HDF5 sampling rate does not match C settings")
        for key in ("rho_x","rho_y"):
            if key not in h5 or h5[key].ndim!=1 or len(h5[key])!=expected:
                raise ValueError(f"Invalid {key} duration/sample count")
        if config.get("require_stochastic_halo",True):
            if h5.attrs.get("model")!="stochastic_gaussian_vendor_shape":
                raise ValueError("This run requires a stochastic halo injection")
        if "signal_frequency_hz" in h5 and "vendor_shape_weights" in h5:
            shape=np.asarray(h5["vendor_shape_weights"][:])
            if not np.isclose(shape.sum(),1.,atol=1e-12,rtol=1e-12):
                raise ValueError("Injected vendor weights must sum to one")
    return source


def run_vendor_mle_stage(frequency, observed, fits, starts, config, vendor):
    """One-bin original vendor optimizer / local candidates / expected UL.

    Inputs are previously computed native-C PS, not Python spectral estimates.
    """
    if config.get("peak_shape_file") is not None:
        raise ValueError("Original vendor MLE pilot supports only one-bin [1]")
    if config.get("n_effective_trials") is not None:
        raise ValueError("No global trials calibration in this pilot")
    options = config["vendor_optimizer"]
    if set(options) != {"min_log10mu", "max_log10mu", "n_seed_lkl"}:
        raise ValueError("Explicit original vendor optimizer options required")
    indices = np.asarray(starts, dtype=int)
    if (indices.ndim != 1 or not indices.size or indices[0] < 0
            or indices[-1] >= len(frequency)
            or np.any(np.diff(indices) <= 0)):
        raise ValueError("Invalid strictly increasing C-grid indices")
    physics = config["physics"]
    prepared = []
    for k in indices:
        prepared.append(prepare_vendor_shape_ape_likelihood(
            observed_ps=observed[:, k:k+1],
            background_ps=np.stack([r["background_psd"][k:k+1] for r in fits]),
            skewnormal_parameters=np.stack(
                [r["skewnormal_loc_scale_alpha"] for r in fits]),
            frequency_bins_hz=frequency[k:k+1],
            vendor_peak_shape=np.array([1.]),
            quadrature_power_fraction=physics["quadrature_power_fraction"],
            reference_coupling_gev_inv=physics["reference_coupling_gev_inv"],
            finesse_s=physics["finesse_s"], tau_s=physics["tau_s"],
            tau1_s=physics["tau1_s"], rho_dm_gev4=physics["rho_dm_gev4"],
            unit_conversion=physics["unit_conversion"]))
    return fit_and_report_single_bin_candidates(
        candidate_frequencies_hz=frequency[indices],
        prepared_windows=prepared, external_root=vendor,
        min_log10mu=options["min_log10mu"],
        max_log10mu=options["max_log10mu"],
        n_seed_lkl=options["n_seed_lkl"],
        local_sigma_threshold=config["local_sigma_threshold"],
        alpha_cl=config.get("alpha_cl", 0.95))


def execute(config):
    """Explicit synthetic-run execution, no self-generating signal or silent defaults."""
    if config.get("synthetic_only") is not True or config.get("execute_c") is not True:
        raise ValueError("Requires explicit synthetic_only and execute_c true")
    # Reject unsafe or pre-existing output paths before inspecting input files.
    # This is a no-overwrite guard, not a change to scientific processing.
    output=Path(config["run_dir"]).expanduser().resolve()
    if output.exists() or not output.name.startswith("ape-pipeline-"):
        raise ValueError("Use a brand-new ape-pipeline-* directory")
    source=validate_input(config)
    if output==source or output in source.parents:
        raise ValueError("Unsafe output directory")
    # Verify executable *before* opening new outputs.
    exe=Path(config["ps_executable"]).expanduser().resolve(strict=True)
    binary_hash=sha256(exe)
    if binary_hash.lower()!=config["ps_sha256"].lower():
        raise ValueError("PS-enabled executable SHA256 mismatch")
    vendor=Path(config["external_vendor_root"]).expanduser().resolve(strict=True)
    if not (vendor/"scalar_dark_matter_lpsd/lvk/stats.py").is_file():
        raise FileNotFoundError("Original vendor stats.py unavailable")
    output.mkdir(parents=True,exist_ok=False)
    paths={}; provenance={}
    for channel in ("rho_x","rho_y"):
        stage=output/f"ape-ps-{channel}"
        plan=plan_configured_ps(
            executable=exe,ps_executable=exe,ps_sha256=binary_hash,
            input_hdf5=source,dataset=channel,fft_hdf5=stage/"fft.h5",
            output_dir=stage,settings=config["c_settings"])
        command_results=run_local(plan,synthetic_approved=True)
        batch=finalize_ps_plan(plan)
        joined=stage/f"{channel}.joined.ps"
        rows=join_ps_sidecars(batch,joined)
        if rows!=int(config["c_settings"]["grid_bins"]):
            raise ValueError("Native C bin count mismatch")
        paths[channel]=str(joined)
        provenance[channel]=dict(commands=command_results,
                                 batch_count=len(batch),frequency_bins=rows)
    x=np.loadtxt(paths["rho_x"],ndmin=2)
    y=np.loadtxt(paths["rho_y"],ndmin=2)
    if x.shape!=y.shape or x.shape[1]!=2 or not np.array_equal(x[:,0],y[:,0]):
        raise ValueError("Native X/Y frequency grids do not match exactly")
    frequency=x[:,0]
    observed=np.stack((x[:,1],y[:,1]))
    if np.any(observed<=0):
        raise ValueError("PS must be strictly positive for log-vendor likelihood")
    fits=[fit_window(frequency,p,vendor) for p in observed]
    if any(not fit["vendor_quality_pass"] for fit in fits):
        raise RuntimeError("External background fitting did not pass quality gate")
    # Original vendor default is precisely one bin when no shape file supplied.
    if config.get("peak_shape_file"):
        peak_shape=load_vendor_peak_shape(config["peak_shape_file"],
                                           frequency_bins_hz=frequency)
        shape_provenance=dict(path=str(config["peak_shape_file"]),
                              sha256=sha256(config["peak_shape_file"]))
    else:
        peak_shape=np.ones(1)
        shape_provenance=dict(vendor_original_default_single_bin=True)
    starts=np.asarray(config["candidate_start_indices"],dtype=int)
    if (starts.ndim!=1 or not starts.size or np.any(starts<0) or
            np.any(np.diff(starts)<=0) or np.any(starts+len(peak_shape)>len(frequency))):
        raise ValueError("Candidate window start indices invalid")
    fwindow=np.stack([frequency[k:k+len(peak_shape)] for k in starts])
    pwindow=np.stack([observed[:,k:k+len(peak_shape)] for k in starts])
    bgwindow=np.stack([np.stack([fit["background_psd"][k:k+len(peak_shape)]
                                  for fit in fits]) for k in starts])
    args=np.stack([fit["skewnormal_loc_scale_alpha"] for fit in fits])
    skew=np.broadcast_to(args,(len(starts),*args.shape)).copy()
    physics=config["physics"]
    mode=config.get("statistical_mode","finite_grid")
    if mode=="vendor_mle_single_bin":
        scan=run_vendor_mle_stage(frequency, observed, fits, starts, config, vendor)
        candidates=scan["selection"]
    elif mode=="physical_multibin_template":
        scan=evaluate_multibin_windows(
            frequency_hz=frequency, observed_ps=observed,
            fitted_backgrounds=fits, candidate_start_indices=starts,
            config=config, external_vendor_root=vendor)
        candidates={"candidates": [], "status": "multibin_diagnostic_no_candidate_inference"}
    elif mode=="finite_grid":
        scan=scan_candidate_windows(
            candidate_frequencies_hz=frequency[starts],
            observed_ps=pwindow,background_ps=bgwindow,
            skewnormal_parameters=skew,frequency_windows_hz=fwindow,
            vendor_peak_shape=peak_shape,
            quadrature_power_fraction=physics["quadrature_power_fraction"],
            reference_coupling_gev_inv=physics["reference_coupling_gev_inv"],
            coupling_grid_gev_inv=config["coupling_grid_gev_inv"],
            external_root=vendor,
            local_sigma_threshold=config["local_sigma_threshold"],
            finesse_s=physics["finesse_s"],tau_s=physics["tau_s"],
            tau1_s=physics["tau1_s"],rho_dm_gev4=physics["rho_dm_gev4"],
            unit_conversion=physics["unit_conversion"],
            n_effective_trials=config.get("n_effective_trials"))
        candidates=scan
    else:
        raise ValueError("Unknown statistical_mode")
    manifest=dict(status="provisional_software_pipeline_not_validated_inference",
                  synthetic_only=True,source_sha256=sha256(source),
                  native_ps_executable_sha256=binary_hash,
                  c_settings=config["c_settings"],c_executions=provenance,
                  outputs=paths,fit_quality=[f["vendor_reduced_chi2"] for f in fits],
                  vendor_sha256=fits[0]["vendor_source_sha256"],
                  shape_provenance=shape_provenance,
                  scientific_warning="Unverified native window response, X/Y correlations, injection recovery and coverage",
                  statistical_mode=mode,
                  candidates=candidates,
                  vendor_mle_report=(scan if mode=="vendor_mle_single_bin" else None),
                  multibin_diagnostic_report=(scan if mode=="physical_multibin_template" else None))
    with (output/"manifest.json").open("x") as handle:
        json.dump(manifest,handle,indent=2)
    with (output/"candidates.csv").open("x",newline="") as handle:
        keys=("frequency_hz","q0","local_sqrt_q0","mu_hat",
              "cluster_start_index","cluster_end_index_exclusive","status")
        writer=csv.DictWriter(handle,fieldnames=keys)
        writer.writeheader(); writer.writerows(candidates["candidates"])
    return manifest


def main():
    parser=argparse.ArgumentParser(description="Approved synthetic-only native C PS -> APE likelihood -> candidate report")
    parser.add_argument("--config",required=True)
    args=parser.parse_args()
    config=json.loads(Path(args.config).read_text())
    result=execute(config)
    print(json.dumps(dict(status=result["status"],
                          candidate_count=len(result["candidates"]["candidates"])),
                     indent=2))


if __name__=="__main__":
    main()
