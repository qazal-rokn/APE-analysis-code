"""Explicit, local, synthetic-only C->PS->vendor fit->APE candidate smoke runner.

One command can execute only when a verified native PS C executable and the
unaltered external vendor package are present. Never runs real data.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import h5py
from ape_c_lpsd.analysis.synthetic_smoke import make_fixture
from ape_c_lpsd.spectral.short_synthetic_ps import plan_short_synthetic
from ape_c_lpsd.backends.local import run_local
from ape_c_lpsd.analysis.ape_pipeline_controller import finalize_ps_plan
from ape_c_lpsd.analysis.vendor_statistics import fit_window
from ape_c_lpsd.analysis.ape_candidate_scan import scan_candidate_windows
from ape_c_lpsd.analysis.vendor_shape_ape import load_vendor_peak_shape


def join_ps_sidecars(batch_results, output_path):
    """Join ordered native two-column PS sidecars without numerical rescaling."""
    path=Path(output_path)
    if path.exists():
        raise FileExistsError(path)
    files=sorted(batch_results,key=lambda a:a["start_index"])
    if not files or files[0]["start_index"]!=0:
        raise ValueError("First frequency partition is missing")
    next_index=0
    previous=0.
    n=0
    with path.open("x") as output:
        for info in files:
            if info["start_index"]!=next_index:
                raise ValueError("Missing/overlapping frequency partition")
            count=0
            with open(info["ps_file"]) as source:
                for line in source:
                    fields=line.split()
                    if len(fields)!=2:
                        raise ValueError("Expected native f PS two-column sidecar")
                    f,p=(float(v) for v in fields)
                    if not np.isfinite(f+p) or f<=previous or p<0:
                        raise ValueError("Invalid native PS or nonmonotonic frequency")
                    output.write(line)
                    previous=f
                    count+=1
            if count!=info["rows"]:
                raise ValueError("Incorrect rows in PS sidecar")
            next_index+=count
            n+=count
    return n


def run(config):
    """Require caller's explicit synthetic approval and never overwrite run directory."""
    if config.get("synthetic_only") is not True or config.get("execute_c") is not True:
        raise ValueError("Synthetic-only execution requires synthetic_only=true and execute_c=true")
    root=Path(config["run_dir"]).expanduser().resolve()
    if root.exists() or not root.name.startswith("ape-synthetic-ps-"):
        raise ValueError("Choose a new ape-synthetic-ps-* run directory")
    root.mkdir(parents=True,exist_ok=False)
    fixture=make_fixture(root/"synthetic_rho_xy.h5",**config["fixture"])
    with h5py.File(fixture["path"]) as h5:
        if h5.attrs.get("synthetic_only") != True:
            raise ValueError("Input lacks synthetic-only guard")
    products={}
    for dataset in ("rho_x","rho_y"):
        folder=root/f"ape-synthetic-ps-{dataset}"
        plan=plan_short_synthetic(
            executable=config["ps_executable"],
            ps_executable=config["ps_executable"],
            ps_sha256=config["ps_sha256"],
            input_hdf5=fixture["path"],dataset=dataset,
            fft_hdf5=folder/"shared_fft.h5",output_dir=folder)
        commands=run_local(plan,synthetic_approved=True)
        batch=finalize_ps_plan(plan)
        ps_path=folder/f"{dataset}.combined.ps"
        n=join_ps_sidecars(batch,ps_path)
        if n!=1301:
            raise RuntimeError(f"Expected 1301 native C PS bins; got {n}")
        products[dataset]=dict(path=str(ps_path),jobs=len(commands),bins=n)
    a=np.loadtxt(products["rho_x"]["path"],ndmin=2)
    b=np.loadtxt(products["rho_y"]["path"],ndmin=2)
    if a.shape!=(1301,2) or b.shape!=(1301,2) or not np.array_equal(a[:,0],b[:,0]):
        raise ValueError("C output X/Y grids differ")
    f=a[:,0]
    observed=np.vstack([a[:,1],b[:,1]])
    if np.any(observed<=0):
        raise ValueError("Nonpositive PS cannot enter logarithmic vendor fitter")
    fitted=[fit_window(f,observed[q],config["external_vendor_root"]) for q in range(2)]
    if not all(x["vendor_quality_pass"] for x in fitted):
        raise RuntimeError("Original vendor background quality gate failed")
    shape=np.array([1.]) if config.get("peak_shape_file") is None else load_vendor_peak_shape(
        config["peak_shape_file"],frequency_bins_hz=f)
    if shape.size>1301:
        raise ValueError("Signal shape exceeds 1301-bin smoke grid")
    indices=np.asarray(config["candidate_start_indices"],dtype=int)
    if (indices.ndim!=1 or indices.size==0 or np.any(indices<0)
            or np.any(indices+len(shape)>len(f)) or np.any(np.diff(indices)<=0)):
        raise ValueError("Candidate indices must be distinct sorted valid windows")
    freq_windows=np.stack([f[i:i+len(shape)] for i in indices])
    obs_windows=np.stack([observed[:,i:i+len(shape)] for i in indices])
    bg_windows=np.stack([np.stack([fit["background_psd"][i:i+len(shape)]
                                      for fit in fitted]) for i in indices])
    skew=np.stack([np.stack([fit["skewnormal_loc_scale_alpha"] for fit in fitted])
                   for _ in indices])
    scan=scan_candidate_windows(
        candidate_frequencies_hz=f[indices],
        observed_ps=obs_windows,background_ps=bg_windows,
        skewnormal_parameters=skew,frequency_windows_hz=freq_windows,
        vendor_peak_shape=shape,quadrature_power_fraction=config["quadrature_power_fraction"],
        reference_coupling_gev_inv=config["reference_coupling_gev_inv"],
        coupling_grid_gev_inv=config["coupling_grid_gev_inv"],
        external_root=config["external_vendor_root"],
        local_sigma_threshold=config["local_sigma_threshold"],
        finesse_s=config["finesse_s"],tau_s=config["tau_s"],
        tau1_s=config["tau1_s"],rho_dm_gev4=config["rho_dm_gev4"],
        unit_conversion=config["unit_conversion"],
        n_effective_trials=config.get("n_effective_trials"))
    manifest=dict(synthetic_only=True,scientific_status="software_smoke_only_not_physical_coverage",
                  source="native_C_PS_not_Python_spectral_estimator",
                  binary_sha256=config["ps_sha256"],
                  input=fixture,products=products,scan=scan,
                  vendor_source_sha256=fitted[0]["vendor_source_sha256"],
                  fit_quality=[x["vendor_reduced_chi2"] for x in fitted],
                  shape_source=config.get("peak_shape_file","original_vendor_single_bin_default"))
    (root/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    with (root/"candidates.csv").open("x",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=["frequency_hz","q0","local_sqrt_q0",
                    "mu_hat","cluster_start_index","cluster_end_index_exclusive","status"])
        writer.writeheader()
        writer.writerows(scan["candidates"])
    return manifest


def main():
    parser=argparse.ArgumentParser(description="Explicit synthetic-only APE C PS and vendor candidate smoke run")
    parser.add_argument("--config",required=True,type=Path)
    args=parser.parse_args()
    config=json.loads(args.config.read_text())
    result=run(config)
    print(json.dumps(dict(result="software_smoke_completed",candidate_count=len(result["scan"]["candidates"]),
                           run_dir=config["run_dir"]),indent=2))


if __name__=="__main__":
    main()
