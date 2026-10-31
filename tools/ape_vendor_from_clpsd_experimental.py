"""Experimental read-only APE PS->vendor workflow coordinator.

The vendor driver is not reimplemented. Original vendor C-LPSD is upstream;
this adapter combines its existing per-observation X/Y native power outputs
and delegates all downstream analysis to the experimental vendor-stage runner.
Inputs are never modified. Outputs must not already exist.

The stage runner is still unvalidated; this is NOT production-ready.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

def main(config):
    cfg = json.loads(Path(config).resolve(strict=True).read_text())
    observations = cfg["observations"]
    if not observations:
        raise ValueError("At least one real observation is required")
    out = Path(cfg["combined_lpsd_directory"]).expanduser().resolve()
    if out.exists():
        raise FileExistsError("Combined directory already exists; no overwrites")
    out.mkdir(parents=True)
    try:
        from ape_vendor_combine_quadrature_power import prepare
        for obs in observations:
            prefix = cfg["data_prefix"]
            t0, t1 = int(obs["t0"]), int(obs["t1"])
            name = f"{prefix}_{t0}_{t1}_APE.h5"
            prepare(
                rho_x_ps=obs["rho_x_ps"], rho_y_ps=obs["rho_y_ps"],
                output_h5=out / name, fs_hz=float(obs["fs_hz"]),
                t0=t0, t1=t1, provenance_json=out / (name + ".json"))
        # Vendor original stages are called from a single adapter entry point.
        from ape_vendor_original_stages_experimental import run
        stage_status = run(config)
        from ape_vendor_single_channel_candidates import select_candidates
        status = json.loads(Path(stage_status).read_text())
        result_dir = Path(cfg["output_directory"]).resolve()
        candidates = select_candidates(
            result_dir / "vendor_mle_result.npy",
            result_dir / "ape_candidates.csv",
            sigma_limit=float(cfg.get("sigma_limit", 5)))
        status["candidate_selection"] = candidates
        status["not_implemented"] = [
            x for x in status["not_implemented"]
            if x != "verified vendor candidate workflow for one APE channel"
            and x != "automatic native C-LPSD to combined vendor HDF5 wiring"
        ]
        status["not_implemented"].extend([
            "C-LPSD execution is external: original C-LPSD PS outputs are inputs",
            "full end-to-end parity and response calibration remain unverified"
        ])
        Path(stage_status).write_text(json.dumps(status, indent=2))
        print(json.dumps(status, indent=2))
    except Exception:
        # Keep partial outputs for audit; never imply successful completion.
        raise

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    main(p.parse_args().config)
