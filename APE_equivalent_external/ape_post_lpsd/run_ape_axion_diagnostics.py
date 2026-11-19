#!/usr/bin/env python3
"""Run the existing selected-frequency APE diagnostic chain without modifying vendor code.

This reproduces the validated 5/8 kHz engineering workflow, NOT a full-band,
physically calibrated or statistically covered axion likelihood scan.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True, help="APE 4.28-hour analysis root containing handoff/ and background/")
    p.add_argument("--output-dir", type=Path, required=True, help="NEW directory for outputs; never overwrite")
    p.add_argument("--finesse-s", type=float, default=1500.)
    p.add_argument("--tau-s", type=float, default=2.5e-9)
    p.add_argument("--tau1-s", type=float, default=1e-8)
    p.add_argument("--rho-dm-gev4", type=float, default=3.1e-42)
    p.add_argument("--unit-conversion", type=float, default=1.519e24)
    p.add_argument("--g-ref-gev-inv", type=float, default=1e-10)
    a = p.parse_args()
    root = a.root.expanduser().resolve()
    out = a.output_dir.expanduser().resolve()
    if out.exists():
        p.error("Output directory already exists. Use a fresh path; no overwrites")
    for channel in ("X", "Y"):
        if not (root / "handoff" / f"APE_0_15435_{channel}.h5").is_file():
            p.error(f"Missing {channel} input HDF5")
        if not (root / "background" / f"background_{channel}_FULL.json").is_file():
            p.error(f"Missing {channel} background JSON")
    out.mkdir(parents=True)
    def run(script, *args):
        cmd = [sys.executable, str(HERE / script), *map(str, args)]
        print("RUN", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True, cwd=str(HERE.parent))
    spectrum_x = root / "handoff" / "APE_0_15435_X.h5"
    spectrum_y = root / "handoff" / "APE_0_15435_Y.h5"
    halo = out / "halo_5k_8k.npz"
    qwp = out / "qwp_5k_8k.npz"
    vendor = out / "vendor_multibin_5k_8k.json"
    run("ape_axion_halo_kernel_diagnostic.py",
        "--spectrum", spectrum_x,
        "--native-kernel-py", HERE / "ape_native_kaiser_kernel.py",
        "--rest-hz", "5000", "--rest-hz", "8000", "--output", halo)
    run("ape_qwp_halo_psd_preview.py",
        "--halo-npz", halo, "--spectrum-x", spectrum_x, "--spectrum-y", spectrum_y,
        "--output", qwp, "--finesse-s", a.finesse_s, "--tau-s", a.tau_s,
        "--tau1-s", a.tau1_s, "--rho-dm-gev4", a.rho_dm_gev4,
        "--unit-conversion", a.unit_conversion, "--g-ref-gev-inv", a.g_ref_gev_inv)
    run("ape_provisional_qwp_multibin_vendor.py",
        "--root", root, "--qwp-npz", qwp, "--output-json", vendor)
    run("ape_signal_weighted_covariance_diagnostic.py",
        "--root", root, "--qwp-npz", qwp,
        "--output-json", out / "covariance_5k_8k.json")
    run("ape_onepass_likelihood_validation.py",
        "--root", root, "--qwp-npz", qwp, "--halo-npz", halo,
        "--output-json", out / "engineering_validation_5k_8k.json")
    run("ape_observed_likelihood_profiles.py",
        "--root", root, "--qwp-npz", qwp,
        "--existing-results-json", vendor,
        "--output-dir", out / "observed_likelihood_profiles")
    print("Complete: selected 5/8 kHz diagnostic chain. NOT PHYSICAL UPPER LIMITS.", flush=True)

if __name__ == "__main__":
    main()
