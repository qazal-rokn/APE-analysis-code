"""Minimal Configuration A: X+Y PSD -> ORIGINAL background fit -> ORIGINAL likelihood.

No vendor C/Python scientific functions are modified. The only experiment
substitutions are the combined APE PSD and its supplied signal-power response.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

from configuration_A.combine_quadrature_psd import combine_psd_hdf5
from ape_post_lpsd.ape_max_likelihood import run_single_quadrature
from ape_post_lpsd.ape_signal_response import save_ape_q_response


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="step", required=True)

    bkg = commands.add_parser("background", help="Sum X/Y PSD and fit with ORIGINAL vendor code")
    bkg.add_argument("--x-h5", required=True)
    bkg.add_argument("--y-h5", required=True)
    bkg.add_argument("--output-dir", required=True)
    bkg.add_argument("--ana-fmin", required=True, type=float)
    bkg.add_argument("--ana-fmax", required=True, type=float)
    bkg.add_argument("--n-processes", type=int)
    bkg.add_argument("--segment-size", type=int)

    response = commands.add_parser("response", help="APE Q beta from existing two-QWP response and ORIGINAL scalar linewidth")
    response.add_argument("--spectrum", required=True)
    response.add_argument("--output-npz", required=True)
    response.add_argument("--finesse-s", type=float, default=1500.0)
    response.add_argument("--tau-s", type=float, default=2.5e-9)
    response.add_argument("--tau1-s", type=float, default=1.0e-8)
    response.add_argument("--rho-dm-gev4", type=float, default=3.1e-42)
    response.add_argument("--unit-conversion", type=float, default=1.519e24)
    response.add_argument("--reference-coupling-gev-inv", type=float, default=1e-10)

    lkl = commands.add_parser("likelihood", help="Fit summed Q spectrum using ORIGINAL vendor optimizer")
    lkl.add_argument("--spectrum", required=True)
    lkl.add_argument("--background-json", required=True)
    lkl.add_argument("--response-npz", required=True)
    lkl.add_argument("--output-npy", required=True)
    lkl.add_argument("--ana-fmin", required=True, type=float)
    lkl.add_argument("--ana-fmax", required=True, type=float)
    lkl.add_argument("--min-log10mu", required=True, type=int)
    lkl.add_argument("--max-log10mu", required=True, type=int)
    lkl.add_argument("--start-index", type=int, default=0)
    lkl.add_argument("--max-frequencies", type=int)
    args = parser.parse_args()

    if args.step == "background":
        if not 0 < args.ana_fmin < args.ana_fmax:
            parser.error("Require 0 < ana-fmin < ana-fmax")
        if args.n_processes is not None and args.n_processes < 1:
            parser.error("n-processes must be positive")
        if args.segment_size is not None and args.segment_size < 5:
            parser.error("segment-size must be at least 5")
        outdir = Path(args.output_dir).expanduser().resolve()
        if (outdir / "background_Q.json").exists():
            parser.error("Refusing to overwrite previous Q background fitting")
        combined = combine_psd_hdf5(args.x_h5, args.y_h5, outdir)
        command = [
            sys.executable, "-m", "scalar_dark_matter_lpsd.lvk.fit_background",
            "--data-path", combined["output_file"],
            "--output-path", str(outdir / "background_Q.json"),
            "--ana-fmin", str(args.ana_fmin),
            "--ana-fmax", str(args.ana_fmax),
        ]
        if args.n_processes is not None:
            command += ["--n-processes", str(args.n_processes)]
        if args.segment_size is not None:
            command += ["--segment-size", str(args.segment_size)]
        subprocess.run(command, check=True)
        print(json.dumps({"combined_psd":combined["output_file"],
                          "background_json":str(outdir / "background_Q.json"),
                          "vendor_background_command":command}, indent=2))
    elif args.step == "response":
        result = save_ape_q_response(
            args.spectrum, args.output_npz,
            finesse_s=args.finesse_s, tau_s=args.tau_s,
            tau1_s=args.tau1_s, rho_dm_gev4=args.rho_dm_gev4,
            unit_conversion=args.unit_conversion,
            reference_coupling_gev_inv=args.reference_coupling_gev_inv,
        )
        print(json.dumps(result, indent=2))
    else:
        if not Path(args.spectrum).name.endswith("_Q.h5"):
            parser.error("Configuration A likelihood requires the summed Q spectrum")
        result = run_single_quadrature(
            args.spectrum, args.background_json, args.response_npz,
            args.output_npy, channel="Q",
            ana_fmin=args.ana_fmin, ana_fmax=args.ana_fmax,
            min_log10mu=args.min_log10mu, max_log10mu=args.max_log10mu,
            start_index=args.start_index, max_frequencies=args.max_frequencies,
        )
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
