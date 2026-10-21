"""APE Constant-Q controller: plan-only safe entry point.

No real data or cluster job is ever executed from this CLI.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from ape_c_lpsd.spectral.constantq import plan_commands, plan_manifest
from ape_c_lpsd.backends.slurm import render_slurm_preview, render_slurm_array_preview


def main():
    parser = argparse.ArgumentParser(description="Read-only APE Constant-Q command planner")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--input-hdf5", type=Path, required=True)
    parser.add_argument("--dataset", choices=("rho_x", "rho_y"), required=True)
    parser.add_argument("--fft-hdf5", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--backend", choices=("local", "slurm"), default="local")
    parser.add_argument("--slurm-array", action="store_true", help="Preview one FFT job plus one Slurm array")
    parser.add_argument("--max-parallel", type=int, default=4, help="Concurrent Slurm array tasks (1-16)")
    args = parser.parse_args()
    plan = plan_commands(args.executable, args.input_hdf5, args.dataset,
                         args.fft_hdf5, args.output_dir)
    result = {"plan": plan_manifest(plan), "executed": False, "submitted": False}
    if args.backend == "slurm":
        result["slurm_preview"] = (render_slurm_array_preview(plan, max_parallel=args.max_parallel)
                                   if args.slurm_array else render_slurm_preview(plan))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
