"""Explicit plan-only APE PS controller entry point; no implicit execution.

PS binary must be verified on the machine where it will run. No raw data
is opened; input path is treated as declared synthetic by the operator.
"""
import argparse
import json
from pathlib import Path
from ape_c_lpsd.analysis.ape_pipeline_controller import make_synthetic_ps_plan
from ape_c_lpsd.spectral.constantq import plan_manifest


def main():
    parser=argparse.ArgumentParser(description="APE PS plan-only, no C execution or submission")
    parser.add_argument("--original-executable",type=Path,required=True)
    parser.add_argument("--ps-executable",type=Path,required=True)
    parser.add_argument("--ps-sha256",required=True)
    parser.add_argument("--synthetic-hdf5",type=Path,required=True)
    parser.add_argument("--dataset",choices=["rho_x","rho_y"],required=True)
    parser.add_argument("--fft-hdf5",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    plan=make_synthetic_ps_plan(
        original_executable=args.original_executable,ps_executable=args.ps_executable,
        ps_binary_sha256=args.ps_sha256,synthetic_hdf5=args.synthetic_hdf5,
        dataset=args.dataset,fft_hdf5=args.fft_hdf5,output_dir=args.output_dir)
    print(json.dumps(dict(plan=plan_manifest(plan),ps_enabled=True,
                          ps_binary_sha256=plan["ps_binary_sha256"],
                          submitted=False,executed=False),indent=2))


if __name__=="__main__":
    main()
