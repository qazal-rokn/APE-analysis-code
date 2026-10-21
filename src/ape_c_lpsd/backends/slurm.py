"""Safe Slurm preview generator. Never submits jobs to the cluster."""
from __future__ import annotations
import shlex
from pathlib import Path
from ape_c_lpsd.spectral.constantq import ConstantQError


def render_slurm_preview(plan, *, modules=("gcc/14", "hdf5-serial/1.14.1", "gsl/2.7")):
    """Return script text; same argv and scientific settings as local backend.

    Separate jobs make the FFT->method2 dependency explicit. No sbatch call.
    """
    hdr = "#!/bin/bash\nset -euo pipefail\n"
    module_setup = "".join(f"module load {shlex.quote(x)}\n" for x in modules)
    config = str(Path(plan["output_dir"])/"lpsd.cfg")
    env = f"export LPSDCFN={shlex.quote(config)}\n"
    fft = hdr + "#SBATCH --mem=32G\n#SBATCH --time=02:00:00\n" + module_setup + env + shlex.join(plan["fft"])+"\n"
    batches = {}
    for i, argv in plan["batches"]:
        batches[i] = (hdr+"#SBATCH --mem=16G\n#SBATCH --time=08:00:00\n"+
                      module_setup+env+shlex.join(argv)+"\n")
    return {"fft": fft, "batches": batches,
            "required_dependency": "submit each batch with --dependency=afterok:<FFT_JOB_ID>",
            "submitted": False}


def write_slurm_preview(plan):
    raise ConstantQError("Slurm preview is returned in memory; no cluster writes/submission authorized")
