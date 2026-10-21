"""Safe Slurm preview generator. Never submits jobs to the cluster."""
from __future__ import annotations
import shlex
from pathlib import Path
from ape_c_lpsd.spectral.constantq import ConstantQError


def render_slurm_preview(plan, *, modules=("gcc/14", "hdf5-serial/1.14.1", "gsl/2.7")):
    """Return dependency-aware previews with identical scientific argv.

    Slurm reads #SBATCH directives only before the first shell command.
    The caller must prepare a fresh run directory and LPSDCFN before submission.
    Nothing here runs sbatch, writes cluster files, or claims Slurm validation.
    """
    module_setup = "".join(f"module load {shlex.quote(x)}\n" for x in modules)
    config = str(Path(plan["output_dir"])/"lpsd.cfg")
    env = f"export LPSDCFN={shlex.quote(config)}\n"
    def script(argv, memory, time):
        return (f"#!/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n"
                f"#SBATCH --mem={memory}\n#SBATCH --time={time}\n"
                "set -euo pipefail\n" + module_setup + env + shlex.join(argv) + "\n")
    fft = script(plan["fft"], "32G", "02:00:00")
    batch_scripts = {i: script(argv, "16G", "08:00:00")
                     for i, argv in plan["batches"]}
    return {"fft": fft, "batches": batch_scripts,
            "required_dependency": "submit each batch with --dependency=afterok:<FFT_JOB_ID>",
            "submitted": False}


def write_slurm_preview(plan):
    raise ConstantQError("Slurm preview is returned in memory; no cluster writes/submission authorized")
