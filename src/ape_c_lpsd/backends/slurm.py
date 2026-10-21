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


def render_slurm_array_preview(plan, *, modules=(), max_parallel=4):
    """One FFT job followed by one Slurm array, not 16 individual submissions.

    Generates scripts and documented sbatch commands. Never calls sbatch.
    Module names are site-specific and must be explicitly supplied on VIPER.
    """
    if not 1 <= max_parallel <= 16:
        raise ConstantQError("max_parallel must be between 1 and 16")
    entries = plan["batches"]
    if len(entries) != 16 or [i for i, _ in entries] != list(range(16)):
        raise ConstantQError("Expected the frozen 16 indexed batch commands")
    output = Path(plan["output_dir"])
    if not output.is_absolute():
        raise ConstantQError("Output directory must be absolute")
    module_lines = "".join("module load " + shlex.quote(name) + "\n" for name in modules)
    config = "export LPSDCFN=" + shlex.quote(str(output / "lpsd.cfg")) + "\n"
    header = "#!/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n"
    fft = (header + "#SBATCH --mem=32G\n#SBATCH --time=02:00:00\n"
           "set -euo pipefail\n" + module_lines + config
           + shlex.join(plan["fft"]) + "\n")
    cases = "".join("  " + str(i) + ") " + shlex.join(argv) + " ;;\n"
                    for i, argv in entries)
    array = (header + "#SBATCH --mem=16G\n#SBATCH --time=08:00:00\n"
             + f"#SBATCH --array=0-15%{max_parallel}\n"
             + "set -euo pipefail\n" + module_lines + config
             + 'case "${SLURM_ARRAY_TASK_ID:?Missing Slurm array index}" in\n'
             + cases + '  *) echo "Invalid array index" >&2; exit 2 ;;\n'
             + "esac\n")
    return {
        "fft_script": fft, "array_script": array,
        "fft_submission": "sbatch --parsable fft.sbatch",
        "array_submission": "sbatch --dependency=afterok:<FFT_JOB_ID> array.sbatch",
        "batch_count": 16, "max_parallel": max_parallel, "submitted": False,
        "site_modules_confirmed": bool(modules),
    }
