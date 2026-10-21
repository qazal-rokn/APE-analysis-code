"""Explicit local Linux execution. No implicit real-data run."""
from __future__ import annotations
import json
import os
import subprocess
from pathlib import Path
from ape_c_lpsd.spectral.constantq import ConstantQError, plan_manifest, sha256_file


def prepare_run(plan):
    out = Path(plan["output_dir"])
    out.mkdir(parents=True, exist_ok=False)
    config = out/"lpsd.cfg"
    config.write_text(plan["c_config"], encoding="utf-8")
    (out/"commands.json").write_text(json.dumps(plan_manifest(plan), indent=2)+"\n")
    return config


def run_local(plan, *, synthetic_approved=False):
    """Execute plan only after explicit synthetic-run authorization.

    Production real-data execution is deliberately excluded from this adapter.
    Input must be verified synthetic by the operator; this function does not
    infer authenticity from a filename.
    """
    if not synthetic_approved:
        raise ConstantQError("REFUSED: explicit synthetic execution permission required")
    executable = Path(plan["fft"][0])
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ConstantQError("External lpsd-exec missing or non-executable")
    cfg = prepare_run(plan)
    env = os.environ.copy()
    env["LPSDCFN"] = str(cfg)
    results = []
    for label, argv in [("fft", plan["fft"]), *[(f"batch_{i}", a) for i, a in plan["batches"]]]:
        process = subprocess.run(argv, cwd=plan["output_dir"], env=env,
                                 text=True, capture_output=True, check=False)
        out = Path(plan["output_dir"])
        (out/f"{label}.stdout.txt").write_text(process.stdout)
        (out/f"{label}.stderr.txt").write_text(process.stderr)
        results.append({"step": label, "return_code": process.returncode, "argv": argv})
        (out/"execution.json").write_text(json.dumps({"binary_sha256": sha256_file(executable),
                         "commands": results}, indent=2))
        if process.returncode != 0:
            raise ConstantQError(f"{label} failed with exit status {process.returncode}")
    return results
