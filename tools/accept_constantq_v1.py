"""One-command synthetic acceptance: real C FFT, 3 batches, merger, reference, Slurm preview.

No production/experimental C jobs; no cluster submission. Runs only when --execute is given.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ape_c_lpsd.spectral.constantq import plan_synthetic_fixture_commands, plan_commands
from ape_c_lpsd.spectral.combine import validate_and_combine
from ape_c_lpsd.backends.slurm import render_slurm_preview

DOCKER = "/Applications/Docker.app/Contents/Resources/bin/docker"
CONTAINER = "lpsd-linux"
EXEC = "/work/lpsd-project/LPSD/install/bin/lpsd-exec"
INPUT = "/tmp/ape_validation_noise_voltage.h5"
REFERENCE = "/tmp/ape_lpsd_validation_q250/noise_constq.lpsd"
PARTITION = (0, 400, 850, 1301)
KNOWN_SHA = "6c47789484cde77d72c959cb807b4dafb9ffe1fdf6f73e0aed36a52f31d3f62c"


def docker(*args, input=None):
    return subprocess.run([DOCKER, *args], text=True, input=input,
                          capture_output=True, check=True).stdout


def spectrum(text, expected):
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        vals = tuple(float(x) for x in line.split())
        if len(vals) < 4 or not all(math.isfinite(v) for v in vals):
            raise AssertionError("Invalid method-2 spectrum columns")
        if vals[0] <= 0 or vals[1] < 0 or vals[2] < 0 or vals[3] < 0:
            raise AssertionError("Invalid spectral values")
        if rows and vals[0] <= rows[-1][0]:
            raise AssertionError("Frequencies not strictly increasing")
        rows.append(vals)
    if len(rows) != expected:
        raise AssertionError(f"Expected {expected} frequency bins, got {len(rows)}")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not Path(DOCKER).is_file():
        raise RuntimeError("Docker CLI missing")
    docker("exec", CONTAINER, "test", "-x", EXEC)
    docker("exec", CONTAINER, "test", "-s", INPUT)
    docker("exec", CONTAINER, "test", "-s", REFERENCE)
    digest = docker("exec", CONTAINER, "sha256sum", INPUT).split()[0]
    if digest != KNOWN_SHA:
        raise RuntimeError("Synthetic fixture checksum mismatch; no C run permitted")

    production = plan_commands(
        EXEC, "/tmp/placeholder-calibrated-ape.h5", "rho_x",
        "/tmp/preview-fft.h5", "/tmp/preview-ape-constantq")
    preview = render_slurm_preview(production)
    if len(preview["batches"]) != 16 or preview["submitted"] is not False:
        raise AssertionError("Slurm plan invalid")
    if "--dependency=afterok:<FFT_JOB_ID>" not in preview["required_dependency"]:
        raise AssertionError("Slurm FFT dependency not declared")
    for job in (preview["fft"], *preview["batches"].values()):
        if job.index("#SBATCH") > job.index("set -euo pipefail"):
            raise AssertionError("Invalid Slurm directive order")

    print("Frozen production profile: preserved; 16 Slurm scripts planned, never submitted", flush=True)
    print("Synthetic fixture SHA256 verified:", digest, flush=True)
    if not args.execute:
        print("PREFLIGHT ONLY. Add --execute for real synthetic FFT + 3 C batches.", flush=True)
        return

    work = docker("exec", CONTAINER, "mktemp", "-d",
                  "/tmp/ape-controller-v1-accept-XXXXXXXX").strip()
    if not re.fullmatch(r"/tmp/ape-controller-v1-accept-[A-Za-z0-9]{8}", work):
        raise RuntimeError("Unsafe synthetic directory")
    plan = plan_synthetic_fixture_commands(
        EXEC, INPUT, work + "/fft_cache.h5", work, endpoints=PARTITION)
    print("New isolated synthetic directory:", work, flush=True)
    steps = [("fft", plan["fft"]), *[(f"batch-{i}", c) for i, c in enumerate(plan["batches"])]]
    with tempfile.TemporaryDirectory(prefix="ape-constantq-accept-") as d:
        local = Path(d)
        batch_paths = []
        for label, argv in steps:
            print("Executing", label, "through shared C planner", flush=True)
            process = subprocess.run([DOCKER, "exec", CONTAINER, *argv],
                                     capture_output=True, text=True)
            for ext, contents in (("stdout", process.stdout), ("stderr", process.stderr)):
                docker("exec", "-i", CONTAINER, "tee", f"{work}/{label}.{ext}.txt", input=contents)
            if process.returncode:
                raise RuntimeError(f"{label} failed: exit {process.returncode}; see {work} logs")
            target = argv[argv.index("-o")+1]
            docker("exec", CONTAINER, "test", "-s", target if label != "fft" else work + "/fft_cache.h5")
            if label != "fft":
                path = local / (label + ".lpsd")
                subprocess.run([DOCKER, "cp", f"{CONTAINER}:{target}", str(path)], check=True)
                batch_paths.append(path)

        merged = local / "combined.lpsd"
        bins = validate_and_combine(batch_paths, merged, endpoints=PARTITION)
        combined = spectrum(merged.read_text(), 1301)
        original = spectrum(docker("exec", CONTAINER, "cat", REFERENCE), 1301)
        freq_match = all(x[0] == y[0] for x, y in zip(combined, original))
        psd_match = all(x[1] == y[1] for x, y in zip(combined, original))
        metadata_match = all(x[2:4] == y[2:4] for x, y in zip(combined, original))
        max_psd_delta = max(abs(x[1]-y[1]) for x, y in zip(combined, original))
        sha_after = docker("exec", CONTAINER, "sha256sum", INPUT).split()[0]
        result = {
            "fft_calls": 1, "method2_batches": 3, "combined_bins": bins,
            "frequency_grid_exact": freq_match, "psd_exact": psd_match,
            "metadata_exact": metadata_match,
            "max_absolute_psd_difference": max_psd_delta,
            "synthetic_input_unchanged": sha_after == digest,
            "slurm_scripts_previewed": 17, "slurm_jobs_submitted": 0,
            "production_ape_data_accessed": False, "logs_in_container": work,
        }
        print(json.dumps(result, indent=2))
        if not all((freq_match, psd_match, metadata_match, sha_after == digest)):
            raise RuntimeError("ACCEPTANCE FAILED: combined output differs from reference")
        print("ACCEPTANCE PASSED: local C planning, 3 batches, combination, reference, Slurm preview")


if __name__ == "__main__":
    main()
