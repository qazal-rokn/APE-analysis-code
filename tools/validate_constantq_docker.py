"""Isolated Docker smoke test of the previously validated 60 s synthetic C workflow.

The synthetic profile stays separate from the frozen production profile,
but both use the same C argument-building function.
No raw/real APE files, production spectra, or original vendor sources are changed.

Default is a READ-ONLY preflight. --execute is explicit authorization to run
the existing synthetic C input in a fresh Docker /tmp directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from ape_c_lpsd.spectral.constantq import plan_synthetic_fixture_commands

DOCKER = Path("/Applications/Docker.app/Contents/Resources/bin/docker")
CONTAINER = "lpsd-linux"
EXEC = "/work/lpsd-project/LPSD/install/bin/lpsd-exec"
INPUT = "/tmp/ape_validation_noise_voltage.h5"
ORIGINAL = "/tmp/ape_lpsd_validation_q250/noise_constq.lpsd"
EXPECTED_BINS = 1301


def docker(*args, capture=True):
    return subprocess.run([str(DOCKER), *args], check=True, text=True,
                          capture_output=capture)


def verify_spectrum(text):
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 4:
            raise RuntimeError("C output missing expected PSD/N_avg/kernel columns")
        try:
            nums = [float(x) for x in cols]
        except ValueError as exc:
            raise RuntimeError("Non-numeric C output") from exc
        import math
        if not all(math.isfinite(x) for x in nums) or nums[0] <= 0 or nums[1] < 0:
            raise RuntimeError("Invalid C spectrum")
        if rows and nums[0] <= rows[-1][0]:
            raise RuntimeError("Nonmonotonic C frequency grid")
        rows.append(nums)
    if len(rows) != EXPECTED_BINS:
        raise RuntimeError(f"Expected {EXPECTED_BINS} bins, got {len(rows)}")
    return rows


def main():
    p = argparse.ArgumentParser(description="APE isolated synthetic Docker C smoke test")
    p.add_argument("--execute", action="store_true",
                   help="Explicitly run on synthetic fixture; otherwise preflight only")
    a = p.parse_args()

    if not DOCKER.is_file():
        raise SystemExit(f"Docker CLI not found: {DOCKER}")
    docker("exec", CONTAINER, "test", "-x", EXEC)
    docker("exec", CONTAINER, "test", "-s", INPUT)
    digest_before = docker("exec", CONTAINER, "sha256sum", INPUT).stdout.split()[0]
    if not re.fullmatch(r"[0-9a-f]{64}", digest_before):
        raise RuntimeError("Could not verify input checksum")
    print("External executable: present")
    print("Synthetic input SHA-256:", digest_before)
    print("Frozen synthetic profile: 60 s, 2048 Hz, 5-900 Hz, 1301 bins")
    if not a.execute:
        print("PREFLIGHT ONLY; no computations run. Use --execute for new synthetic run.")
        return

    # Use a unique directory. Never reuse the historical results directory.
    dirname = docker("exec", CONTAINER, "mktemp", "-d",
                     "/tmp/ape-controller-v1-XXXXXXXX").stdout.strip()
    if not re.fullmatch(r"/tmp/ape-controller-v1-[A-Za-z0-9]{8}", dirname):
        raise RuntimeError("Unexpected mktemp directory")
    fft = f"{dirname}/fft_cache.h5"
    first = f"{dirname}/fft_stage.lpsd"
    second = f"{dirname}/noise_constq.lpsd"
    plan = plan_synthetic_fixture_commands(EXEC, INPUT, fft, dirname)
    steps = [(3, plan["fft"]), (2, plan["batch"])]
    print("NEW synthetic work directory:", dirname, flush=True)
    for method, argv in steps:
        print("Running original C method", method, flush=True)
        proc = subprocess.run([str(DOCKER), "exec", CONTAINER, *argv],
                              capture_output=True, text=True)
        # Logs are preserved only inside the new synthetic directory.
        for label, value in (("stdout", proc.stdout), ("stderr", proc.stderr)):
            subprocess.run([str(DOCKER), "exec", "-i", CONTAINER,
                            "tee", f"{dirname}/method{method}.{label}.txt"],
                           input=value, text=True, check=True, capture_output=True)
        if proc.returncode:
            raise RuntimeError(f"C method {method} failed (exit {proc.returncode}); "
                               f"logs in {dirname}")
        docker("exec", CONTAINER, "test", "-s", fft if method == 3 else second)

    actual_text = docker("exec", CONTAINER, "cat", second).stdout
    actual = verify_spectrum(actual_text)
    original_available = subprocess.run([str(DOCKER), "exec", CONTAINER,
                                        "test", "-s", ORIGINAL]).returncode == 0
    report = {"synthetic_directory": dirname, "input_sha256": digest_before,
              "c_spectrum_bins": len(actual), "comparison_available": original_available,
              "original_c_unchanged": True}
    if original_available:
        old = verify_spectrum(docker("exec", CONTAINER, "cat", ORIGINAL).stdout)
        report["frequency_grid_exact"] = all(x[0] == y[0] for x, y in zip(actual, old))
        report["psd_exact"] = all(x[1] == y[1] for x, y in zip(actual, old))
        report["max_absolute_psd_difference"] = max(abs(x[1] - y[1])
                                                    for x, y in zip(actual, old))
    digest_after = docker("exec", CONTAINER, "sha256sum", INPUT).stdout.split()[0]
    if digest_after != digest_before:
        raise RuntimeError("Synthetic HDF5 unexpectedly changed")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
