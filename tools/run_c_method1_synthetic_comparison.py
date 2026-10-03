#!/usr/bin/env python3
"""Synthetic-only direct C-LPSD versus C method-1 FFT-approximation study.

This tool deliberately creates new generated HDF5 files in a new results
directory.  It never reads an APE cache, raw Dewesoft file, or real export.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
import uuid

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ape_c_lpsd.c_output_validation import validate_frequency_rows
from ape_c_lpsd.core import sha256


FS_HZ = 512.0
SAMPLES = 500  # deliberately not a power of two; no real-data relevance
JDES = 24
EPSILON_PERCENT = 10.0
PROFILE = {"window": -2, "psll_db": 200, "linear_regression": 1,
           "overlap_percent": 50, "desired_averages": 100,
           "minimum_averages": 1, "input_scale": 1, "max_memory_power": 20}
REGIONS = {"low": (64.0, 128.0), "middle": (128.0, 192.0), "high": (192.0, 240.0)}


def generated_nodes(fmin: float, fmax: float) -> np.ndarray:
    return fmin * np.exp(np.arange(JDES) * math.log(fmax / fmin) / (JDES - 1))


def write_input(path: Path, values: np.ndarray, scenario: str) -> str:
    with h5py.File(path, "x") as handle:
        dataset = handle.create_dataset("/strain", data=values, dtype="f8")
        dataset.attrs["unit"] = "synthetic_arbitrary_unit"
        dataset.attrs["synthetic_only"] = True
        handle.attrs["sample_rate_hz"] = FS_HZ
        handle.attrs["scenario"] = scenario
        handle.attrs["no_real_ape_data"] = True
    return sha256(path)


def run_method(*, exe: Path, directory: Path, input_path: Path, region: str,
               scenario: str, method: int, fmin: float, fmax: float) -> dict:
    stem = f"{region}-{scenario}-h{method}"
    output = directory / f"{stem}.txt"
    config = directory / f"{stem}.cfg"
    config.write_text('GNUTERM "APE synthetic comparison"\nFORMAT "fDNRb"\nEND\n')
    command = [str(exe), "-d", "-D", "/strain", "-f", repr(FS_HZ), "-h", str(method),
               "-i", str(input_path), "-o", str(output), "-b", "0",
               "-e", repr((SAMPLES - 1) / FS_HZ), "-s", repr(fmin), "-t", repr(fmax),
               "-J", str(JDES), "-n", str(JDES), "-N", "0", "-w", "-2", "-p", "200",
               "-r", "1", "-l", "50", "-a", "100", "-m", "1", "-x", "1", "-M", "20"]
    if method == 1:
        command.extend(["-E", repr(EPSILON_PERCENT)])
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    started = time.perf_counter()
    result = subprocess.run(command, text=True, capture_output=True,
                            env={"LPSDCFN": str(config), **os.environ})
    elapsed = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    (directory / f"{stem}.stdout.txt").write_text(result.stdout)
    (directory / f"{stem}.stderr.txt").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"{stem} returned {result.returncode}: {result.stderr or result.stdout}")
    expected = generated_nodes(fmin, fmax)
    raw_rows = []
    for line in output.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            fields = line.split()
            if len(fields) != 5:
                raise RuntimeError(f"{stem}: non-five-column C output row: {line!r}")
            raw_rows.append([float(value) for value in fields])
    rows = np.asarray(raw_rows, dtype=float)
    try:
        _, frequency_validation = validate_frequency_rows(output.read_text().splitlines(), expected)
        frequency_validation["strict_pass"] = True
    except Exception as exc:
        frequency_validation = {"strict_pass": False, "error": str(exc)}
    np.save(directory / f"{stem}.npy", rows)
    return {"method": method, "command": command, "elapsed_s": elapsed,
            "child_maxrss_kb_before": before, "child_maxrss_kb_after": after,
            "output": output.name, "output_sha256": sha256(output),
            "stdout": f"{stem}.stdout.txt", "stderr": f"{stem}.stderr.txt",
            "frequency_hz": rows[:, 0].tolist(), "psd": rows[:, 1].tolist(),
            "averages": rows[:, 2].astype(int).tolist(),
            "resolution_bandwidth_hz": rows[:, 3].tolist(),
            "bin_number": rows[:, 4].astype(int).tolist(), "frequency_text_validation": frequency_validation}


def relative_difference(reference: np.ndarray, other: np.ndarray) -> dict:
    scale = np.maximum(np.abs(reference), np.finfo(float).tiny)
    values = np.abs(other - reference) / scale
    return {"max": float(values.max()), "median": float(np.median(values)),
            "p95": float(np.percentile(values, 95))}


def comparison(direct: dict, approx: dict, tone_hz: float | None) -> dict:
    d = np.asarray(direct["frequency_hz"]); a = np.asarray(approx["frequency_hz"])
    dp = np.asarray(direct["psd"]); ap = np.asarray(approx["psd"])
    same_nodes = bool(np.array_equal(d, a))
    answer = {"same_returned_nodes_exactly": same_nodes,
              "max_frequency_difference_hz": float(np.max(np.abs(d - a))),
              "psd_relative_difference": relative_difference(dp, ap),
              "runtime_speedup_direct_over_approx": float(direct["elapsed_s"] / approx["elapsed_s"]),
              "direct_peak_frequency_hz": float(d[np.argmax(dp)]),
              "approx_peak_frequency_hz": float(a[np.argmax(ap)])}
    if tone_hz is not None:
        answer["injected_tone_hz"] = tone_hz
        answer["direct_peak_offset_hz"] = answer["direct_peak_frequency_hz"] - tone_hz
        answer["approx_peak_offset_hz"] = answer["approx_peak_frequency_hz"] - tone_hz
    return answer


def main() -> None:
    exe = Path("/Users/user/Downloads/Scalar-Dark-Matter-LPSD-1.0.0/LPSD/manual-build/lpsd-exec")
    if not exe.is_file():
        raise SystemExit(f"missing lpsd-exec: {exe}")
    directory = ROOT / "results" / "synthetic-method1-comparisons" / f"run-{uuid.uuid4().hex}"
    directory.mkdir(parents=True, exist_ok=False)
    t = np.arange(SAMPLES, dtype=np.float64) / FS_HZ
    all_results: dict[str, object] = {"synthetic_only": True, "no_real_ape_data": True,
      "executable": str(exe), "executable_sha256": sha256(exe), "sample_rate_hz": FS_HZ,
      "sample_count": SAMPLES, "jdes": JDES, "method1_epsilon_percent": EPSILON_PERCENT,
      "matched_profile": PROFILE, "regions": {}, "run_directory": str(directory)}
    for index, (region, (fmin, fmax)) in enumerate(REGIONS.items()):
        nodes = generated_nodes(fmin, fmax)
        tone_hz = float(nodes[JDES // 2])
        noise = np.random.default_rng(7300 + index).normal(0.0, 0.25, SAMPLES)
        scenarios = {"tone": np.sin(2.0 * np.pi * tone_hz * t),
                     "white_noise": noise,
                     "mixed": np.sin(2.0 * np.pi * tone_hz * t) + noise}
        region_result: dict[str, object] = {"limits_hz": [fmin, fmax], "expected_nodes_hz": nodes.tolist(),
                                             "tone_hz": tone_hz, "scenarios": {}}
        for scenario, values in scenarios.items():
            input_path = directory / f"{region}-{scenario}.h5"
            input_hash = write_input(input_path, values, scenario)
            direct = run_method(exe=exe, directory=directory, input_path=input_path, region=region,
                                scenario=scenario, method=0, fmin=fmin, fmax=fmax)
            approx = run_method(exe=exe, directory=directory, input_path=input_path, region=region,
                                scenario=scenario, method=1, fmin=fmin, fmax=fmax)
            region_result["scenarios"][scenario] = {"synthetic_input": input_path.name,
                "synthetic_input_sha256": input_hash, "direct_h0": direct, "fft_approx_h1": approx,
                "comparison": comparison(direct, approx, tone_hz if scenario != "white_noise" else None)}
        all_results["regions"][region] = region_result
    (directory / "run_record.json").write_text(json.dumps(all_results, indent=2))
    print(json.dumps(all_results, indent=2))


if __name__ == "__main__":
    main()
