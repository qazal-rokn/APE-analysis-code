"""Immutable historical APE Constant-Q command profile and batch plan.

This module DOES NOT estimate a spectrum in Python. Executable commands use the
original, unchanged C lpsd-exec; no real-data execution is authorized here.
"""
from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

ENDPOINTS = (0, 3649978, 3996383, 4199059, 4342872, 4454427,
             4545576, 4622643, 4689403, 4748290, 4800966,
             4848618, 4892122, 4932141, 4969193, 5003688, 5035956)
FREQUENCIES = 5035956


class ConstantQError(ValueError):
    """The proposed C run violates a recorded APE invariant."""


@dataclass(frozen=True)
class ConstantQProfile:
    fs: float = 20000.0
    t_start: float = 0.0
    t_stop: float = 15435.0
    f_min: float = 65.0
    f_max: float = 10000.0
    grid_bins: int = FREQUENCIES
    overlap: int = 20
    psll_db: float = 238.13
    memory_power: int = 26
    window: int = -2
    regression: int = 0
    scale: int = 1
    epsilon: int = 10
    kernel_threshold: float = 1e-15
    a_option: int = 2

    def validate_reference(self):
        if self != ConstantQProfile():
            raise ConstantQError("Production profile differs from frozen historical settings")
        if len(ENDPOINTS) != 17 or ENDPOINTS[0] != 0 or ENDPOINTS[-1] != self.grid_bins:
            raise ConstantQError("Incorrect Constant-Q partition")
        if any(b <= a for a, b in zip(ENDPOINTS, ENDPOINTS[1:])):
            raise ConstantQError("Nonincreasing Constant-Q batch endpoints")
        return self


def batches(profile=ConstantQProfile()):
    profile.validate_reference()
    return tuple((i, ENDPOINTS[i], ENDPOINTS[i+1] - ENDPOINTS[i]) for i in range(16))


def c_format() -> str:
    # C output method-2 columns: frequency, PSD, N_avg, N_kernel, real, imaginary.
    return 'GNUTERM "APE Constant-Q"\nFORMAT "fDNbri"\nEND\n'


def plan_commands(executable, input_hdf5, dataset, fft_hdf5, output_dir,
                  profile=ConstantQProfile()):
    """Generate one shared FFT and sixteen Constant-Q commands for ONE quadrature.

    Paths must be absolute, and the caller must create a new, isolated run dir.
    The FFT file is shared by the 16 batches; it is NOT recomputed per batch.
    The historical flags are intentionally explicit and frozen.
    """
    profile.validate_reference()
    if dataset not in ("rho_x", "rho_y"):
        raise ConstantQError("Dataset must be rho_x or rho_y; no time-domain magnitude")
    exe, inp, fft, out = map(lambda p: Path(p).expanduser(),
                              (executable, input_hdf5, fft_hdf5, output_dir))
    if not all(p.is_absolute() for p in (exe, inp, fft, out)):
        raise ConstantQError("Executable/input/FFT/output paths must be absolute")
    if out == inp or out == fft or out in inp.parents or out in exe.parents:
        raise ConstantQError("Unsafe output directory")
    p = profile
    common = [str(exe), "-d", "-D", "/"+dataset,
              "-f", str(p.fs), "-i", str(inp),
              "-b", str(p.t_start), "-e", str(p.t_stop),
              "-s", str(p.f_min), "-t", str(p.f_max),
              "-J", str(p.grid_bins), "-A", str(p.a_option),
              "-l", str(p.overlap), "-r", str(p.regression),
              "-w", str(p.window), "-p", str(p.psll_db),
              "-x", str(p.scale), "-E", str(p.epsilon),
              "-M", str(p.memory_power), "-R", str(p.kernel_threshold),
              "-T", "-O", str(fft)]
    # A new FFT is generated precisely once and then reused by method 2.
    fft_command = [*common, "-h", "3", "-N", "0", "-n", str(p.grid_bins),
                   "-o", str(out/"fft_stage.lpsd.0")]
    constq = [(i, [*common, "-h", "2", "-N", str(start), "-n", str(n),
                   "-o", str(out/f"{dataset}.lpsd.{i}")])
              for i, start, n in batches(profile)]
    return {"dataset": dataset, "fft": fft_command, "batches": constq,
            "output_dir": str(out), "fft_hdf5": str(fft),
            "c_config": c_format(), "format": "fDNbri"}


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def plan_manifest(plan):
    return {"dataset": plan["dataset"], "fft": plan["fft"],
            "batches": [{"index": i, "argv": argv} for i, argv in plan["batches"]],
            "output_dir": plan["output_dir"], "format": plan["format"]}
