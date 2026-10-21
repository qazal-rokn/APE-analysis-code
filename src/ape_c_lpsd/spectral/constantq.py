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



def _c_argv(exe, inp, dataset, fft, output, *, fs, t_stop, f_min, f_max,
            grid_bins, method, start, count, overlap=20, regression=0,
            window=-2, psll=238.13, scale=1, epsilon=10,
            memory_power=26, kernel_threshold=1e-15, a_option=2, debug=False):
    """Shared production/synthetic C argument assembly, with explicit flags."""
    if method not in (2, 3) or start < 0 or count <= 0 or start + count > grid_bins:
        raise ConstantQError("Invalid C method or frequency index interval")
    return [str(exe), *(["-d"] if debug else []), "-A", str(a_option), "-b", "0", "-e", str(t_stop),
            "-f", str(fs), "-h", str(method), "-i", str(inp), "-D", str(dataset),
            "-l", str(overlap), "-n", str(count), "-o", str(output),
            "-r", str(regression), "-s", str(f_min), "-t", str(f_max),
            "-T", "-w", str(window), "-p", str(psll), "-x", str(scale),
            "-N", str(start), "-J", str(grid_bins), "-E", str(epsilon),
            "-M", str(memory_power), "-R", str(kernel_threshold), "-O", str(fft)]


def plan_synthetic_fixture_commands(executable, input_hdf5, fft_hdf5, output_dir, *, endpoints=(0, 1301)):
    """Only the established 60 s VOLTAGE fixture; never callable with APE HDF5.

    An absolute, exact fixture path and a single 1301-bin method-2 job prevent
    this test adapter from becoming an alternative production configuration.
    """
    inp = Path(input_hdf5)
    if str(inp) != "/tmp/ape_validation_noise_voltage.h5":
        raise ConstantQError("Synthetic fixture requires exact known /tmp voltage input")
    exe, fft, out = map(Path, (executable, fft_hdf5, output_dir))
    if not all(p.is_absolute() for p in (inp, exe, fft, out)):
        raise ConstantQError("Synthetic paths must be absolute")
    if not str(out).startswith("/tmp/ape-controller-v1-"):
        raise ConstantQError("Refusing non-isolated synthetic output location")
    if out == fft or out == inp:
        raise ConstantQError("Unsafe synthetic output path")
    if endpoints not in ((0, 1301), (0, 400, 850, 1301)):
        raise ConstantQError("Only the fixed synthetic partitions are authorized")
    args = dict(exe=exe, inp=inp, dataset="voltage", fft=fft,
                fs=2048.0, t_stop=60.0, f_min=5.0, f_max=900.0,
                grid_bins=1301)
    fft_command = _c_argv(output=out/"fft_stage.lpsd", method=3,
                          start=0, count=1301, **args)
    groups = [_c_argv(output=out/f"noise_constq.{i}.lpsd", method=2,
                      start=start, count=end-start, **args)
              for i, (start, end) in enumerate(zip(endpoints, endpoints[1:]))]
    return {"fft": fft_command, "batch": _c_argv(
                output=out/"noise_constq.lpsd", method=2,
                start=0, count=1301, **args),
            "batches": groups, "endpoints": endpoints}


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
    args = dict(exe=exe, inp=inp, dataset="/"+dataset, fft=fft,
                fs=p.fs, t_stop=p.t_stop, f_min=p.f_min, f_max=p.f_max,
                grid_bins=p.grid_bins, overlap=p.overlap, regression=p.regression,
                window=p.window, psll=p.psll_db, scale=p.scale,
                epsilon=p.epsilon, memory_power=p.memory_power,
                kernel_threshold=p.kernel_threshold, a_option=p.a_option, debug=True)
    fft_command = _c_argv(output=out/"fft_stage.lpsd.0", method=3,
                          start=0, count=p.grid_bins, **args)
    constq = [(i, _c_argv(output=out/f"{dataset}.lpsd.{i}", method=2,
                          start=start, count=n, **args))
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
