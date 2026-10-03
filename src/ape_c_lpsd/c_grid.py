"""Exact no-data planner for the unchanged C LPSD method-0 grid."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any
import uuid

import numpy as np


class CGridError(RuntimeError):
    pass


@dataclass(frozen=True)
class CNativeGridPlan:
    frequency_hz: np.ndarray
    required_dft_samples: np.ndarray
    fmin_hz: float
    fmax_hz: float
    jdes: int
    relative_step: float
    segment_samples: int
    sample_rate_hz: float


def _round_positive_c(value: float) -> int:
    """C round() for the positive values used by lpsd.c."""
    return math.floor(value + 0.5)


def c_frequency(index: int, fmin_hz: float, fmax_hz: float, jdes: int) -> float:
    """Mirror lpsd.c calc_params: fmin*exp(i*g/(Jdes-1))."""
    g = math.log(fmax_hz / fmin_hz)
    return fmin_hz * math.exp(index * g / (jdes - 1.0))


def c_required_dft_samples(frequency_hz: float, relative_step: float, sample_rate_hz: float) -> int:
    """Mirror method-0 calc_params: round(fs/(f*expm1(g/(J-1))))."""
    return _round_positive_c(sample_rate_hz / (frequency_hz * relative_step))


def _first_required(fmin_hz: float, fmax_hz: float, jdes: int, sample_rate_hz: float) -> int:
    step = math.exp(math.log(fmax_hz / fmin_hz) / (jdes - 1.0)) - 1.0
    return c_required_dft_samples(fmin_hz, step, sample_rate_hz)


def densest_feasible_jdes(*, fmin_hz: float, fmax_hz: float, segment_samples: int, sample_rate_hz: float) -> int:
    """Largest C Jdes for which its first (longest) DFT fits the segment."""
    if not (0.0 < fmin_hz < fmax_hz < sample_rate_hz / 2.0):
        raise CGridError('frequency bounds must satisfy 0 < fmin < fmax < Nyquist')
    if segment_samples < 2 or sample_rate_hz <= 0.0:
        raise CGridError('segment_samples and sample_rate_hz must be positive')
    estimate = max(2, math.floor(math.log(fmax_hz / fmin_hz) / math.log1p(sample_rate_hz / (segment_samples * fmin_hz))) + 1)
    while _first_required(fmin_hz, fmax_hz, estimate + 1, sample_rate_hz) <= segment_samples:
        estimate += 1
    while _first_required(fmin_hz, fmax_hz, estimate, sample_rate_hz) > segment_samples:
        estimate -= 1
    return estimate


def plan_c_native_grid(*, fmin_hz: float, fmax_hz: float, segment_samples: int, sample_rate_hz: float) -> CNativeGridPlan:
    jdes = densest_feasible_jdes(fmin_hz=fmin_hz, fmax_hz=fmax_hz, segment_samples=segment_samples, sample_rate_hz=sample_rate_hz)
    g = math.log(fmax_hz / fmin_hz)
    step = math.exp(g / (jdes - 1.0)) - 1.0
    # Python math.exp and C exp both call the platform libm; retain scalar order to mirror lpsd.c.
    freq = np.fromiter((fmin_hz * math.exp(i * g / (jdes - 1.0)) for i in range(jdes)), dtype=np.float64, count=jdes)
    required = np.fromiter((c_required_dft_samples(float(x), step, sample_rate_hz) for x in freq), dtype=np.int64, count=jdes)
    if not np.all(np.isfinite(freq)) or not np.all(np.diff(freq) > 0.0):
        raise CGridError('generated C frequency grid is non-finite or unordered')
    bad = np.flatnonzero(required > segment_samples)
    if bad.size:
        i = int(bad[0])
        raise CGridError(f'C-native grid unexpectedly infeasible at node {i}: required={required[i]}, available={segment_samples}')
    if _first_required(fmin_hz, fmax_hz, jdes + 1, sample_rate_hz) <= segment_samples:
        raise CGridError('planner did not select the densest feasible Jdes')
    return CNativeGridPlan(freq, required, fmin_hz, fmax_hz, jdes, step, segment_samples, sample_rate_hz)


def write_manifest(*, plan: CNativeGridPlan, output_directory: Path, c_options: dict[str, Any], reference_grid: Path | None = None) -> Path:
    """Write one immutable NPZ/JSON pair; never overwrite a prior plan."""
    output_directory.mkdir(parents=True, exist_ok=True)
    stamp = uuid.uuid4().hex
    base = output_directory / f'c-native-grid-{stamp}'
    npz = base.with_suffix('.npz')
    if npz.exists():
        raise FileExistsError(f'refusing to overwrite grid manifest: {npz}')
    np.savez_compressed(npz, frequency_hz=plan.frequency_hz, required_dft_samples=plan.required_dft_samples)
    digest = hashlib.sha256(npz.read_bytes()).hexdigest()
    comparison: dict[str, Any] = {'reference_grid_path': None, 'reference_grid_sha256': None, 'reference_grid_node_count': None}
    if reference_grid is not None:
        comparison.update(reference_grid_path=str(reference_grid.resolve()), reference_grid_sha256=hashlib.sha256(reference_grid.read_bytes()).hexdigest())
        with np.load(reference_grid, allow_pickle=False) as ref:
            comparison['reference_grid_node_count'] = int(len(ref['frequency_hz']))
            comparison['reference_grid_first_hz'] = float(ref['frequency_hz'][0])
            comparison['reference_grid_last_hz'] = float(ref['frequency_hz'][-1])
    meta = {
        'kind': 'C-native method-0 grid manifest; no signal data or C execution',
        'frequency_formula_source': 'lpsd.c calc_params: f_j=fmin*exp(j*log(fmax/fmin)/(Jdes-1))',
        'required_length_formula_source': 'lpsd.c calc_params: round(fs/(f_j*(exp(g/(Jdes-1))-1)))',
        'segment_samples': plan.segment_samples, 'sample_rate_hz': plan.sample_rate_hz,
        'fmin_hz': plan.fmin_hz, 'fmax_hz': plan.fmax_hz, 'jdes': plan.jdes,
        'relative_step': plan.relative_step, 'first_frequency_hz': float(plan.frequency_hz[0]),
        'last_frequency_hz': float(plan.frequency_hz[-1]), 'maximum_required_dft_samples': int(plan.required_dft_samples.max()),
        'maximum_required_dft_index': int(plan.required_dft_samples.argmax()),
        'nodes_exceeding_segment': int(np.count_nonzero(plan.required_dft_samples > plan.segment_samples)),
        'manifest_npz_sha256': digest, 'c_command_options': c_options, 'reference_comparison': comparison,
        'batch_note': 'Exact native nodes can be calculated in one C invocation with -n Jdes -N 0. Fixed -n/-N batching is exact only when n divides Jdes because lpsd.c starts at nspec*iter.',
    }
    json_path = base.with_suffix('.json')
    if json_path.exists():
        raise FileExistsError(f'refusing to overwrite grid manifest: {json_path}')
    json_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + '\n')
    return json_path
