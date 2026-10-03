"""APE reduced-block quality cleaning adapted from the validated earlier workflow.

Cavity stage:
- reduce cavity transmission into explicit non-overlapping blocks;
- use block average to build a centered physical-time local q90 reference;
- test block minimum against 90% of that reference;
- refine the reference three times using only provisionally good blocks.

Eta0 stage:
- compute one reduced-record eta0 from block-averaged X2f, Y2f, and cavity
  transmission using the approved calibration equation;
- use cavity-good eta0 records to build the same centered local q90 reference;
- test eta0 against 90% of that reference;
- refine three times using only provisionally eta0-good, cavity-good records.

Masks are expanded back to raw-sample indices. No interpolation, smoothing,
gap filling, sample deletion, or X/Y science-quadrature combination occurs.
"""
from __future__ import annotations

from bisect import bisect_left, insort
from dataclasses import dataclass
from typing import Callable, Iterator

import numpy as np

from .ingestion.config import CalibrationConfig, CleaningConfig
from .ingestion.models import APEDXDChunk


class CleaningError(ValueError):
    pass


@dataclass(frozen=True)
class APECleaningResult:
    cavity_good: np.ndarray
    eta0_preliminary: np.ndarray
    eta_good: np.ndarray
    final_good: np.ndarray
    cavity_reference: np.ndarray
    eta0_reference: np.ndarray
    metadata: dict[str, object]


@dataclass(frozen=True)
class APECleaningReferencePlan:
    """Full-record reduced-block quality plan."""

    sample_rate_hz: float
    total_samples: int
    block_samples: int
    block_start_indices: np.ndarray
    block_stop_indices: np.ndarray
    block_midpoint_s: np.ndarray
    cavity_average_blocks: np.ndarray
    cavity_minimum_blocks: np.ndarray
    cavity_reference_blocks: np.ndarray
    cavity_good_blocks: np.ndarray
    eta0_blocks: np.ndarray
    eta0_reference_blocks: np.ndarray
    eta0_good_blocks: np.ndarray
    final_good_blocks: np.ndarray
    metadata: dict[str, object]


def _validate_cleaning(cleaning: CleaningConfig) -> None:
    if not np.isfinite(cleaning.block_size_seconds) or cleaning.block_size_seconds <= 0:
        raise CleaningError("block_size_seconds must be positive and finite")
    if (
        not np.isfinite(cleaning.reference_window_seconds)
        or cleaning.reference_window_seconds <= 0
    ):
        raise CleaningError("reference_window_seconds must be positive and finite")
    if not (0.0 < cleaning.reference_quantile <= 1.0):
        raise CleaningError("reference_quantile must be in (0, 1]")
    if cleaning.refinement_iterations < 0:
        raise CleaningError("refinement_iterations must be non-negative")
    if not (0.0 < cleaning.cavity_lower_fraction <= 1.0):
        raise CleaningError("cavity_lower_fraction must be in (0, 1]")
    if not (0.0 < cleaning.eta0_lower_fraction <= 1.0):
        raise CleaningError("eta0_lower_fraction must be in (0, 1]")


def _samples_per_block(sample_rate_hz: float, block_size_seconds: float) -> int:
    value = int(round(sample_rate_hz * block_size_seconds))
    if value <= 0:
        raise CleaningError("cleaning block must contain at least one sample")
    actual = value / sample_rate_hz
    if not np.isclose(actual, block_size_seconds, rtol=0.0, atol=1e-12):
        raise CleaningError(
            "block_size_seconds must map to an integer number of raw samples"
        )
    return value


def _positive_finite_calibration(calibration: CalibrationConfig) -> tuple[float, float, float, float]:
    values = (
        float(calibration.cavity_transmission_w_per_v),
        float(calibration.responsivity_a_per_w),
        float(calibration.transimpedance_gain_v_per_a),
        float(calibration.second_harmonic_pd_response),
    )
    names = (
        "cavity_transmission_w_per_v",
        "responsivity_a_per_w",
        "transimpedance_gain_v_per_a",
        "second_harmonic_pd_response",
    )
    for name, value in zip(names, values):
        if not np.isfinite(value) or value <= 0:
            raise CleaningError(f"{name} must be positive and finite")
    return values


def _linear_quantile_from_sorted(values: list[float], probability: float) -> float:
    if not values:
        raise CleaningError("cannot calculate a quantile from an empty window")
    position = (len(values) - 1) * probability
    lower = int(np.floor(position))
    upper = int(np.ceil(position))
    if lower == upper:
        return float(values[lower])
    fraction = position - lower
    return float(values[lower] + fraction * (values[upper] - values[lower]))


def _rolling_quantile_reference(
    values: np.ndarray,
    timestamps_s: np.ndarray,
    eligible: np.ndarray,
    *,
    window_seconds: float,
    quantile: float,
) -> np.ndarray:
    """Exact centered physical-time rolling quantile with no interpolation/fallback.

    Window membership follows the earlier workflow's scientific predicate:
    abs(t_j - t_i) <= window_seconds / 2.
    """
    values = np.asarray(values, dtype=np.float64)
    timestamps_s = np.asarray(timestamps_s, dtype=np.float64)
    eligible = np.asarray(eligible, dtype=bool)
    if values.ndim != 1 or timestamps_s.shape != values.shape or eligible.shape != values.shape:
        raise CleaningError("rolling-reference arrays must be aligned one-dimensional arrays")
    if values.size and np.any(np.diff(timestamps_s) <= 0):
        raise CleaningError("block timestamps must be strictly increasing")

    n = values.size
    refs = np.full(n, np.nan, dtype=np.float64)
    half_window = window_seconds / 2.0
    ordered: list[float] = []
    left = 0
    right = 0

    for i in range(n):
        target = timestamps_s[i]
        while right < n and abs(timestamps_s[right] - target) <= half_window:
            if eligible[right] and np.isfinite(values[right]):
                insort(ordered, float(values[right]))
            right += 1
        while left < right and abs(timestamps_s[left] - target) > half_window:
            if eligible[left] and np.isfinite(values[left]):
                value = float(values[left])
                pos = bisect_left(ordered, value)
                if pos >= len(ordered) or ordered[pos] != value:
                    raise CleaningError("rolling quantile window bookkeeping failed")
                ordered.pop(pos)
            left += 1
        if ordered:
            refs[i] = _linear_quantile_from_sorted(ordered, quantile)
    return refs


def _refined_mask(
    reference_values: np.ndarray,
    decision_values: np.ndarray,
    timestamps_s: np.ndarray,
    base_eligible: np.ndarray,
    *,
    window_seconds: float,
    quantile: float,
    lower_fraction: float,
    refinement_iterations: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce the earlier reference-refinement pattern for one quality variable."""
    reference_values = np.asarray(reference_values, dtype=np.float64)
    decision_values = np.asarray(decision_values, dtype=np.float64)
    base_eligible = np.asarray(base_eligible, dtype=bool)

    eligible = base_eligible & np.isfinite(reference_values)
    refs = _rolling_quantile_reference(
        reference_values,
        timestamps_s,
        eligible,
        window_seconds=window_seconds,
        quantile=quantile,
    )

    for _ in range(refinement_iterations):
        provisional = (
            base_eligible
            & np.isfinite(decision_values)
            & np.isfinite(refs)
            & (decision_values >= lower_fraction * refs)
        )
        refs = _rolling_quantile_reference(
            reference_values,
            timestamps_s,
            provisional,
            window_seconds=window_seconds,
            quantile=quantile,
        )

    good = (
        base_eligible
        & np.isfinite(decision_values)
        & np.isfinite(refs)
        & (decision_values >= lower_fraction * refs)
    )
    return refs, good


def _eta0_from_reduced_averages(
    x2_average: np.ndarray,
    y2_average: np.ndarray,
    cavity_average: np.ndarray,
    calibration: CalibrationConfig,
) -> np.ndarray:
    c_trans, q, gain, h2 = _positive_finite_calibration(calibration)
    i0 = c_trans * np.asarray(cavity_average, dtype=np.float64)
    r2f = np.hypot(
        np.asarray(x2_average, dtype=np.float64),
        np.asarray(y2_average, dtype=np.float64),
    )
    eta0 = np.full(i0.shape, np.nan, dtype=np.float64)
    valid = (
        np.isfinite(i0)
        & np.isfinite(r2f)
        & (i0 > 0.0)
        & (r2f > 0.0)
    )
    arg = np.full(i0.shape, np.nan, dtype=np.float64)
    arg[valid] = 2.0 * np.sqrt(2.0) * r2f[valid] / (h2 * q * gain * i0[valid])
    valid &= np.isfinite(arg) & (arg > 0.0)
    eta0[valid] = np.sqrt(arg[valid])
    return eta0


def _reduce_chunk_blocks(
    chunk: APEDXDChunk,
    block_samples: int,
) -> tuple[np.ndarray, ...]:
    if chunk.start_index % block_samples != 0:
        raise CleaningError("chunk starts must align to reduced cleaning blocks")
    n_blocks = (chunk.n_samples + block_samples - 1) // block_samples

    starts = chunk.start_index + np.arange(n_blocks, dtype=np.int64) * block_samples
    stops = np.minimum(starts + block_samples, chunk.stop_index)
    # Physical midpoint of the half-open block interval [start/fs, stop/fs).
    midpoints = (starts + stops) / (2.0 * chunk.sample_rate_hz)

    cavity_avg = np.full(n_blocks, np.nan, dtype=np.float64)
    cavity_min = np.full(n_blocks, np.nan, dtype=np.float64)
    x2_avg = np.full(n_blocks, np.nan, dtype=np.float64)
    y2_avg = np.full(n_blocks, np.nan, dtype=np.float64)

    v = np.asarray(chunk.cavity_transmission, dtype=np.float64)
    x2 = np.asarray(chunk.x_2f, dtype=np.float64)
    y2 = np.asarray(chunk.y_2f, dtype=np.float64)

    for b in range(n_blocks):
        lo = b * block_samples
        hi = min(chunk.n_samples, lo + block_samples)
        vb = v[lo:hi]
        xb = x2[lo:hi]
        yb = y2[lo:hi]
        if vb.size and np.all(np.isfinite(vb)):
            cavity_avg[b] = float(np.mean(vb))
            cavity_min[b] = float(np.min(vb))
        if xb.size and np.all(np.isfinite(xb)):
            x2_avg[b] = float(np.mean(xb))
        if yb.size and np.all(np.isfinite(yb)):
            y2_avg[b] = float(np.mean(yb))

    return starts, stops, midpoints, cavity_avg, cavity_min, x2_avg, y2_avg


def build_cleaning_reference_plan(
    chunk_factory: Callable[[], Iterator[APEDXDChunk]],
    cleaning: CleaningConfig,
    calibration: CalibrationConfig,
) -> APECleaningReferencePlan:
    """Build the complete reduced-block quality plan in one raw-data pass."""
    _validate_cleaning(cleaning)
    _positive_finite_calibration(calibration)

    starts_parts: list[np.ndarray] = []
    stops_parts: list[np.ndarray] = []
    time_parts: list[np.ndarray] = []
    cavity_avg_parts: list[np.ndarray] = []
    cavity_min_parts: list[np.ndarray] = []
    x2_avg_parts: list[np.ndarray] = []
    y2_avg_parts: list[np.ndarray] = []

    fs: float | None = None
    block_samples: int | None = None
    expected_start = 0
    saw_chunk = False

    for chunk in chunk_factory():
        saw_chunk = True
        if fs is None:
            fs = chunk.sample_rate_hz
            block_samples = _samples_per_block(fs, cleaning.block_size_seconds)
        if chunk.sample_rate_hz != fs:
            raise CleaningError("sample rate changes across cleaning chunks")
        if chunk.start_index != expected_start:
            raise CleaningError("cleaning chunks must be contiguous and ordered")
        if chunk.start_index % block_samples != 0:
            raise CleaningError("cleaning chunk boundary is not block-aligned")

        reduced = _reduce_chunk_blocks(chunk, block_samples)
        starts_parts.append(reduced[0])
        stops_parts.append(reduced[1])
        time_parts.append(reduced[2])
        cavity_avg_parts.append(reduced[3])
        cavity_min_parts.append(reduced[4])
        x2_avg_parts.append(reduced[5])
        y2_avg_parts.append(reduced[6])
        expected_start = chunk.stop_index

    if not saw_chunk or fs is None or block_samples is None:
        raise CleaningError("no chunks supplied by chunk_factory")

    starts = np.concatenate(starts_parts)
    stops = np.concatenate(stops_parts)
    timestamps = np.concatenate(time_parts)
    cavity_avg = np.concatenate(cavity_avg_parts)
    cavity_min = np.concatenate(cavity_min_parts)
    x2_avg = np.concatenate(x2_avg_parts)
    y2_avg = np.concatenate(y2_avg_parts)

    # Refuse overlapping or duplicated reduced blocks. A partial block may exist
    # only at the end of the complete record.
    if starts.size:
        if np.any(starts[1:] < stops[:-1]):
            raise CleaningError("reduced cleaning blocks overlap")
        if np.any(starts[1:] != stops[:-1]):
            raise CleaningError("reduced cleaning blocks contain a gap")

    cavity_base = np.isfinite(cavity_avg)
    cavity_reference, cavity_good = _refined_mask(
        cavity_avg,
        cavity_min,
        timestamps,
        cavity_base,
        window_seconds=cleaning.reference_window_seconds,
        quantile=cleaning.reference_quantile,
        lower_fraction=cleaning.cavity_lower_fraction,
        refinement_iterations=cleaning.refinement_iterations,
    )

    eta0 = _eta0_from_reduced_averages(x2_avg, y2_avg, cavity_avg, calibration)
    eta_base = cavity_good & np.isfinite(eta0)
    eta_reference, eta_good = _refined_mask(
        eta0,
        eta0,
        timestamps,
        eta_base,
        window_seconds=cleaning.reference_window_seconds,
        quantile=cleaning.reference_quantile,
        lower_fraction=cleaning.eta0_lower_fraction,
        refinement_iterations=cleaning.refinement_iterations,
    )
    final_good = cavity_good & eta_good

    return APECleaningReferencePlan(
        sample_rate_hz=fs,
        total_samples=expected_start,
        block_samples=block_samples,
        block_start_indices=starts,
        block_stop_indices=stops,
        block_midpoint_s=timestamps,
        cavity_average_blocks=cavity_avg,
        cavity_minimum_blocks=cavity_min,
        cavity_reference_blocks=cavity_reference,
        cavity_good_blocks=cavity_good,
        eta0_blocks=eta0,
        eta0_reference_blocks=eta_reference,
        eta0_good_blocks=eta_good,
        final_good_blocks=final_good,
        metadata={
            "method": "ape_reduced_q90_refined_cleaning_v2",
            "source_inspiration": "validated earlier APE cavity_lock / pem_modulation workflow",
            "block_size_seconds": float(cleaning.block_size_seconds),
            "reference_window_seconds": float(cleaning.reference_window_seconds),
            "reference_quantile": float(cleaning.reference_quantile),
            "refinement_iterations": int(cleaning.refinement_iterations),
            "cavity_lower_fraction": float(cleaning.cavity_lower_fraction),
            "eta0_lower_fraction": float(cleaning.eta0_lower_fraction),
            "cavity_reference_statistic": "reduced-block average",
            "cavity_decision_statistic": "reduced-block minimum",
            "eta0_reference_statistic": "reduced-record eta0",
            "eta0_decision_statistic": "reduced-record eta0",
            "window_definition": "abs(t_j - t_i) <= reference_window_seconds / 2",
            "quantile_method": "linear interpolation at (n - 1) * q",
            "reference_interpolation": "none",
            "reference_fallback": "none",
            "smoothing": "none",
            "gap_filling": "none",
            "indices_preserved": True,
            "rho_not_computed_during_cleaning": True,
            "raw_record_retained_in_memory": False,
        },
    )


def _blocks_for_chunk(
    plan: APECleaningReferencePlan,
    chunk: APEDXDChunk,
) -> np.ndarray:
    first = int(np.searchsorted(plan.block_stop_indices, chunk.start_index, side="right"))
    last = int(np.searchsorted(plan.block_start_indices, chunk.stop_index, side="left"))
    indices = np.arange(first, last, dtype=np.int64)
    if indices.size == 0:
        raise CleaningError("cleaning plan does not cover this chunk")
    if plan.block_start_indices[indices[0]] != chunk.start_index:
        raise CleaningError("chunk start is not aligned to cleaning plan")
    if plan.block_stop_indices[indices[-1]] != chunk.stop_index:
        raise CleaningError("chunk stop is not aligned to cleaning plan")
    return indices


def _expand_blocks_for_chunk(
    block_values: np.ndarray,
    plan: APECleaningReferencePlan,
    chunk: APEDXDChunk,
    indices: np.ndarray,
    *,
    dtype=None,
) -> np.ndarray:
    pieces = []
    for idx in indices:
        count = int(plan.block_stop_indices[idx] - plan.block_start_indices[idx])
        pieces.append(np.full(count, block_values[idx], dtype=dtype))
    out = np.concatenate(pieces)
    if out.shape != (chunk.n_samples,):
        raise CleaningError("expanded cleaning plan does not match chunk length")
    return out


def clean_chunk_with_plan(
    chunk: APEDXDChunk,
    plan: APECleaningReferencePlan,
    cleaning: CleaningConfig,
    calibration: CalibrationConfig,
) -> APECleaningResult:
    """Expand the precomputed reduced-block quality decisions onto one raw chunk."""
    _validate_cleaning(cleaning)
    _positive_finite_calibration(calibration)
    if chunk.sample_rate_hz != plan.sample_rate_hz:
        raise CleaningError("chunk sample rate differs from cleaning plan")

    indices = _blocks_for_chunk(plan, chunk)
    cavity_good = _expand_blocks_for_chunk(
        plan.cavity_good_blocks, plan, chunk, indices, dtype=bool
    )
    eta_good = _expand_blocks_for_chunk(
        plan.eta0_good_blocks, plan, chunk, indices, dtype=bool
    )
    final_good = _expand_blocks_for_chunk(
        plan.final_good_blocks, plan, chunk, indices, dtype=bool
    )
    eta0 = _expand_blocks_for_chunk(
        plan.eta0_blocks, plan, chunk, indices, dtype=np.float64
    )
    cavity_reference = _expand_blocks_for_chunk(
        plan.cavity_reference_blocks, plan, chunk, indices, dtype=np.float64
    )
    eta0_reference = _expand_blocks_for_chunk(
        plan.eta0_reference_blocks, plan, chunk, indices, dtype=np.float64
    )

    eta0[~cavity_good] = np.nan
    eta0_reference[~cavity_good] = np.nan

    return APECleaningResult(
        cavity_good=cavity_good,
        eta0_preliminary=eta0,
        eta_good=eta_good,
        final_good=final_good,
        cavity_reference=cavity_reference,
        eta0_reference=eta0_reference,
        metadata=dict(plan.metadata),
    )


def clean_chunk(
    chunk: APEDXDChunk,
    cleaning: CleaningConfig,
    calibration: CalibrationConfig,
) -> APECleaningResult:
    """Convenience single-chunk cleaning path used by controlled tests."""
    plan = build_cleaning_reference_plan(
        lambda: iter((chunk,)),
        cleaning,
        calibration,
    )
    return clean_chunk_with_plan(chunk, plan, cleaning, calibration)
