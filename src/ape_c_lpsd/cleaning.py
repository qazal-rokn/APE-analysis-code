"""APE reduced-block quality cleaning adapted from the validated earlier workflow.

Cavity stage (unchanged):
- reduce cavity transmission into explicit non-overlapping blocks;
- use block average to build a centered physical-time local q90 reference;
- test block minimum against 90% of that reference;
- refine the reference three times using only provisionally good blocks.

Eta0 stage:
- calibrate eta0 at the original raw-sample rate from X2f, Y2f, and cavity
  transmission using the approved calibration equation;
- reduce sample-level eta0 into a block mean and block minimum;
- build a past-only 60 s q90 reference from previously accepted cavity-good
  eta0 block means;
- test the eta0 block minimum against the configured lower fraction;
- rejected eta0 blocks never update the future reference;
- if the accepted-history window empties inside a cavity-good segment, hold the last healthy reference;\n- reset eta0 reference state at the start of each cavity-good segment.

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
class APECavityReferencePlan:
    """Full-record cavity-only reduced-block quality plan."""

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
    metadata: dict[str, object]


@dataclass(frozen=True)
class APEEtaSegmentPlan:
    """Eta0 quality plan for one cavity-good candidate segment."""

    sample_rate_hz: float
    start_index: int
    stop_index: int
    block_samples: int
    block_start_indices: np.ndarray
    block_stop_indices: np.ndarray
    block_midpoint_s: np.ndarray
    eta0_average_blocks: np.ndarray
    eta0_minimum_blocks: np.ndarray
    eta0_reference_blocks: np.ndarray
    eta0_good_blocks: np.ndarray
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
    eta0_minimum_blocks: np.ndarray
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


def _protected_past_quantile_mask(
    reference_values: np.ndarray,
    decision_values: np.ndarray,
    timestamps_s: np.ndarray,
    base_eligible: np.ndarray,
    *,
    window_seconds: float,
    quantile: float,
    lower_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Past-only q90 reference built only from previously accepted records.

    The current record is never included in its own reference. Rejected records
    never update the future reference. If a sustained failure (or cavity-bad gap)
    leaves no accepted records inside the past window, the last healthy reference
    is held until a record passes again.
    """
    reference_values = np.asarray(reference_values, dtype=np.float64)
    decision_values = np.asarray(decision_values, dtype=np.float64)
    timestamps_s = np.asarray(timestamps_s, dtype=np.float64)
    base_eligible = np.asarray(base_eligible, dtype=bool)

    if (
        reference_values.ndim != 1
        or decision_values.shape != reference_values.shape
        or timestamps_s.shape != reference_values.shape
        or base_eligible.shape != reference_values.shape
    ):
        raise CleaningError("protected-reference arrays must be aligned one-dimensional arrays")
    if reference_values.size and np.any(np.diff(timestamps_s) <= 0):
        raise CleaningError("block timestamps must be strictly increasing")

    refs = np.full(reference_values.shape, np.nan, dtype=np.float64)
    good = np.zeros(reference_values.shape, dtype=bool)

    history: list[tuple[float, float]] = []
    ordered: list[float] = []
    left = 0
    last_valid_reference = np.nan

    for k, target in enumerate(timestamps_s):
        while left < len(history) and history[left][0] < target - window_seconds:
            old_value = history[left][1]
            pos = bisect_left(ordered, old_value)
            if pos >= len(ordered) or ordered[pos] != old_value:
                raise CleaningError("protected quantile window bookkeeping failed")
            ordered.pop(pos)
            left += 1

        if ordered:
            ref = _linear_quantile_from_sorted(ordered, quantile)
            last_valid_reference = ref
        elif np.isfinite(last_valid_reference):
            ref = float(last_valid_reference)
        elif (
            base_eligible[k]
            and np.isfinite(reference_values[k])
            and np.isfinite(decision_values[k])
        ):
            # Initial bootstrap: the acquisition is assumed to begin in its
            # normal operating state, matching the validated diagnostic.
            ref = float(reference_values[k])
            last_valid_reference = ref
        else:
            continue

        refs[k] = ref

        if (
            base_eligible[k]
            and np.isfinite(reference_values[k])
            and np.isfinite(decision_values[k])
            and decision_values[k] >= lower_fraction * ref
        ):
            good[k] = True
            value = float(reference_values[k])
            history.append((float(target), value))
            insort(ordered, value)

        # Periodically compact expired history so long records do not retain
        # an ever-growing Python list of old timestamps.
        if left > 4096 and left > len(history) // 2:
            history = history[left:]
            left = 0

    return refs, good



def _protected_past_quantile_mask_by_segments(
    reference_values: np.ndarray,
    decision_values: np.ndarray,
    timestamps_s: np.ndarray,
    base_eligible: np.ndarray,
    segment_mask: np.ndarray,
    *,
    window_seconds: float,
    quantile: float,
    lower_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the protected eta0 reference independently in each cavity-good run.

    This implements the researcher-approved efficient-search rule: eta0 state is
    reset at the beginning of every cavity-good candidate segment. Within a
    segment, rejected/invalid eta0 blocks do not update the reference and the
    last healthy reference is held.
    """
    reference_values = np.asarray(reference_values, dtype=np.float64)
    decision_values = np.asarray(decision_values, dtype=np.float64)
    timestamps_s = np.asarray(timestamps_s, dtype=np.float64)
    base_eligible = np.asarray(base_eligible, dtype=bool)
    segment_mask = np.asarray(segment_mask, dtype=bool)

    if (
        reference_values.ndim != 1
        or decision_values.shape != reference_values.shape
        or timestamps_s.shape != reference_values.shape
        or base_eligible.shape != reference_values.shape
        or segment_mask.shape != reference_values.shape
    ):
        raise CleaningError("segmented protected-reference arrays must be aligned")

    refs = np.full(reference_values.shape, np.nan, dtype=np.float64)
    good = np.zeros(reference_values.shape, dtype=bool)

    idx = np.flatnonzero(segment_mask)
    if idx.size == 0:
        return refs, good

    split = np.flatnonzero(np.diff(idx) > 1)
    run_first = np.concatenate(([0], split + 1))
    run_last = np.concatenate((split, [idx.size - 1]))

    for first_pos, last_pos in zip(run_first, run_last):
        lo = int(idx[first_pos])
        hi = int(idx[last_pos]) + 1
        sub_refs, sub_good = _protected_past_quantile_mask(
            reference_values[lo:hi],
            decision_values[lo:hi],
            timestamps_s[lo:hi],
            base_eligible[lo:hi],
            window_seconds=window_seconds,
            quantile=quantile,
            lower_fraction=lower_fraction,
        )
        refs[lo:hi] = sub_refs
        good[lo:hi] = sub_good

    return refs, good


def _eta0_from_samples(
    x2: np.ndarray,
    y2: np.ndarray,
    cavity_transmission: np.ndarray,
    calibration: CalibrationConfig,
) -> np.ndarray:
    c_trans, q, gain, h2 = _positive_finite_calibration(calibration)
    i0 = c_trans * np.asarray(cavity_transmission, dtype=np.float64)
    r2f = np.hypot(
        np.asarray(x2, dtype=np.float64),
        np.asarray(y2, dtype=np.float64),
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
    calibration: CalibrationConfig,
) -> tuple[np.ndarray, ...]:
    """Reduce one aligned raw chunk without Python loops over full blocks.

    This is numerically equivalent to the earlier per-block implementation:
    a block statistic is valid only when every sample in that block is finite.
    The final partial block, when present, is handled separately.
    """
    if chunk.start_index % block_samples != 0:
        raise CleaningError("chunk starts must align to reduced cleaning blocks")
    n_blocks = (chunk.n_samples + block_samples - 1) // block_samples

    starts = chunk.start_index + np.arange(n_blocks, dtype=np.int64) * block_samples
    stops = np.minimum(starts + block_samples, chunk.stop_index)
    # Physical midpoint of the half-open block interval [start/fs, stop/fs).
    midpoints = (starts + stops) / (2.0 * chunk.sample_rate_hz)

    cavity_avg = np.full(n_blocks, np.nan, dtype=np.float64)
    cavity_min = np.full(n_blocks, np.nan, dtype=np.float64)
    eta0_avg = np.full(n_blocks, np.nan, dtype=np.float64)
    eta0_min = np.full(n_blocks, np.nan, dtype=np.float64)

    v = np.asarray(chunk.cavity_transmission, dtype=np.float64)
    x2 = np.asarray(chunk.x_2f, dtype=np.float64)
    y2 = np.asarray(chunk.y_2f, dtype=np.float64)
    eta0_samples = _eta0_from_samples(x2, y2, v, calibration)

    n_full = chunk.n_samples // block_samples
    full_stop = n_full * block_samples

    if n_full:
        vb = v[:full_stop].reshape(n_full, block_samples)
        eb = eta0_samples[:full_stop].reshape(n_full, block_samples)

        valid_v = np.all(np.isfinite(vb), axis=1)
        valid_eta = np.all(np.isfinite(eb), axis=1)

        v_mean = np.mean(vb, axis=1)
        v_min = np.min(vb, axis=1)
        eta_mean = np.mean(eb, axis=1)
        eta_minimum = np.min(eb, axis=1)

        cavity_avg[:n_full] = np.where(valid_v, v_mean, np.nan)
        cavity_min[:n_full] = np.where(valid_v, v_min, np.nan)
        eta0_avg[:n_full] = np.where(valid_eta, eta_mean, np.nan)
        eta0_min[:n_full] = np.where(valid_eta, eta_minimum, np.nan)

    if full_stop < chunk.n_samples:
        vb = v[full_stop:]
        eb = eta0_samples[full_stop:]
        if vb.size and np.all(np.isfinite(vb)):
            cavity_avg[-1] = float(np.mean(vb))
            cavity_min[-1] = float(np.min(vb))
        if eb.size and np.all(np.isfinite(eb)):
            eta0_avg[-1] = float(np.mean(eb))
            eta0_min[-1] = float(np.min(eb))

    return starts, stops, midpoints, cavity_avg, cavity_min, eta0_avg, eta0_min


def _reduce_chunk_cavity_blocks(
    chunk: APEDXDChunk,
    block_samples: int,
) -> tuple[np.ndarray, ...]:
    """Vectorized cavity-only block reduction; no eta0 calibration."""
    if chunk.start_index % block_samples != 0:
        raise CleaningError("chunk starts must align to reduced cleaning blocks")
    n_blocks = (chunk.n_samples + block_samples - 1) // block_samples

    starts = chunk.start_index + np.arange(n_blocks, dtype=np.int64) * block_samples
    stops = np.minimum(starts + block_samples, chunk.stop_index)
    midpoints = (starts + stops) / (2.0 * chunk.sample_rate_hz)

    cavity_avg = np.full(n_blocks, np.nan, dtype=np.float64)
    cavity_min = np.full(n_blocks, np.nan, dtype=np.float64)
    v = np.asarray(chunk.cavity_transmission, dtype=np.float64)

    n_full = chunk.n_samples // block_samples
    full_stop = n_full * block_samples
    if n_full:
        vb = v[:full_stop].reshape(n_full, block_samples)
        valid = np.all(np.isfinite(vb), axis=1)
        cavity_avg[:n_full] = np.where(valid, np.mean(vb, axis=1), np.nan)
        cavity_min[:n_full] = np.where(valid, np.min(vb, axis=1), np.nan)

    if full_stop < chunk.n_samples:
        vb = v[full_stop:]
        if vb.size and np.all(np.isfinite(vb)):
            cavity_avg[-1] = float(np.mean(vb))
            cavity_min[-1] = float(np.min(vb))

    return starts, stops, midpoints, cavity_avg, cavity_min


def build_cavity_reference_plan(
    chunk_factory: Callable[[], Iterator[APEDXDChunk]],
    cleaning: CleaningConfig,
) -> APECavityReferencePlan:
    """Build the approved cavity mask without calculating sample-level eta0."""
    _validate_cleaning(cleaning)

    starts_parts: list[np.ndarray] = []
    stops_parts: list[np.ndarray] = []
    time_parts: list[np.ndarray] = []
    cavity_avg_parts: list[np.ndarray] = []
    cavity_min_parts: list[np.ndarray] = []

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

        reduced = _reduce_chunk_cavity_blocks(chunk, block_samples)
        starts_parts.append(reduced[0])
        stops_parts.append(reduced[1])
        time_parts.append(reduced[2])
        cavity_avg_parts.append(reduced[3])
        cavity_min_parts.append(reduced[4])
        expected_start = chunk.stop_index

    if not saw_chunk or fs is None or block_samples is None:
        raise CleaningError("no chunks supplied by chunk_factory")

    starts = np.concatenate(starts_parts)
    stops = np.concatenate(stops_parts)
    timestamps = np.concatenate(time_parts)
    cavity_avg = np.concatenate(cavity_avg_parts)
    cavity_min = np.concatenate(cavity_min_parts)

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

    return APECavityReferencePlan(
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
        metadata={
            "method": "ape_cavity_only_reduced_q90_refined_v1",
            "block_size_seconds": float(cleaning.block_size_seconds),
            "reference_window_seconds": float(cleaning.reference_window_seconds),
            "reference_quantile": float(cleaning.reference_quantile),
            "refinement_iterations": int(cleaning.refinement_iterations),
            "cavity_lower_fraction": float(cleaning.cavity_lower_fraction),
            "cavity_reference_statistic": "reduced-block average",
            "cavity_decision_statistic": "reduced-block minimum",
            "eta0_calculated": False,
        },
    )


def build_eta0_segment_plan(
    chunk_factory: Callable[[], Iterator[APEDXDChunk]],
    cleaning: CleaningConfig,
    calibration: CalibrationConfig,
) -> APEEtaSegmentPlan:
    """Clean eta0 inside one already cavity-good candidate segment.

    The protected eta0 history is intentionally reset at the segment start,
    using the first valid eta0 block mean as its bootstrap reference. This is
    the researcher-approved efficient-search rule.
    """
    _validate_cleaning(cleaning)
    _positive_finite_calibration(calibration)

    starts_parts: list[np.ndarray] = []
    stops_parts: list[np.ndarray] = []
    time_parts: list[np.ndarray] = []
    eta_avg_parts: list[np.ndarray] = []
    eta_min_parts: list[np.ndarray] = []

    fs: float | None = None
    block_samples: int | None = None
    expected_start: int | None = None
    segment_start: int | None = None
    segment_stop: int | None = None

    for chunk in chunk_factory():
        if fs is None:
            fs = chunk.sample_rate_hz
            block_samples = _samples_per_block(fs, cleaning.block_size_seconds)
            expected_start = chunk.start_index
            segment_start = chunk.start_index
        if chunk.sample_rate_hz != fs:
            raise CleaningError("sample rate changes across eta0 segment chunks")
        if chunk.start_index != expected_start:
            raise CleaningError("eta0 segment chunks must be contiguous and ordered")
        if chunk.start_index % block_samples != 0:
            raise CleaningError("eta0 segment chunk boundary is not block-aligned")

        reduced = _reduce_chunk_blocks(chunk, block_samples, calibration)
        starts_parts.append(reduced[0])
        stops_parts.append(reduced[1])
        time_parts.append(reduced[2])
        eta_avg_parts.append(reduced[5])
        eta_min_parts.append(reduced[6])
        expected_start = chunk.stop_index
        segment_stop = chunk.stop_index

    if fs is None or block_samples is None or segment_start is None or segment_stop is None:
        raise CleaningError("no chunks supplied for eta0 candidate segment")

    starts = np.concatenate(starts_parts)
    stops = np.concatenate(stops_parts)
    timestamps = np.concatenate(time_parts)
    eta_avg = np.concatenate(eta_avg_parts)
    eta_min = np.concatenate(eta_min_parts)

    if np.any(starts[1:] != stops[:-1]):
        raise CleaningError("eta0 candidate reduced blocks contain a gap")

    base = np.isfinite(eta_avg) & np.isfinite(eta_min)
    eta_ref, eta_good = _protected_past_quantile_mask(
        eta_avg,
        eta_min,
        timestamps,
        base,
        window_seconds=cleaning.reference_window_seconds,
        quantile=cleaning.reference_quantile,
        lower_fraction=cleaning.eta0_lower_fraction,
    )

    return APEEtaSegmentPlan(
        sample_rate_hz=fs,
        start_index=segment_start,
        stop_index=segment_stop,
        block_samples=block_samples,
        block_start_indices=starts,
        block_stop_indices=stops,
        block_midpoint_s=timestamps,
        eta0_average_blocks=eta_avg,
        eta0_minimum_blocks=eta_min,
        eta0_reference_blocks=eta_ref,
        eta0_good_blocks=eta_good,
        metadata={
            "method": "ape_eta0_candidate_segment_reset_v1",
            "segment_start_reference_policy": "reset and bootstrap from first valid eta0 block mean",
            "reference_window_seconds": float(cleaning.reference_window_seconds),
            "reference_quantile": float(cleaning.reference_quantile),
            "eta0_lower_fraction": float(cleaning.eta0_lower_fraction),
        },
    )


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
    eta0_avg_parts: list[np.ndarray] = []
    eta0_min_parts: list[np.ndarray] = []

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

        reduced = _reduce_chunk_blocks(chunk, block_samples, calibration)
        starts_parts.append(reduced[0])
        stops_parts.append(reduced[1])
        time_parts.append(reduced[2])
        cavity_avg_parts.append(reduced[3])
        cavity_min_parts.append(reduced[4])
        eta0_avg_parts.append(reduced[5])
        eta0_min_parts.append(reduced[6])
        expected_start = chunk.stop_index

    if not saw_chunk or fs is None or block_samples is None:
        raise CleaningError("no chunks supplied by chunk_factory")

    starts = np.concatenate(starts_parts)
    stops = np.concatenate(stops_parts)
    timestamps = np.concatenate(time_parts)
    cavity_avg = np.concatenate(cavity_avg_parts)
    cavity_min = np.concatenate(cavity_min_parts)
    eta0_avg = np.concatenate(eta0_avg_parts)
    eta0_min = np.concatenate(eta0_min_parts)

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

    eta_base = cavity_good & np.isfinite(eta0_avg) & np.isfinite(eta0_min)
    eta_reference, eta_good = _protected_past_quantile_mask_by_segments(
        eta0_avg,
        eta0_min,
        timestamps,
        eta_base,
        cavity_good,
        window_seconds=cleaning.reference_window_seconds,
        quantile=cleaning.reference_quantile,
        lower_fraction=cleaning.eta0_lower_fraction,
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
        eta0_blocks=eta0_avg,
        eta0_minimum_blocks=eta0_min,
        eta0_reference_blocks=eta_reference,
        eta0_good_blocks=eta_good,
        final_good_blocks=final_good,
        metadata={
            "method": "ape_reduced_q90_protected_eta_cleaning_v3",
            "source_inspiration": "validated earlier APE cavity_lock / pem_modulation workflow",
            "block_size_seconds": float(cleaning.block_size_seconds),
            "reference_window_seconds": float(cleaning.reference_window_seconds),
            "reference_quantile": float(cleaning.reference_quantile),
            "refinement_iterations": int(cleaning.refinement_iterations),
            "cavity_lower_fraction": float(cleaning.cavity_lower_fraction),
            "eta0_lower_fraction": float(cleaning.eta0_lower_fraction),
            "cavity_reference_statistic": "reduced-block average",
            "cavity_decision_statistic": "reduced-block minimum",
            "eta0_reference_statistic": "0.05 s mean of sample-level eta0",
            "eta0_decision_statistic": "0.05 s minimum of sample-level eta0",
            "cavity_window_definition": "abs(t_j - t_i) <= reference_window_seconds / 2",
            "eta0_window_definition": "accepted history with t_k - window_seconds <= t_j < t_k",
            "eta0_reference_update": "accepted cavity-good and eta0-good block means only",
            "eta0_reference_hold": "hold last healthy reference when accepted-history window is empty inside a cavity-good segment",
            "eta0_segment_reset": "reset protected eta0 history at each cavity-good segment start",
            "eta0_bootstrap": "first valid eta0 block mean of each cavity-good segment",
            "quantile_method": "linear interpolation at (n - 1) * q",
            "reference_interpolation": "none",
            "reference_fallback": "eta0 holds last healthy reference only",
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
