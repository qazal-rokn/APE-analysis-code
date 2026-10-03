"""Approved staged APE cleaning.

Stage 1: cavity transmission one-sided lower cut against a local median.
Stage 2: eta0 one-sided lower cut against a local median, evaluated only on
cavity-good samples.

Local references are formed from short block medians and then a centered
rolling median across those blocks. Raw sample indices are always preserved.
"""
from __future__ import annotations

from dataclasses import dataclass
from bisect import bisect_left, insort
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
    """Full-record local-median references stored at block resolution."""

    sample_rate_hz: float
    total_samples: int
    block_samples: int
    cavity_reference_blocks: np.ndarray
    eta0_reference_blocks: np.ndarray
    metadata: dict[str, object]


def _samples_per_block(sample_rate_hz: float, block_size_seconds: float) -> int:
    value = int(round(sample_rate_hz * block_size_seconds))
    if value <= 0:
        raise CleaningError("cleaning block must contain at least one sample")
    if not np.isclose(value / sample_rate_hz, block_size_seconds, rtol=0.0, atol=1e-12):
        raise CleaningError("block_size_seconds must map to an integer sample count")
    return value


def _block_medians(values: np.ndarray, block_samples: int, valid: np.ndarray | None = None) -> np.ndarray:
    n = values.size
    n_blocks = (n + block_samples - 1) // block_samples
    out = np.full(n_blocks, np.nan, dtype=np.float64)
    for b in range(n_blocks):
        start = b * block_samples
        stop = min(n, start + block_samples)
        block = values[start:stop]
        if valid is not None:
            block = block[valid[start:stop]]
        block = block[np.isfinite(block)]
        if block.size:
            out[b] = float(np.median(block))
    return out


def _rolling_block_median(block_values: np.ndarray, window_blocks: int) -> np.ndarray:
    """Centered rolling median using an incrementally maintained sorted window."""
    if window_blocks <= 0:
        raise CleaningError("reference window must contain at least one block")
    values = np.asarray(block_values, dtype=np.float64)
    n = values.size
    refs = np.full(n, np.nan, dtype=np.float64)
    half_left = (window_blocks - 1) // 2
    half_right = window_blocks // 2
    ordered: list[float] = []
    previous_lo = 0
    previous_hi = 0

    for i in range(n):
        lo = max(0, i - half_left)
        hi = min(n, i + half_right + 1)

        while previous_lo < lo:
            value = values[previous_lo]
            if np.isfinite(value):
                pos = bisect_left(ordered, float(value))
                if pos >= len(ordered) or ordered[pos] != float(value):
                    raise CleaningError("rolling-median window bookkeeping failed")
                ordered.pop(pos)
            previous_lo += 1

        while previous_hi < hi:
            value = values[previous_hi]
            if np.isfinite(value):
                insort(ordered, float(value))
            previous_hi += 1

        if ordered:
            m = len(ordered)
            if m % 2:
                refs[i] = ordered[m // 2]
            else:
                refs[i] = 0.5 * (ordered[m // 2 - 1] + ordered[m // 2])

    return refs


def _expand_block_values(block_values: np.ndarray, n_samples: int, block_samples: int) -> np.ndarray:
    return np.repeat(block_values, block_samples)[:n_samples]


def _eta0_on_mask(
    chunk: APEDXDChunk,
    calibration: CalibrationConfig,
    mask: np.ndarray,
) -> np.ndarray:
    c_trans = float(calibration.cavity_transmission_w_per_v)
    q = float(calibration.responsivity_a_per_w)
    gain = float(calibration.transimpedance_gain_v_per_a)
    h2 = float(calibration.second_harmonic_pd_response)
    for name, value in {
        "cavity_transmission_w_per_v": c_trans,
        "responsivity_a_per_w": q,
        "transimpedance_gain_v_per_a": gain,
        "second_harmonic_pd_response": h2,
    }.items():
        if not np.isfinite(value) or value <= 0:
            raise CleaningError(f"{name} must be positive and finite")

    i0 = c_trans * np.asarray(chunk.cavity_transmission, dtype=np.float64)
    r2f = np.hypot(
        np.asarray(chunk.x_2f, dtype=np.float64),
        np.asarray(chunk.y_2f, dtype=np.float64),
    )
    valid = (
        np.asarray(mask, dtype=bool)
        & np.isfinite(i0)
        & np.isfinite(r2f)
        & (i0 > 0.0)
        & (r2f > 0.0)
    )
    eta0 = np.full(chunk.n_samples, np.nan, dtype=np.float64)
    arg = np.full(chunk.n_samples, np.nan, dtype=np.float64)
    arg[valid] = 2.0 * np.sqrt(2.0) * r2f[valid] / (h2 * q * gain * i0[valid])
    valid &= np.isfinite(arg) & (arg > 0.0)
    eta0[valid] = np.sqrt(arg[valid])
    return eta0


def _reference_for_chunk(
    reference_blocks: np.ndarray,
    chunk: APEDXDChunk,
    block_samples: int,
) -> np.ndarray:
    if chunk.start_index % block_samples != 0:
        raise CleaningError(
            "production cleaning requires chunk starts aligned to cleaning blocks"
        )
    first_block = chunk.start_index // block_samples
    n_blocks = (chunk.n_samples + block_samples - 1) // block_samples
    selected = reference_blocks[first_block:first_block + n_blocks]
    if selected.size != n_blocks:
        raise CleaningError("cleaning reference does not cover this chunk")
    return _expand_block_values(selected, chunk.n_samples, block_samples)


def _validate_chunk_sequence(
    chunks: list[APEDXDChunk],
    block_samples: int,
) -> tuple[float, int]:
    if not chunks:
        raise CleaningError("no chunks supplied")
    fs = chunks[0].sample_rate_hz
    expected_start = 0
    for i, chunk in enumerate(chunks):
        if chunk.sample_rate_hz != fs:
            raise CleaningError("sample rate changes across cleaning chunks")
        if chunk.start_index != expected_start:
            raise CleaningError("cleaning chunks must be contiguous and ordered")
        if i < len(chunks) - 1 and chunk.n_samples % block_samples != 0:
            raise CleaningError(
                "all non-final cleaning chunks must end on a cleaning-block boundary"
            )
        expected_start = chunk.stop_index
    return fs, expected_start


def build_cleaning_reference_plan(
    chunk_factory: Callable[[], Iterator[APEDXDChunk]],
    cleaning: CleaningConfig,
    calibration: CalibrationConfig,
) -> APECleaningReferencePlan:
    """Build full-record cavity and eta0 local references in two bounded-memory passes.

    Only block medians/references are retained for the complete record. Raw
    science channels remain chunked.
    """
    first_pass_chunks = list(chunk_factory())
    if not first_pass_chunks:
        raise CleaningError("no chunks supplied by chunk_factory")
    fs = first_pass_chunks[0].sample_rate_hz
    block_samples = _samples_per_block(fs, cleaning.block_size_seconds)
    fs, total_samples = _validate_chunk_sequence(first_pass_chunks, block_samples)

    window_blocks_float = cleaning.reference_window_seconds / cleaning.block_size_seconds
    window_blocks = int(round(window_blocks_float))
    if not np.isclose(window_blocks, window_blocks_float, rtol=0.0, atol=1e-12):
        raise CleaningError("reference_window_seconds must be an integer number of blocks")

    cavity_block_parts = [
        _block_medians(
            np.asarray(chunk.cavity_transmission, dtype=np.float64),
            block_samples,
        )
        for chunk in first_pass_chunks
    ]
    cavity_blocks = np.concatenate(cavity_block_parts)
    cavity_reference_blocks = _rolling_block_median(cavity_blocks, window_blocks)

    # Release first-pass raw chunks before the eta0 pass.
    del first_pass_chunks

    eta_block_parts: list[np.ndarray] = []
    expected_start = 0
    for chunk in chunk_factory():
        if chunk.start_index != expected_start:
            raise CleaningError("second cleaning pass changed chunk ordering")
        cavity_reference = _reference_for_chunk(
            cavity_reference_blocks, chunk, block_samples
        )
        vtrans = np.asarray(chunk.cavity_transmission, dtype=np.float64)
        cavity_good = (
            np.isfinite(vtrans)
            & np.isfinite(cavity_reference)
            & (cavity_reference > 0.0)
            & (vtrans >= cleaning.cavity_lower_fraction * cavity_reference)
        )
        eta0 = _eta0_on_mask(chunk, calibration, cavity_good)
        eta_block_parts.append(
            _block_medians(eta0, block_samples, valid=cavity_good)
        )
        expected_start = chunk.stop_index

    if expected_start != total_samples:
        raise CleaningError("second cleaning pass sample count changed")
    eta_blocks = np.concatenate(eta_block_parts)
    if eta_blocks.size != cavity_reference_blocks.size:
        raise CleaningError("cavity and eta0 cleaning block counts disagree")
    eta0_reference_blocks = _rolling_block_median(eta_blocks, window_blocks)

    return APECleaningReferencePlan(
        sample_rate_hz=fs,
        total_samples=total_samples,
        block_samples=block_samples,
        cavity_reference_blocks=cavity_reference_blocks,
        eta0_reference_blocks=eta0_reference_blocks,
        metadata={
            "method": "ape_full_record_local_median_reference_plan_v1",
            "block_size_seconds": float(cleaning.block_size_seconds),
            "reference_window_seconds": float(cleaning.reference_window_seconds),
            "cavity_lower_fraction": float(cleaning.cavity_lower_fraction),
            "eta0_lower_fraction": float(cleaning.eta0_lower_fraction),
            "raw_record_retained_in_memory": False,
        },
    )


def clean_chunk_with_plan(
    chunk: APEDXDChunk,
    plan: APECleaningReferencePlan,
    cleaning: CleaningConfig,
    calibration: CalibrationConfig,
) -> APECleaningResult:
    """Apply full-record references to one raw chunk."""
    if chunk.sample_rate_hz != plan.sample_rate_hz:
        raise CleaningError("chunk sample rate differs from cleaning plan")
    cavity_reference = _reference_for_chunk(
        plan.cavity_reference_blocks, chunk, plan.block_samples
    )
    vtrans = np.asarray(chunk.cavity_transmission, dtype=np.float64)
    cavity_good = (
        np.isfinite(vtrans)
        & np.isfinite(cavity_reference)
        & (cavity_reference > 0.0)
        & (vtrans >= cleaning.cavity_lower_fraction * cavity_reference)
    )

    eta0 = _eta0_on_mask(chunk, calibration, cavity_good)
    eta0_reference = _reference_for_chunk(
        plan.eta0_reference_blocks, chunk, plan.block_samples
    )
    eta_good = (
        cavity_good
        & np.isfinite(eta0)
        & np.isfinite(eta0_reference)
        & (eta0_reference > 0.0)
        & (eta0 >= cleaning.eta0_lower_fraction * eta0_reference)
    )
    final_good = cavity_good & eta_good

    return APECleaningResult(
        cavity_good=cavity_good,
        eta0_preliminary=eta0,
        eta_good=eta_good,
        final_good=final_good,
        cavity_reference=cavity_reference,
        eta0_reference=eta0_reference,
        metadata={
            **dict(plan.metadata),
            "indices_preserved": True,
            "one_sided_lower_cuts_only": True,
            "rho_not_computed_during_cleaning": True,
        },
    )


def clean_chunk(
    chunk: APEDXDChunk,
    cleaning: CleaningConfig,
    calibration: CalibrationConfig,
) -> APECleaningResult:
    """Apply the approved cavity-first, eta0-second cleaning to one chunk."""
    if not (0.0 < cleaning.cavity_lower_fraction <= 1.0):
        raise CleaningError("cavity_lower_fraction must be in (0, 1]")
    if not (0.0 < cleaning.eta0_lower_fraction <= 1.0):
        raise CleaningError("eta0_lower_fraction must be in (0, 1]")

    block_samples = _samples_per_block(chunk.sample_rate_hz, cleaning.block_size_seconds)
    window_blocks_float = cleaning.reference_window_seconds / cleaning.block_size_seconds
    window_blocks = int(round(window_blocks_float))
    if not np.isclose(window_blocks, window_blocks_float, rtol=0.0, atol=1e-12):
        raise CleaningError("reference_window_seconds must be an integer number of blocks")

    vtrans = np.asarray(chunk.cavity_transmission, dtype=np.float64)
    cavity_blocks = _block_medians(vtrans, block_samples)
    cavity_ref_blocks = _rolling_block_median(cavity_blocks, window_blocks)
    cavity_reference = _expand_block_values(cavity_ref_blocks, chunk.n_samples, block_samples)

    cavity_good = (
        np.isfinite(vtrans)
        & np.isfinite(cavity_reference)
        & (cavity_reference > 0.0)
        & (vtrans >= cleaning.cavity_lower_fraction * cavity_reference)
    )

    eta0 = _eta0_on_mask(chunk, calibration, cavity_good)
    eta_blocks = _block_medians(eta0, block_samples, valid=cavity_good)
    eta_ref_blocks = _rolling_block_median(eta_blocks, window_blocks)
    eta0_reference = _expand_block_values(eta_ref_blocks, chunk.n_samples, block_samples)

    eta_good = (
        cavity_good
        & np.isfinite(eta0)
        & np.isfinite(eta0_reference)
        & (eta0_reference > 0.0)
        & (eta0 >= cleaning.eta0_lower_fraction * eta0_reference)
    )
    final_good = cavity_good & eta_good

    metadata = {
        "method": "ape_local_median_two_stage_cleaning_v1",
        "block_size_seconds": float(cleaning.block_size_seconds),
        "reference_window_seconds": float(cleaning.reference_window_seconds),
        "cavity_lower_fraction": float(cleaning.cavity_lower_fraction),
        "eta0_lower_fraction": float(cleaning.eta0_lower_fraction),
        "one_sided_lower_cuts_only": True,
        "slow_drift_preserved_by_local_reference": True,
        "indices_preserved": True,
        "rho_not_computed_during_cleaning": True,
    }

    return APECleaningResult(
        cavity_good=cavity_good,
        eta0_preliminary=eta0,
        eta_good=eta_good,
        final_good=final_good,
        cavity_reference=cavity_reference,
        eta0_reference=eta0_reference,
        metadata=metadata,
    )
