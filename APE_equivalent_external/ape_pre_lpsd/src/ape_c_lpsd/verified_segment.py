"""Bounded, read-only validation of the approved three-part voltage segment."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from time import monotonic
from typing import Iterable

import numpy as np

from ape_analysis.calibration import PolarizationCalibrationParameters, RawPemQuadratureChannel, build_polarization_calibration
from ape_analysis.io import DewesoftRecordedDatasetSource, DewesoftReducedDatasetSource
from ape_analysis.quality import CavityLockParameters, ExternallyReducedCavityRecord, PemModulationParameters, ReducedAverageChannel, build_cavity_lock_mask, build_pem_modulation_quality

SAMPLES_PER_REDUCED_RECORD = 1_000
EXPECTED_SAMPLE_RATE_HZ = 20_000.0
EXPECTED_TOTAL_SAMPLES = 308_717_000


@dataclass(frozen=True)
class Part:
    path: Path
    raw_start: int
    raw_stop: int


@dataclass(frozen=True)
class ValidationReport:
    sample_count: int
    chunks: int
    elapsed_seconds: float
    part_counts: tuple[int, ...]


def _reduced(path: Path, name: str, count: int, *, dewesoft_paths=None):
    kwargs = {} if dewesoft_paths is None else dewesoft_paths.source_kwargs()
    return DewesoftReducedDatasetSource(path, reduced_record_start=0, reduced_record_count=count, **kwargs).read_reduced_channel(name)


def _quality(path: Path, count: int, *, dewesoft_paths=None):
    """Unchanged Module 2A/2B sequence used by the reference preview."""
    cavity, pem2_x, pem2_y = (_reduced(path, n, count, dewesoft_paths=dewesoft_paths) for n in ("cavity_trans", "pem2_X", "pem2_Y"))
    source, start = str(path), float(cavity.metadata["start_measure_time"])
    lock = build_cavity_lock_mask(tuple(ExternallyReducedCavityRecord(timestamp=r.dewesoft_timestamp, cavity_average=r.ave, cavity_minimum=r.min, source_file=source, start_measure_time=start) for r in cavity.records), CavityLockParameters())

    def channel(name, values):
        return ReducedAverageChannel(recorded_name=name, timestamps=tuple(r.dewesoft_timestamp for r in values.records), averages=tuple(r.ave for r in values.records), block_size_seconds=values.block_size_seconds, source_file=source, start_measure_time=start)

    return build_pem_modulation_quality(channel("pem2_X", pem2_x), channel("pem2_Y", pem2_y), channel("cavity_trans", cavity), lock.cavity_good, PemModulationParameters())


def _read(path: Path, start: int, count: int, name: str, *, dewesoft_paths=None):
    kwargs = {} if dewesoft_paths is None else dewesoft_paths.source_kwargs()
    return DewesoftRecordedDatasetSource(path, sample_start=start, sample_count=count, **kwargs).read_recorded_channel(name)


def _check_metadata(x, y, part: Part):
    if x.metadata.get("unit") != "V" or y.metadata.get("unit") != "V":
        raise RuntimeError(f"unit refusal {part.path.name}: pem1_x={x.metadata.get('unit')!r}, pem1_Y={y.metadata.get('unit')!r}")
    if x.sample_rate != EXPECTED_SAMPLE_RATE_HZ or y.sample_rate != EXPECTED_SAMPLE_RATE_HZ:
        raise RuntimeError(f"sample-rate refusal in {part.path.name}")


def _mask(x: np.ndarray, y: np.ndarray, raw_start: int, quality) -> np.ndarray:
    """Vectorized evaluation of the frozen Module 3 validity definition."""
    p = PolarizationCalibrationParameters()
    parent = np.arange(raw_start, raw_start + len(x), dtype=np.int64) // SAMPLES_PER_REDUCED_RECORD
    i0 = np.asarray([r.i0_w for r in quality.records], dtype=np.float64)[parent]
    eta = np.asarray([r.eta0 for r in quality.records], dtype=np.float64)[parent]
    good = np.asarray(quality.combined_good, dtype=bool)[parent]
    magnitude = np.hypot(x, y)
    denominator = p.first_harmonic_rms_factor * i0 * eta * p.responsivity_a_per_w * p.transimpedance_gain_v_per_a
    rho = magnitude / denominator
    return good & np.isfinite(x) & np.isfinite(y) & np.isfinite(magnitude) & np.isfinite(i0) & np.isfinite(eta) & (i0 > 0) & (eta > 0) & np.isfinite(denominator) & (denominator > 0) & np.isfinite(rho)


def validate_reference_equivalence(part: Part, *, sample_count: int = 100_000) -> None:
    """Compare the chunk mask to unmodified Module 3 on an origin-zero prefix."""
    if part.raw_start != 0:
        raise ValueError("reference equivalence requires an origin-zero part")
    count = min(sample_count, part.raw_stop) // 1000 * 1000
    quality = _quality(part.path, count // 1000)
    raw_x, raw_y = _read(part.path, 0, count, "pem1_x"), _read(part.path, 0, count, "pem1_Y")
    _check_metadata(raw_x, raw_y, part)
    x, y = np.asarray(raw_x.samples, dtype=np.float64), np.asarray(raw_y.samples, dtype=np.float64)
    reference = build_polarization_calibration(RawPemQuadratureChannel("pem1_x", x, EXPECTED_SAMPLE_RATE_HZ, str(part.path), 0.0), RawPemQuadratureChannel("pem1_Y", y, EXPECTED_SAMPLE_RATE_HZ, str(part.path), 0.0), tuple(r.i0_w for r in quality.records), tuple(r.eta0 for r in quality.records), quality.combined_good, PolarizationCalibrationParameters())
    expected = np.fromiter((r.sample_good and r.calibration_valid for r in reference.records), dtype=bool, count=count)
    if not np.array_equal(_mask(x, y, 0, quality), expected):
        raise AssertionError("bounded validity mask differs from frozen Module 3")


def validate(parts: Iterable[Part], *, chunk_samples: int = 1_000_000) -> ValidationReport:
    """Validate all retained samples without materializing the full arrays."""
    if chunk_samples <= 0 or chunk_samples % SAMPLES_PER_REDUCED_RECORD:
        raise ValueError("chunk_samples must be positive and divisible by 1000")
    started, total, chunks, counts = monotonic(), 0, 0, []
    for part in tuple(parts):
        if part.raw_start < 0 or part.raw_stop <= part.raw_start:
            raise ValueError(f"invalid slice for {part.path.name}")
        quality = _quality(part.path, ((part.raw_stop + 999) // 1000))
        part_count = 0
        for start in range(part.raw_start, part.raw_stop, chunk_samples):
            count = min(chunk_samples, part.raw_stop - start)
            raw_x, raw_y = _read(part.path, start, count, "pem1_x"), _read(part.path, start, count, "pem1_Y")
            _check_metadata(raw_x, raw_y, part)
            x, y = np.asarray(raw_x.samples, dtype=np.float64), np.asarray(raw_y.samples, dtype=np.float64)
            valid = _mask(x, y, start, quality)
            if not np.all(valid):
                raise RuntimeError(f"invalid selected sample in {part.path.name} at raw index {start + int(np.flatnonzero(~valid)[0])}")
            total += count
            part_count += count
            chunks += 1
        counts.append(part_count)
    if total != EXPECTED_TOTAL_SAMPLES:
        raise RuntimeError(f"approved segment count mismatch: {total}")
    return ValidationReport(total, chunks, monotonic() - started, tuple(counts))
