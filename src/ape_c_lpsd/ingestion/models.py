from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class APEDXDData:
    """Canonical raw-data object produced by the APE Dewesoft reader.

    Arrays are raw recorded quantities. No calibration, filtering, cleaning,
    interpolation, resampling, trimming, or gap filling has been applied.
    """

    source_file: Path
    sample_rate_hz: float | None
    n_samples: int
    start_time: datetime | None
    x_1f: np.ndarray
    y_1f: np.ndarray
    x_2f: np.ndarray
    y_2f: np.ndarray
    cavity_transmission: np.ndarray
    channel_names: dict[str, str]
    channel_units: dict[str, str]
    regular_sampling: bool
    relative_timestamps_s: np.ndarray | None
    metadata: dict[str, Any]

    def __post_init__(self) -> None:
        if self.n_samples <= 0:
            raise ValueError("n_samples must be positive")
        arrays = (
            self.x_1f,
            self.y_1f,
            self.x_2f,
            self.y_2f,
            self.cavity_transmission,
        )
        if any(np.asarray(values).shape != (self.n_samples,) for values in arrays):
            raise ValueError("all five APE channels must have shape (n_samples,)")
        if self.regular_sampling and (
            self.sample_rate_hz is None or self.sample_rate_hz <= 0
        ):
            raise ValueError("regular sampling requires a positive sample_rate_hz")
        if not self.regular_sampling and self.relative_timestamps_s is None:
            raise ValueError("irregular sampling requires preserved timestamps")


@dataclass(frozen=True)
class OrderedSource:
    path: Path
    ordinal: int
    chronology_evidence: str
    boundary_to_next_verified: bool
    source_sha256: str | None = None


@dataclass(frozen=True)
class ManualPart:
    path: Path
    raw_start: int
    raw_stop: int


@dataclass(frozen=True)
class ValidRun:
    source_ordinal: int
    raw_start: int
    raw_stop: int
    sample_rate_hz: float
    boundary_from_previous_verified: bool
    source_path: Path

    def __post_init__(self):
        if self.raw_start < 0 or self.raw_stop <= self.raw_start or self.sample_rate_hz <= 0:
            raise ValueError("invalid valid-run interval")

    @property
    def sample_count(self) -> int:
        return self.raw_stop - self.raw_start
