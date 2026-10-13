from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class APEDXDMetadata:
    """Lightweight file-level metadata for one canonical APE .dxd source."""

    source_file: Path
    sample_rate_hz: float
    n_samples: int
    duration_s: float
    start_time: datetime | None
    channel_names: dict[str, str]
    channel_units: dict[str, str]
    dewesoft_channel_count: int
    metadata: dict[str, Any]


@dataclass(frozen=True)
class APEDXDInspection:
    """Result of full-file bounded-memory validation."""

    source_file: Path
    sample_rate_hz: float
    n_samples: int
    duration_s: float
    start_time: datetime | None
    channel_names: dict[str, str]
    channel_units: dict[str, str]
    dewesoft_channel_count: int
    chunk_size_samples: int
    chunk_count: int
    metadata: dict[str, Any]


@dataclass(frozen=True)
class APEDXDChunk:
    """One synchronized raw chunk from the five required APE channels."""

    source_file: Path
    start_index: int
    stop_index: int
    sample_rate_hz: float
    x_1f: np.ndarray
    y_1f: np.ndarray
    x_2f: np.ndarray
    y_2f: np.ndarray
    cavity_transmission: np.ndarray
    first_timestamp_s: float
    last_timestamp_s: float

    def __post_init__(self) -> None:
        count = self.stop_index - self.start_index
        if self.start_index < 0 or count <= 0:
            raise ValueError("invalid chunk interval")
        arrays = (
            self.x_1f,
            self.y_1f,
            self.x_2f,
            self.y_2f,
            self.cavity_transmission,
        )
        if any(np.asarray(values).shape != (count,) for values in arrays):
            raise ValueError("all five chunk channels must match the chunk interval")
        if self.sample_rate_hz <= 0:
            raise ValueError("chunk sample rate must be positive")

    @property
    def n_samples(self) -> int:
        return self.stop_index - self.start_index


@dataclass(frozen=True)
class APECavityChunk:
    """One bounded-memory chunk containing only cavity transmission."""

    source_file: Path
    start_index: int
    stop_index: int
    sample_rate_hz: float
    cavity_transmission: np.ndarray
    first_timestamp_s: float
    last_timestamp_s: float

    def __post_init__(self) -> None:
        count = self.stop_index - self.start_index
        if self.start_index < 0 or count <= 0:
            raise ValueError("invalid cavity chunk interval")
        if np.asarray(self.cavity_transmission).shape != (count,):
            raise ValueError("cavity transmission must match the chunk interval")
        if self.sample_rate_hz <= 0:
            raise ValueError("cavity chunk sample rate must be positive")

    @property
    def n_samples(self) -> int:
        return self.stop_index - self.start_index


@dataclass(frozen=True)
class APEEtaChunk:
    """One bounded-memory chunk containing only channels needed for eta0."""

    source_file: Path
    start_index: int
    stop_index: int
    sample_rate_hz: float
    x_2f: np.ndarray
    y_2f: np.ndarray
    cavity_transmission: np.ndarray
    first_timestamp_s: float
    last_timestamp_s: float

    def __post_init__(self) -> None:
        count = self.stop_index - self.start_index
        if self.start_index < 0 or count <= 0:
            raise ValueError("invalid eta0 chunk interval")
        arrays = (self.x_2f, self.y_2f, self.cavity_transmission)
        if any(np.asarray(values).shape != (count,) for values in arrays):
            raise ValueError("eta0 chunk channels must match the chunk interval")
        if self.sample_rate_hz <= 0:
            raise ValueError("eta0 chunk sample rate must be positive")

    @property
    def n_samples(self) -> int:
        return self.stop_index - self.start_index


@dataclass(frozen=True)
class APECalibratedChunk:
    """Calibrated observables for one raw APE chunk."""

    source_file: Path
    start_index: int
    stop_index: int
    sample_rate_hz: float
    i0_w: np.ndarray
    eta0: np.ndarray
    rho_x: np.ndarray
    rho_y: np.ndarray
    valid_mask: np.ndarray
    calibration_metadata: dict[str, Any]

    def __post_init__(self) -> None:
        count = self.stop_index - self.start_index
        if self.start_index < 0 or count <= 0:
            raise ValueError("invalid calibrated chunk interval")
        arrays = (self.i0_w, self.eta0, self.rho_x, self.rho_y, self.valid_mask)
        if any(np.asarray(values).shape != (count,) for values in arrays):
            raise ValueError("all calibrated arrays must match the chunk interval")
        if self.sample_rate_hz <= 0:
            raise ValueError("calibrated chunk sample rate must be positive")
        if np.asarray(self.valid_mask).dtype != np.bool_:
            raise ValueError("valid_mask must be boolean")

    @property
    def n_samples(self) -> int:
        return self.stop_index - self.start_index


@dataclass(frozen=True)
class APEDXDData:
    """Compatibility full-memory raw-data object for small controlled inputs."""

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
