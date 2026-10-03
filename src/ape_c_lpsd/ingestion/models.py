from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

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
    def sample_count(self) -> int: return self.raw_stop - self.raw_start
