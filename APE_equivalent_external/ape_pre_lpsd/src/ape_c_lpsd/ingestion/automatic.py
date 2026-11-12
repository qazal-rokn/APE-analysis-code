"""Future bounded automatic runner; chronology and reference workflow are mandatory."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import numpy as np
from .chronology import read_ordered_source_manifest
from .models import ValidRun
from .reference_workflow import ReferenceWorkflow

@dataclass(frozen=True)
class FileMaskSummary:
    source_path: Path; raw_sample_count: int; valid_count: int; sample_rate_hz: float

class AutomaticModeRefused(RuntimeError): pass

def contiguous_runs(mask: np.ndarray, *, source_ordinal:int, source_path:Path, sample_rate_hz:float, boundary_from_previous_verified:bool) -> tuple[ValidRun,...]:
    """Split at every false raw sample; never creates data or bridges a false value."""
    if mask.ndim!=1: raise ValueError("valid mask must be one-dimensional")
    padded=np.concatenate(([False],np.asarray(mask,dtype=bool),[False])); edges=np.flatnonzero(padded[1:]!=padded[:-1]); result=[]
    for start,stop in zip(edges[::2],edges[1::2]): result.append(ValidRun(source_ordinal,int(start),int(stop),sample_rate_hz,boundary_from_previous_verified,source_path))
    return tuple(result)

def automatic_mode_contract(*, source_manifest:Path, input_root:Path, workflow:ReferenceWorkflow) -> tuple:
    """Validate mandatory chronology and frozen-workflow availability before any file read."""
    return read_ordered_source_manifest(source_manifest,input_root=input_root)

def run_automatic_ingestion(*args, **kwargs):
    """Reserved for a separately approved execution turn; never implicitly runs."""
    raise AutomaticModeRefused("automatic ingestion execution requires separate approval after synthetic validation")

def validate_source_metadata(*, raw_sample_count:int, observed_raw_sample_count:int, sample_rate_hz:float, expected_sample_rate_hz:float, source_path:Path) -> None:
    """Fail closed on missing observations, inconsistent counts, or rate changes."""
    if raw_sample_count<=0 or observed_raw_sample_count!=raw_sample_count:
        raise AutomaticModeRefused(f"source observation-count mismatch: {source_path}")
    if sample_rate_hz!=expected_sample_rate_hz:
        raise AutomaticModeRefused(f"source sample-rate mismatch: {source_path}")
