"""Pure no-gap valid-run selection and immutable segment-manifest writing."""
from __future__ import annotations
import hashlib,json,uuid
from pathlib import Path
from .models import ManualPart, ValidRun

class SegmentSelectionError(RuntimeError): pass

def choose_global_longest(runs: tuple[ValidRun,...]) -> ValidRun:
    if not runs: raise SegmentSelectionError("no valid runs")
    # Inputs are already split at every invalid sample, unverified boundary, gap, or overlap.
    return max(runs,key=lambda r:(r.sample_count,-r.source_ordinal,-r.raw_start))

def freeze_manual_segment(parts: tuple[ManualPart,...], *, sample_rate_hz: float, unit: str="V") -> dict:
    if not parts or sample_rate_hz<=0 or unit!="V": raise SegmentSelectionError("manual manifest requires non-empty V parts and positive fs")
    for p in parts:
        if p.raw_start<0 or p.raw_stop<=p.raw_start: raise SegmentSelectionError("invalid manual part")
    return {"kind":"frozen_validated_segment","selection_mode":"manual","sample_rate_hz":sample_rate_hz,"unit":unit,"parts":[{"path":str(p.path),"raw_start":p.raw_start,"raw_stop":p.raw_stop,"sample_count":p.raw_stop-p.raw_start} for p in parts],"sample_count":sum(p.raw_stop-p.raw_start for p in parts),"continuity":"researcher-supplied explicit ordered slices; no interpolation/filling"}

def write_frozen_segment_manifest(manifest: dict, output_directory: Path) -> tuple[Path,Path]:
    output_directory.mkdir(parents=True,exist_ok=True); token=uuid.uuid4().hex; y=output_directory/f"frozen-segment-{token}.json"
    if y.exists(): raise FileExistsError(f"refusing to overwrite {y}")
    text=json.dumps(manifest,indent=2,sort_keys=True)+"\n"; y.write_text(text)
    digest=hashlib.sha256(y.read_bytes()).hexdigest(); provenance=y.with_suffix(".sha256.json")
    provenance.write_text(json.dumps({"manifest":str(y),"sha256":digest,"method":"manual_or_automatic_frozen_validated_segment"},indent=2)+"\n")
    return y,provenance

from dataclasses import dataclass
@dataclass(frozen=True)
class GlobalSegment:
    parts: tuple[ValidRun,...]
    @property
    def sample_count(self)->int: return sum(p.sample_count for p in self.parts)

def build_global_segments(runs: tuple[ValidRun,...], *, source_sample_counts: dict[int,int]) -> tuple[GlobalSegment,...]:
    """Join only endpoint-to-endpoint runs across an explicitly verified boundary."""
    ordered=sorted(runs,key=lambda r:(r.source_ordinal,r.raw_start)); result=[]; active=[]
    for run in ordered:
        if not active: active=[run]; continue
        previous=active[-1]; can_join=(run.source_ordinal==previous.source_ordinal+1 and previous.raw_stop==source_sample_counts[previous.source_ordinal] and run.raw_start==0 and run.boundary_from_previous_verified)
        if can_join: active.append(run)
        else: result.append(GlobalSegment(tuple(active))); active=[run]
    if active: result.append(GlobalSegment(tuple(active)))
    return tuple(result)

def choose_longest_global_segment(segments: tuple[GlobalSegment,...]) -> GlobalSegment:
    if not segments: raise SegmentSelectionError("no global valid segments")
    return max(segments,key=lambda s:(s.sample_count,-s.parts[0].source_ordinal,-s.parts[0].raw_start))
