"""C backend input gate: only frozen validated-segment manifests are acceptable."""
from __future__ import annotations
import json
from pathlib import Path

def require_frozen_validated_segment_manifest(path:str|Path)->dict:
    p=Path(path)
    with p.open() as f: manifest=json.load(f)
    if manifest.get("kind")!="frozen_validated_segment": raise ValueError("C-LPSD backend refuses non-frozen segment input")
    if manifest.get("unit")!="V" or not manifest.get("parts") or int(manifest.get("sample_count",0))<=0: raise ValueError("invalid frozen validated-segment manifest")
    return manifest
