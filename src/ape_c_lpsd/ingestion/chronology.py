"""Manifest-only chronology gate; never infers order from names or mtimes."""
from __future__ import annotations
import hashlib
from pathlib import Path
import yaml
from .models import OrderedSource

class ChronologyError(RuntimeError): pass

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for c in iter(lambda:f.read(1<<20),b""): h.update(c)
    return h.hexdigest()

def read_ordered_source_manifest(path: str|Path, *, input_root: Path) -> tuple[OrderedSource,...]:
    with Path(path).open() as f: raw=yaml.safe_load(f)
    if not isinstance(raw,dict) or raw.get("chronology_basis") != "researcher_provided_or_explicit_manifest":
        raise ChronologyError("source manifest lacks explicit chronology basis")
    entries=raw.get("sources")
    if not isinstance(entries,list) or not entries: raise ChronologyError("source manifest needs ordered sources")
    result=[]
    for ordinal,item in enumerate(entries):
        if not isinstance(item,dict) or not item.get("chronology_evidence"):
            raise ChronologyError(f"source {ordinal} lacks chronology evidence")
        p=input_root/str(item["relative_path"])
        if not p.is_file(): raise ChronologyError(f"missing source file: {p}")
        claimed=item.get("sha256")
        if claimed and sha256(p)!=claimed: raise ChronologyError(f"source checksum mismatch: {p}")
        result.append(OrderedSource(p,ordinal,str(item["chronology_evidence"]),bool(item.get("boundary_to_next_verified",False)),claimed))
    return tuple(result)
