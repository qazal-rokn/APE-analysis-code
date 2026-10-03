"""Strict YAML loading for manual and chronology-evidenced automatic modes."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

class IngestionConfigError(ValueError): pass

@dataclass(frozen=True)
class IngestionConfig:
    mode: str
    input_root: Path
    source_manifest: Path | None
    reference_project: Path
    output_directory: Path
    channel_mapping: dict[str, str]
    diagnostic_channels: tuple[str, ...]
    raw_chunk_samples: int
    display_max_points: int
    cleaning_settings: dict[str, Any]
    plotting_settings: dict[str, Any]

def _dict(raw: dict, name: str) -> dict:
    value=raw.get(name)
    if not isinstance(value,dict): raise IngestionConfigError(f"{name} must be a mapping")
    return value

def load_ingestion_config(path: str | Path) -> IngestionConfig:
    path=Path(path).expanduser().resolve()
    with path.open() as f: raw=yaml.safe_load(f)
    if not isinstance(raw,dict): raise IngestionConfigError("configuration must be a YAML mapping")
    mode=str(raw.get("mode")); input_cfg=_dict(raw,"input"); review=_dict(raw,"review")
    if mode not in {"manual","automatic"}: raise IngestionConfigError("mode must be manual or automatic")
    source_manifest=input_cfg.get("source_manifest")
    if mode=="automatic" and not source_manifest: raise IngestionConfigError("automatic mode requires an ordered source_manifest")
    if not bool(input_cfg.get("require_explicit_chronology_evidence",False)):
        raise IngestionConfigError("chronology evidence must be explicitly required")
    cleaning=_dict(raw,"cleaning_settings")
    plotting=_dict(raw,"plotting")
    if int(plotting["display_max_points"]) <= 0 or int(plotting["zoom_max_samples"]) <= 0 or int(plotting["dpi"]) <= 0: raise IngestionConfigError("plotting limits must be positive")
    return IngestionConfig(mode, Path(input_cfg["root"]).expanduser(), None if not source_manifest else Path(source_manifest).expanduser(), Path(raw["reference_project"]).expanduser(), Path(raw["output_directory"]).expanduser(), dict(_dict(raw,"channel_mapping")), tuple(str(x) for x in review["diagnostic_channels"]), int(review["raw_chunk_samples"]), int(review["display_max_points"]), cleaning, plotting)
