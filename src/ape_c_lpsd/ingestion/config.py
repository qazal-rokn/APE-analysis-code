"""Configuration loading for APE Dewesoft ingestion and legacy review modes."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


class IngestionConfigError(ValueError):
    pass


@dataclass(frozen=True)
class APEChannelNames:
    x_1f: str
    y_1f: str
    x_2f: str
    y_2f: str
    cavity_transmission: str

    def as_mapping(self) -> dict[str, str]:
        return {
            "x_1f": self.x_1f,
            "y_1f": self.y_1f,
            "x_2f": self.x_2f,
            "y_2f": self.y_2f,
            "cavity_transmission": self.cavity_transmission,
        }


@dataclass(frozen=True)
class CalibrationConfig:
    """Run-specific calibration parameters.

    These values are configuration, not universal APE constants. This ingestion
    stage records them but does not apply calibration.
    """

    cavity_transmission_w_per_v: float
    responsivity_a_per_w: float
    transimpedance_gain_v_per_a: float
    second_harmonic_pd_response: float


@dataclass(frozen=True)
class DXDIngestionConfig:
    input_file: Path
    reference_project: Path
    output_directory: Path | None
    channels: APEChannelNames
    calibration: CalibrationConfig


@dataclass(frozen=True)
class IngestionConfig:
    """Legacy chronology/quality-review configuration kept for cross-checks."""

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
    value = raw.get(name)
    if not isinstance(value, dict):
        raise IngestionConfigError(f"{name} must be a mapping")
    return value


def _required_text(mapping: dict, key: str, section: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise IngestionConfigError(f"{section}.{key} must be a non-empty string")
    return value.strip()


def _required_positive_float(mapping: dict, key: str, section: str) -> float:
    if key not in mapping:
        raise IngestionConfigError(f"{section}.{key} is required")
    try:
        value = float(mapping[key])
    except (TypeError, ValueError) as exc:
        raise IngestionConfigError(f"{section}.{key} must be numeric") from exc
    if value <= 0:
        raise IngestionConfigError(f"{section}.{key} must be positive")
    return value


def load_dxd_ingestion_config(path: str | Path) -> DXDIngestionConfig:
    """Load the canonical first-stage .dxd ingestion configuration."""
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise IngestionConfigError("configuration must be a YAML mapping")

    input_file = raw.get("input_file")
    reference_project = raw.get("reference_project")
    if not isinstance(input_file, str) or not input_file.strip():
        raise IngestionConfigError("input_file must be a non-empty path string")
    if not isinstance(reference_project, str) or not reference_project.strip():
        raise IngestionConfigError("reference_project must be a non-empty path string")

    channels_raw = _dict(raw, "channels")
    channels = APEChannelNames(
        x_1f=_required_text(channels_raw, "x_1f", "channels"),
        y_1f=_required_text(channels_raw, "y_1f", "channels"),
        x_2f=_required_text(channels_raw, "x_2f", "channels"),
        y_2f=_required_text(channels_raw, "y_2f", "channels"),
        cavity_transmission=_required_text(
            channels_raw, "cavity_transmission", "channels"
        ),
    )
    mapped = list(channels.as_mapping().values())
    if len(set(mapped)) != len(mapped):
        raise IngestionConfigError(
            "the five APE roles must map to five distinct Dewesoft channels"
        )

    calibration_raw = _dict(raw, "calibration")
    calibration = CalibrationConfig(
        cavity_transmission_w_per_v=_required_positive_float(
            calibration_raw, "cavity_transmission_w_per_v", "calibration"
        ),
        responsivity_a_per_w=_required_positive_float(
            calibration_raw, "responsivity_a_per_w", "calibration"
        ),
        transimpedance_gain_v_per_a=_required_positive_float(
            calibration_raw, "transimpedance_gain_v_per_a", "calibration"
        ),
        second_harmonic_pd_response=_required_positive_float(
            calibration_raw, "second_harmonic_pd_response", "calibration"
        ),
    )

    output = raw.get("output_directory")
    if output is not None and (not isinstance(output, str) or not output.strip()):
        raise IngestionConfigError(
            "output_directory must be omitted/null or a non-empty path string"
        )

    return DXDIngestionConfig(
        input_file=Path(input_file).expanduser(),
        reference_project=Path(reference_project).expanduser(),
        output_directory=None if output is None else Path(output).expanduser(),
        channels=channels,
        calibration=calibration,
    )


def load_ingestion_config(path: str | Path) -> IngestionConfig:
    """Load the legacy chronology-gated ingestion/review configuration."""
    path = Path(path).expanduser().resolve()
    with path.open() as f:
        raw = yaml.safe_load(f)
    if not isinstance(raw, dict):
        raise IngestionConfigError("configuration must be a YAML mapping")
    mode = str(raw.get("mode"))
    input_cfg = _dict(raw, "input")
    review = _dict(raw, "review")
    if mode not in {"manual", "automatic"}:
        raise IngestionConfigError("mode must be manual or automatic")
    source_manifest = input_cfg.get("source_manifest")
    if mode == "automatic" and not source_manifest:
        raise IngestionConfigError("automatic mode requires an ordered source_manifest")
    if not bool(input_cfg.get("require_explicit_chronology_evidence", False)):
        raise IngestionConfigError("chronology evidence must be explicitly required")
    cleaning = _dict(raw, "cleaning_settings")
    plotting = _dict(raw, "plotting")
    if (
        int(plotting["display_max_points"]) <= 0
        or int(plotting["zoom_max_samples"]) <= 0
        or int(plotting["dpi"]) <= 0
    ):
        raise IngestionConfigError("plotting limits must be positive")
    return IngestionConfig(
        mode,
        Path(input_cfg["root"]).expanduser(),
        None if not source_manifest else Path(source_manifest).expanduser(),
        Path(raw["reference_project"]).expanduser(),
        Path(raw["output_directory"]).expanduser(),
        dict(_dict(raw, "channel_mapping")),
        tuple(str(x) for x in review["diagnostic_channels"]),
        int(review["raw_chunk_samples"]),
        int(review["display_max_points"]),
        cleaning,
        plotting,
    )
