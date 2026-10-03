"""Configuration contract for the independent full-record FFT extension."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

class FFTConfigurationError(ValueError):
    pass

@dataclass(frozen=True)
class FFTAnalysisOptions:
    requested_display_min_hz: float
    requested_display_max_hz: float
    fft_min_hz: float
    fft_max_hz: float
    c_lpsd_handoff_hz: float
    c_lpsd_min_hz: float
    c_lpsd_max_hz: float
    detrend: str
    window_policy: str
    kaiser_beta: float | None
    one_sided: bool
    full_record_only: bool
    zero_padding_factor: int
    memory_limit_bytes: int

@dataclass(frozen=True)
class FFTArtifactOptions:
    output_directory: Path
    write_npz: bool
    write_json_provenance: bool
    allow_combined_display: bool

@dataclass(frozen=True)
class FFTInputContract:
    sequence_metadata_path: Path | None
    expected_sample_rate_hz: float | None
    required_unit: str
    accepted_channels: tuple[str, ...]

@dataclass(frozen=True)
class FFTExtensionConfig:
    input_contract: FFTInputContract
    analysis: FFTAnalysisOptions
    artifacts: FFTArtifactOptions

def _positive(value: Any, label: str) -> float:
    value = float(value)
    if not value > 0.0:
        raise FFTConfigurationError(f"{label} must be positive")
    return value

def load_fft_extension_config(path: str | Path) -> FFTExtensionConfig:
    """Load only configuration; it never opens a signal sequence."""
    path = Path(path).expanduser().resolve()
    with path.open() as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise FFTConfigurationError("FFT configuration must be a YAML mapping")
    try:
        i, a, o = raw["input_contract"], raw["analysis"], raw["artifacts"]
        expected_fs = i.get("expected_sample_rate_hz")
        config = FFTExtensionConfig(
            input_contract=FFTInputContract(
                sequence_metadata_path=(None if i.get("sequence_metadata_path") in (None, "") else Path(i["sequence_metadata_path"]).expanduser()),
                expected_sample_rate_hz=(None if expected_fs in (None, "") else _positive(expected_fs, "expected_sample_rate_hz")),
                required_unit=str(i["required_unit"]),
                accepted_channels=tuple(str(x) for x in i["accepted_channels"]),
            ),
            analysis=FFTAnalysisOptions(
                requested_display_min_hz=_positive(a["requested_display_range_hz"][0], "requested display minimum"),
                requested_display_max_hz=_positive(a["requested_display_range_hz"][1], "requested display maximum"),
                fft_min_hz=_positive(a["fft_range_hz"][0], "FFT minimum"),
                fft_max_hz=_positive(a["fft_range_hz"][1], "FFT maximum"),
                c_lpsd_handoff_hz=_positive(a["c_lpsd_handoff_hz"], "C-LPSD handoff"),
                c_lpsd_min_hz=_positive(a["c_lpsd_calculable_range_hz"][0], "C-LPSD minimum"),
                c_lpsd_max_hz=_positive(a["c_lpsd_calculable_range_hz"][1], "C-LPSD maximum"),
                detrend=str(a["detrend"]), window_policy=str(a["window_policy"]),
                kaiser_beta=(None if a.get("kaiser_beta") is None else _positive(a["kaiser_beta"], "kaiser_beta")),
                one_sided=bool(a["one_sided"]), full_record_only=bool(a["full_record_only"]),
                zero_padding_factor=int(a["zero_padding_factor"]), memory_limit_bytes=int(a["memory_limit_bytes"]),
            ),
            artifacts=FFTArtifactOptions(
                output_directory=Path(o["output_directory"]).expanduser(), write_npz=bool(o["write_npz"]),
                write_json_provenance=bool(o["write_json_provenance"]), allow_combined_display=bool(o["allow_combined_display"]),
            ),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise FFTConfigurationError(f"invalid FFT extension configuration: {error}") from error
    a = config.analysis
    if not (a.requested_display_min_hz < a.requested_display_max_hz and a.fft_min_hz < a.fft_max_hz):
        raise FFTConfigurationError("all configured ranges must be increasing")
    if not (a.c_lpsd_min_hz <= a.c_lpsd_handoff_hz <= a.c_lpsd_max_hz):
        raise FFTConfigurationError("handoff must lie within configured C-LPSD calculable range")
    if a.detrend not in {"none", "mean", "linear"}:
        raise FFTConfigurationError("detrend must be none, mean, or linear")
    if a.window_policy not in {"rectangular", "numpy_kaiser_to_be_validated"}:
        raise FFTConfigurationError("unsupported window_policy")
    if a.window_policy == "numpy_kaiser_to_be_validated" and a.kaiser_beta is None:
        raise FFTConfigurationError("numpy_kaiser_to_be_validated requires kaiser_beta")
    if not a.one_sided or not a.full_record_only or a.zero_padding_factor != 1:
        raise FFTConfigurationError("approved extension requires one-sided, full-record, unpadded FFT")
    if a.memory_limit_bytes <= 0 or not config.input_contract.accepted_channels:
        raise FFTConfigurationError("memory limit and accepted channels are required")
    if config.input_contract.required_unit != "V":
        raise FFTConfigurationError("approved FFT extension is voltage-only: required_unit must be V")
    return config
