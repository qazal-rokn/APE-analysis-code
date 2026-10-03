"""Canonical APE Dewesoft .dxd ingestion.

This module reads exactly the five raw APE channels needed by the first
analysis stage. It performs no calibration, filtering, cleaning, resampling,
interpolation, trimming, or gap filling.

The implementation calls Dewesoft's vendored DWDataReader library through the
official Python header shipped with the reference APE project.
"""
from __future__ import annotations

from ctypes import byref, c_char_p, c_double, c_int, c_int64
from datetime import datetime, timedelta
import importlib.util
from pathlib import Path
from typing import Mapping

import numpy as np

from .models import APEDXDData


class DXDReaderError(RuntimeError):
    """The requested .dxd data cannot be represented without hidden changes."""


_REQUIRED_KEYS = (
    "x_1f",
    "y_1f",
    "x_2f",
    "y_2f",
    "cavity_transmission",
)


def _decode_c_string(value) -> str:
    raw = bytes(value)
    return raw.split(b"\x00", 1)[0].decode("utf-8", errors="replace")


def _load_header(wrapper_path: str | Path):
    path = Path(wrapper_path).expanduser().resolve()
    if not path.is_file():
        raise DXDReaderError(f"Dewesoft Python header not found: {path}")
    spec = importlib.util.spec_from_file_location("ape_dewesoft_header", path)
    if spec is None or spec.loader is None:
        raise DXDReaderError(f"cannot import Dewesoft Python header: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _ole_datetime(days: float) -> datetime | None:
    if not np.isfinite(days) or days <= 0:
        return None
    # Dewesoft DWFileInfo.start_store_time uses the OLE Automation epoch.
    # The SDK value does not encode a timezone; do not invent one here.
    return datetime(1899, 12, 30) + timedelta(days=float(days))


def _regular_rate(timestamps_s: np.ndarray) -> float | None:
    if timestamps_s.size < 2:
        return None
    dt = np.diff(timestamps_s)
    if not np.all(np.isfinite(dt)) or np.any(dt <= 0):
        return None
    step = float(np.median(dt))
    tolerance = max(1e-12, abs(step) * 1e-6)
    if not np.allclose(dt, step, rtol=0.0, atol=tolerance):
        return None
    return 1.0 / step


def _assemble_ape_data(
    *,
    source_file: Path,
    arrays: Mapping[str, np.ndarray],
    timestamps: Mapping[str, np.ndarray],
    units: Mapping[str, str],
    channel_names: Mapping[str, str],
    start_time: datetime | None,
    metadata: dict,
) -> APEDXDData:
    """Validate synchronized raw channels and construct the canonical object."""
    missing = [key for key in _REQUIRED_KEYS if key not in arrays or key not in timestamps]
    if missing:
        raise DXDReaderError(f"missing required APE channel data: {missing}")

    values = {key: np.asarray(arrays[key], dtype=np.float64) for key in _REQUIRED_KEYS}
    times = {key: np.asarray(timestamps[key], dtype=np.float64) for key in _REQUIRED_KEYS}

    for key in _REQUIRED_KEYS:
        if values[key].ndim != 1 or times[key].ndim != 1:
            raise DXDReaderError(f"{key} must be one-dimensional")
        if values[key].size != times[key].size:
            raise DXDReaderError(f"{key} data/timestamp length mismatch")
        if not np.all(np.isfinite(values[key])):
            raise DXDReaderError(f"{key} contains non-finite raw samples")
        if not np.all(np.isfinite(times[key])):
            raise DXDReaderError(f"{key} contains non-finite timestamps")

    lengths = {key: values[key].size for key in _REQUIRED_KEYS}
    if len(set(lengths.values())) != 1:
        raise DXDReaderError(f"selected channel lengths disagree: {lengths}")
    n_samples = next(iter(lengths.values()))
    if n_samples == 0:
        raise DXDReaderError("selected channels contain no samples")

    reference_key = _REQUIRED_KEYS[0]
    reference_t = times[reference_key]
    for key in _REQUIRED_KEYS[1:]:
        candidate = times[key]
        if candidate.shape != reference_t.shape or not np.allclose(
            candidate, reference_t, rtol=0.0, atol=1e-9
        ):
            raise DXDReaderError(
                f"selected channels are not synchronized: {reference_key} and {key} timestamps differ"
            )

    sample_rate_hz = _regular_rate(reference_t)
    regular_sampling = sample_rate_hz is not None
    relative_timestamps_s = None
    if not regular_sampling:
        # Preserve the actual timing instead of pretending the record is uniform.
        relative_timestamps_s = reference_t - reference_t[0]

    return APEDXDData(
        source_file=source_file,
        sample_rate_hz=sample_rate_hz,
        n_samples=n_samples,
        start_time=start_time,
        x_1f=values["x_1f"],
        y_1f=values["y_1f"],
        x_2f=values["x_2f"],
        y_2f=values["y_2f"],
        cavity_transmission=values["cavity_transmission"],
        channel_names={key: str(channel_names[key]) for key in _REQUIRED_KEYS},
        channel_units={key: str(units.get(key, "")) for key in _REQUIRED_KEYS},
        regular_sampling=regular_sampling,
        relative_timestamps_s=relative_timestamps_s,
        metadata=dict(metadata),
    )


def read_ape_dxd(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
) -> APEDXDData:
    """Read the five configured raw APE channels from one Dewesoft .dxd file.

    Parameters are intentionally explicit. The caller supplies the physical
    channel mapping; this function never guesses channel names.
    """
    source = Path(input_file).expanduser().resolve()
    if not source.is_file():
        raise DXDReaderError(f"input .dxd file not found: {source}")
    if source.suffix.lower() != ".dxd":
        raise DXDReaderError(f"expected a .dxd input file, got: {source.name}")

    keys = set(channel_names)
    missing = [key for key in _REQUIRED_KEYS if key not in keys]
    extra = sorted(keys.difference(_REQUIRED_KEYS))
    if missing or extra:
        raise DXDReaderError(
            f"channel mapping must contain exactly {_REQUIRED_KEYS}; missing={missing}, extra={extra}"
        )
    requested_names = [str(channel_names[key]) for key in _REQUIRED_KEYS]
    if len(set(requested_names)) != len(requested_names):
        raise DXDReaderError("the five APE roles must map to five distinct Dewesoft channels")

    header = _load_header(wrapper_path)
    library = Path(library_path).expanduser().resolve()
    if not library.is_file():
        raise DXDReaderError(f"Dewesoft native library not found: {library}")

    import ctypes

    dll = ctypes.cdll.LoadLibrary(str(library))
    initialized = False
    opened = False
    try:
        if int(dll.DWInit()) != 0:
            raise DXDReaderError("Dewesoft DWInit() failed")
        initialized = True
        if int(dll.DWAddReader()) != 0:
            raise DXDReaderError("Dewesoft DWAddReader() failed")

        file_info = header.DWFileInfo(0, 0, 0)
        if int(dll.DWOpenDataFile(c_char_p(str(source).encode()), byref(file_info))) != 0:
            raise DXDReaderError(f"Dewesoft could not open: {source}")
        opened = True

        channel_count = int(dll.DWGetChannelListCount())
        if channel_count < 0:
            raise DXDReaderError("Dewesoft DWGetChannelListCount() failed")
        channel_list = (header.DWChannel * channel_count)()
        if int(dll.DWGetChannelList(byref(channel_list))) != 0:
            raise DXDReaderError("Dewesoft DWGetChannelList() failed")

        by_name: dict[str, object] = {}
        duplicates: set[str] = set()
        for channel in channel_list:
            name = _decode_c_string(channel.name)
            if name in by_name:
                duplicates.add(name)
            by_name[name] = channel
        requested_duplicates = sorted(set(requested_names).intersection(duplicates))
        if requested_duplicates:
            raise DXDReaderError(
                f"requested Dewesoft channel names are duplicated in the file: {requested_duplicates}"
            )

        absent = [name for name in requested_names if name not in by_name]
        if absent:
            raise DXDReaderError(
                f"requested Dewesoft channels are absent: {absent}; "
                f"available channels={sorted(by_name)}"
            )

        arrays: dict[str, np.ndarray] = {}
        timestamps: dict[str, np.ndarray] = {}
        units: dict[str, str] = {}
        first_sdk_timestamp_s: dict[str, float] = {}

        for role in _REQUIRED_KEYS:
            recorded_name = str(channel_names[role])
            channel = by_name[recorded_name]
            if int(channel.array_size) != 1:
                raise DXDReaderError(
                    f"{recorded_name} has array_size={channel.array_size}; only scalar channels are supported"
                )

            channel_index = c_int(int(channel.index))
            count = int(dll.DWGetScaledSamplesCount(channel_index))
            if count < 0:
                raise DXDReaderError(f"DWGetScaledSamplesCount() failed for {recorded_name}")

            data_buffer = (c_double * count)()
            time_buffer = (c_double * count)()
            status = int(
                dll.DWGetScaledSamples(
                    channel_index,
                    c_int64(0),
                    count,
                    data_buffer,
                    time_buffer,
                )
            )
            if status != 0:
                raise DXDReaderError(f"DWGetScaledSamples() failed for {recorded_name}")

            arrays[role] = np.ctypeslib.as_array(data_buffer).copy()
            timestamps[role] = np.ctypeslib.as_array(time_buffer).copy()
            units[role] = _decode_c_string(channel.unit)
            first_sdk_timestamp_s[role] = (
                float(timestamps[role][0]) if timestamps[role].size else float("nan")
            )

        metadata = {
            "dewesoft_file_sample_rate_hz": float(file_info.sample_rate),
            "dewesoft_duration_s": float(file_info.duration),
            "dewesoft_start_store_time_ole_days": float(file_info.start_store_time),
            "dewesoft_first_sdk_timestamp_s": first_sdk_timestamp_s,
            "dewesoft_channel_count": channel_count,
            "reader": "DWDataReader",
        }
        return _assemble_ape_data(
            source_file=source,
            arrays=arrays,
            timestamps=timestamps,
            units=units,
            channel_names=channel_names,
            start_time=_ole_datetime(float(file_info.start_store_time)),
            metadata=metadata,
        )
    finally:
        if opened:
            close_status = int(dll.DWCloseDataFile())
            if close_status != 0:
                raise DXDReaderError("Dewesoft DWCloseDataFile() failed")
        if initialized:
            deinit_status = int(dll.DWDeInit())
            if deinit_status != 0:
                raise DXDReaderError("Dewesoft DWDeInit() failed")
