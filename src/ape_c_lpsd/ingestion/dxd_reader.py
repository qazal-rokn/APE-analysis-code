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
from typing import Iterator, Mapping

import numpy as np

from .models import APEDXDChunk, APEDXDData, DXDInspection


class DXDReaderError(RuntimeError):
    """The requested .dxd data cannot be represented without hidden changes."""


_REQUIRED_KEYS = (
    "x_1f",
    "y_1f",
    "x_2f",
    "y_2f",
    "cavity_transmission",
)

_DEFAULT_CHUNK_SAMPLES = 1_000_000


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


def _timestamps_match_declared_rate(
    timestamps_s: np.ndarray, declared_sample_rate_hz: float
) -> bool:
    """Return True only when all timestamp steps agree with the declared rate."""
    if timestamps_s.size < 2:
        return False
    if not np.isfinite(declared_sample_rate_hz) or declared_sample_rate_hz <= 0:
        return False
    dt = np.diff(timestamps_s)
    if not np.all(np.isfinite(dt)) or np.any(dt <= 0):
        return False
    expected_step = 1.0 / declared_sample_rate_hz
    tolerance = max(1e-12, abs(expected_step) * 1e-6)
    return bool(np.allclose(dt, expected_step, rtol=0.0, atol=tolerance))


def _assemble_ape_data(
    *,
    source_file: Path,
    arrays: Mapping[str, np.ndarray],
    timestamps: Mapping[str, np.ndarray],
    units: Mapping[str, str],
    channel_names: Mapping[str, str],
    start_time: datetime | None,
    metadata: dict,
    declared_sample_rate_hz: float | None = None,
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

    timestamp_rate_hz = _regular_rate(reference_t)
    regular_sampling = timestamp_rate_hz is not None
    relative_timestamps_s = None

    if declared_sample_rate_hz is not None:
        if not np.isfinite(declared_sample_rate_hz) or declared_sample_rate_hz <= 0:
            raise DXDReaderError(
                f"invalid Dewesoft declared sample rate: {declared_sample_rate_hz}"
            )
        if not _timestamps_match_declared_rate(reference_t, declared_sample_rate_hz):
            raise DXDReaderError(
                "selected-channel timestamps are inconsistent with the Dewesoft "
                f"declared sample rate {declared_sample_rate_hz} Hz"
            )
        sample_rate_hz = float(declared_sample_rate_hz)
        regular_sampling = True
    else:
        sample_rate_hz = timestamp_rate_hz
        if not regular_sampling:
            # Without authoritative file-rate metadata, preserve actual timing
            # rather than pretending the record is uniform.
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



def _validate_channel_mapping(channel_names: Mapping[str, str]) -> list[str]:
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
    return requested_names


def iter_ape_dxd_chunks(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_samples: int = _DEFAULT_CHUNK_SAMPLES,
) -> Iterator[APEDXDChunk]:
    """Yield synchronized raw chunks without retaining the full recording in RAM."""
    if chunk_samples <= 1:
        raise DXDReaderError("chunk_samples must be greater than 1")

    source = Path(input_file).expanduser().resolve()
    if not source.is_file():
        raise DXDReaderError(f"input .dxd file not found: {source}")
    if source.suffix.lower() != ".dxd":
        raise DXDReaderError(f"expected a .dxd input file, got: {source.name}")

    requested_names = _validate_channel_mapping(channel_names)
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

        declared_rate = float(file_info.sample_rate)
        if not np.isfinite(declared_rate) or declared_rate <= 0:
            raise DXDReaderError(f"invalid Dewesoft declared sample rate: {declared_rate}")

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

        channel_indices: dict[str, c_int] = {}
        channel_units: dict[str, str] = {}
        counts: dict[str, int] = {}
        for role in _REQUIRED_KEYS:
            recorded_name = str(channel_names[role])
            channel = by_name[recorded_name]
            if int(channel.array_size) != 1:
                raise DXDReaderError(
                    f"{recorded_name} has array_size={channel.array_size}; only scalar channels are supported"
                )
            index = c_int(int(channel.index))
            count = int(dll.DWGetScaledSamplesCount(index))
            if count <= 0:
                raise DXDReaderError(
                    f"invalid sample count for {recorded_name}: {count}"
                )
            channel_indices[role] = index
            channel_units[role] = _decode_c_string(channel.unit)
            counts[role] = count

        if len(set(counts.values())) != 1:
            raise DXDReaderError(f"selected channel lengths disagree: {counts}")
        total_count = next(iter(counts.values()))

        previous_last_timestamp: float | None = None
        expected_step = 1.0 / declared_rate
        boundary_tolerance = max(1e-12, abs(expected_step) * 1e-6)

        for raw_start in range(0, total_count, chunk_samples):
            raw_stop = min(raw_start + chunk_samples, total_count)
            count = raw_stop - raw_start
            arrays: dict[str, np.ndarray] = {}
            reference_timestamps: np.ndarray | None = None

            for role in _REQUIRED_KEYS:
                data_buffer = (c_double * count)()
                time_buffer = (c_double * count)()
                status = int(
                    dll.DWGetScaledSamples(
                        channel_indices[role],
                        c_int64(raw_start),
                        count,
                        data_buffer,
                        time_buffer,
                    )
                )
                if status != 0:
                    raise DXDReaderError(
                        f"DWGetScaledSamples() failed for {channel_names[role]} "
                        f"at raw samples [{raw_start}, {raw_stop})"
                    )

                values = np.ctypeslib.as_array(data_buffer).copy()
                times = np.ctypeslib.as_array(time_buffer).copy()
                if not np.all(np.isfinite(values)):
                    raise DXDReaderError(
                        f"{role} contains non-finite raw samples in "
                        f"[{raw_start}, {raw_stop})"
                    )
                if not np.all(np.isfinite(times)):
                    raise DXDReaderError(
                        f"{role} contains non-finite timestamps in "
                        f"[{raw_start}, {raw_stop})"
                    )

                if reference_timestamps is None:
                    reference_timestamps = times
                elif not np.allclose(times, reference_timestamps, rtol=0.0, atol=1e-9):
                    raise DXDReaderError(
                        f"selected channels are not synchronized in "
                        f"[{raw_start}, {raw_stop})"
                    )
                arrays[role] = values

            assert reference_timestamps is not None
            if not _timestamps_match_declared_rate(reference_timestamps, declared_rate):
                raise DXDReaderError(
                    "selected-channel timestamps are inconsistent with the Dewesoft "
                    f"declared sample rate {declared_rate} Hz in "
                    f"[{raw_start}, {raw_stop})"
                )

            first_timestamp = float(reference_timestamps[0])
            last_timestamp = float(reference_timestamps[-1])
            if previous_last_timestamp is not None:
                observed_boundary_step = first_timestamp - previous_last_timestamp
                if not np.isclose(
                    observed_boundary_step,
                    expected_step,
                    rtol=0.0,
                    atol=boundary_tolerance,
                ):
                    raise DXDReaderError(
                        "timestamp discontinuity between chunks at raw sample "
                        f"{raw_start}: observed step={observed_boundary_step}, "
                        f"expected={expected_step}"
                    )
            previous_last_timestamp = last_timestamp

            yield APEDXDChunk(
                source_file=source,
                raw_start=raw_start,
                raw_stop=raw_stop,
                sample_rate_hz=declared_rate,
                x_1f=arrays["x_1f"],
                y_1f=arrays["y_1f"],
                x_2f=arrays["x_2f"],
                y_2f=arrays["y_2f"],
                cavity_transmission=arrays["cavity_transmission"],
                channel_names={key: str(channel_names[key]) for key in _REQUIRED_KEYS},
                channel_units=dict(channel_units),
                first_timestamp_s=first_timestamp,
                last_timestamp_s=last_timestamp,
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


def inspect_ape_dxd(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_samples: int = _DEFAULT_CHUNK_SAMPLES,
) -> DXDInspection:
    """Validate an entire .dxd file while retaining only one chunk at a time."""
    source = Path(input_file).expanduser().resolve()
    total = 0
    chunk_count = 0
    sample_rate_hz: float | None = None
    units: dict[str, str] | None = None
    first_timestamp_s: float | None = None
    last_timestamp_s: float | None = None

    for chunk in iter_ape_dxd_chunks(
        source,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
        chunk_samples=chunk_samples,
    ):
        if chunk_count == 0:
            sample_rate_hz = chunk.sample_rate_hz
            units = dict(chunk.channel_units)
            first_timestamp_s = chunk.first_timestamp_s
        total += chunk.raw_stop - chunk.raw_start
        chunk_count += 1
        last_timestamp_s = chunk.last_timestamp_s

    if chunk_count == 0 or sample_rate_hz is None or units is None:
        raise DXDReaderError("Dewesoft file produced no chunks")

    duration_s = total / sample_rate_hz
    start_time = None
    # Read only file-level metadata through the first chunk's validated source.
    # The public inspection contract reports duration from N/fs so it is tied
    # to the exact selected-channel sample count.
    return DXDInspection(
        source_file=source,
        sample_rate_hz=sample_rate_hz,
        n_samples=total,
        duration_s=duration_s,
        start_time=start_time,
        channel_names={key: str(channel_names[key]) for key in _REQUIRED_KEYS},
        channel_units=units,
        chunk_count=chunk_count,
        chunk_samples=chunk_samples,
        metadata={
            "reader": "DWDataReader",
            "first_selected_timestamp_s": first_timestamp_s,
            "last_selected_timestamp_s": last_timestamp_s,
            "full_record_retained_in_memory": False,
        },
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
            if count == 0:
                raise DXDReaderError(f"{recorded_name} contains no samples")

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
            declared_sample_rate_hz=float(file_info.sample_rate),
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
