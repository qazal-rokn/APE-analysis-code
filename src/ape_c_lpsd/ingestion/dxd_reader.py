"""Canonical APE Dewesoft .dxd ingestion.

The production path is streaming: metadata are inspected once and raw data are
read in bounded synchronized chunks. No calibration, filtering, cleaning,
resampling, interpolation, trimming, or gap filling is performed here.
"""
from __future__ import annotations

from ctypes import byref, c_char_p, c_double, c_int, c_int64
from datetime import datetime, timedelta
import importlib.util
from pathlib import Path
from typing import Iterator, Mapping

import numpy as np

from .models import APEDXDChunk, APEDXDData, APEDXDMetadata


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
    if timestamps_s.size < 2:
        return True
    if not np.isfinite(declared_sample_rate_hz) or declared_sample_rate_hz <= 0:
        return False
    dt = np.diff(timestamps_s)
    if not np.all(np.isfinite(dt)) or np.any(dt <= 0):
        return False
    expected_step = 1.0 / declared_sample_rate_hz
    tolerance = max(1e-12, abs(expected_step) * 1e-6)
    return bool(np.allclose(dt, expected_step, rtol=0.0, atol=tolerance))


def _validate_mapping(channel_names: Mapping[str, str]) -> list[str]:
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


def _validate_source(input_file: str | Path) -> Path:
    source = Path(input_file).expanduser().resolve()
    if not source.is_file():
        raise DXDReaderError(f"input .dxd file not found: {source}")
    if source.suffix.lower() != ".dxd":
        raise DXDReaderError(f"expected a .dxd input file, got: {source.name}")
    return source


def _open_context(
    source: Path, wrapper_path: str | Path, library_path: str | Path
):
    """Open DWDataReader and return (header, dll, file_info, channel_list)."""
    import ctypes

    header = _load_header(wrapper_path)
    library = Path(library_path).expanduser().resolve()
    if not library.is_file():
        raise DXDReaderError(f"Dewesoft native library not found: {library}")
    dll = ctypes.cdll.LoadLibrary(str(library))
    if int(dll.DWInit()) != 0:
        raise DXDReaderError("Dewesoft DWInit() failed")
    try:
        if int(dll.DWAddReader()) != 0:
            raise DXDReaderError("Dewesoft DWAddReader() failed")
        file_info = header.DWFileInfo(0, 0, 0)
        if int(dll.DWOpenDataFile(c_char_p(str(source).encode()), byref(file_info))) != 0:
            raise DXDReaderError(f"Dewesoft could not open: {source}")
        channel_count = int(dll.DWGetChannelListCount())
        if channel_count < 0:
            raise DXDReaderError("Dewesoft DWGetChannelListCount() failed")
        channel_list = (header.DWChannel * channel_count)()
        if int(dll.DWGetChannelList(byref(channel_list))) != 0:
            raise DXDReaderError("Dewesoft DWGetChannelList() failed")
        return header, dll, file_info, channel_list
    except Exception:
        dll.DWDeInit()
        raise


def _close_context(dll) -> None:
    close_status = int(dll.DWCloseDataFile())
    deinit_status = int(dll.DWDeInit())
    if close_status != 0:
        raise DXDReaderError("Dewesoft DWCloseDataFile() failed")
    if deinit_status != 0:
        raise DXDReaderError("Dewesoft DWDeInit() failed")


def _selected_channels(channel_list, channel_names: Mapping[str, str]):
    requested_names = _validate_mapping(channel_names)
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
            f"requested Dewesoft channels are absent: {absent}; available channels={sorted(by_name)}"
        )
    selected = {role: by_name[str(channel_names[role])] for role in _REQUIRED_KEYS}
    for role, channel in selected.items():
        if int(channel.array_size) != 1:
            raise DXDReaderError(
                f"{channel_names[role]} has array_size={channel.array_size}; only scalar channels are supported"
            )
    return selected


def inspect_ape_dxd(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
) -> APEDXDMetadata:
    """Inspect file/channel metadata without allocating full channel arrays."""
    source = _validate_source(input_file)
    _validate_mapping(channel_names)
    header, dll, file_info, channel_list = _open_context(source, wrapper_path, library_path)
    try:
        selected = _selected_channels(channel_list, channel_names)
        counts: dict[str, int] = {}
        units: dict[str, str] = {}
        for role, channel in selected.items():
            count = int(dll.DWGetScaledSamplesCount(c_int(int(channel.index))))
            if count <= 0:
                raise DXDReaderError(
                    f"DWGetScaledSamplesCount() failed or returned no samples for {channel_names[role]}"
                )
            counts[role] = count
            units[role] = _decode_c_string(channel.unit)
        if len(set(counts.values())) != 1:
            raise DXDReaderError(f"selected channel lengths disagree: {counts}")
        fs = float(file_info.sample_rate)
        if not np.isfinite(fs) or fs <= 0:
            raise DXDReaderError(f"invalid Dewesoft declared sample rate: {fs}")
        n_samples = next(iter(counts.values()))
        metadata = {
            "dewesoft_file_sample_rate_hz": fs,
            "dewesoft_duration_s": float(file_info.duration),
            "dewesoft_start_store_time_ole_days": float(file_info.start_store_time),
            "dewesoft_channel_count": len(channel_list),
            "reader": "DWDataReader",
            "read_mode": "streaming",
        }
        return APEDXDMetadata(
            source_file=source,
            sample_rate_hz=fs,
            n_samples=n_samples,
            duration_s=float(file_info.duration),
            start_time=_ole_datetime(float(file_info.start_store_time)),
            channel_names={key: str(channel_names[key]) for key in _REQUIRED_KEYS},
            channel_units=units,
            dewesoft_channel_count=len(channel_list),
            metadata=metadata,
        )
    finally:
        _close_context(dll)


def _validate_chunk_arrays(
    *,
    arrays: Mapping[str, np.ndarray],
    timestamps: Mapping[str, np.ndarray],
    sample_rate_hz: float,
) -> tuple[float, float]:
    lengths = {key: len(arrays[key]) for key in _REQUIRED_KEYS}
    if len(set(lengths.values())) != 1 or next(iter(lengths.values())) <= 0:
        raise DXDReaderError(f"chunk channel lengths disagree: {lengths}")
    reference = np.asarray(timestamps[_REQUIRED_KEYS[0]], dtype=np.float64)
    if not np.all(np.isfinite(reference)):
        raise DXDReaderError("chunk contains non-finite timestamps")
    for role in _REQUIRED_KEYS:
        values = np.asarray(arrays[role], dtype=np.float64)
        times = np.asarray(timestamps[role], dtype=np.float64)
        if not np.all(np.isfinite(values)):
            raise DXDReaderError(f"{role} contains non-finite raw samples")
        if times.shape != reference.shape or not np.allclose(
            times, reference, rtol=0.0, atol=1e-9
        ):
            raise DXDReaderError(
                f"selected channels are not synchronized within chunk: {_REQUIRED_KEYS[0]} vs {role}"
            )
    if not _timestamps_match_declared_rate(reference, sample_rate_hz):
        raise DXDReaderError(
            f"chunk timestamps are inconsistent with Dewesoft declared sample rate {sample_rate_hz} Hz"
        )
    return float(reference[0]), float(reference[-1])


def _validate_chunk_boundary(
    *, previous_last_timestamp_s: float, current_first_timestamp_s: float, sample_rate_hz: float
) -> None:
    expected = previous_last_timestamp_s + 1.0 / sample_rate_hz
    tolerance = max(1e-12, (1.0 / sample_rate_hz) * 1e-6)
    if not np.isclose(current_first_timestamp_s, expected, rtol=0.0, atol=tolerance):
        raise DXDReaderError(
            "timestamp discontinuity between chunks: "
            f"expected {expected}, observed {current_first_timestamp_s}"
        )


def iter_ape_dxd_chunks(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_samples: int = 1_000_000,
) -> Iterator[APEDXDChunk]:
    """Yield validated synchronized chunks without retaining the full file in RAM."""
    if chunk_samples <= 0:
        raise DXDReaderError("chunk_samples must be positive")
    metadata = inspect_ape_dxd(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    )
    source = metadata.source_file
    header, dll, file_info, channel_list = _open_context(source, wrapper_path, library_path)
    del header, file_info
    previous_last: float | None = None
    try:
        selected = _selected_channels(channel_list, channel_names)
        for start in range(0, metadata.n_samples, chunk_samples):
            count = min(chunk_samples, metadata.n_samples - start)
            arrays: dict[str, np.ndarray] = {}
            timestamps: dict[str, np.ndarray] = {}
            for role, channel in selected.items():
                data_buffer = (c_double * count)()
                time_buffer = (c_double * count)()
                status = int(
                    dll.DWGetScaledSamples(
                        c_int(int(channel.index)),
                        c_int64(start),
                        count,
                        data_buffer,
                        time_buffer,
                    )
                )
                if status != 0:
                    raise DXDReaderError(
                        f"DWGetScaledSamples() failed for {channel_names[role]} "
                        f"at start={start}, count={count}"
                    )
                arrays[role] = np.ctypeslib.as_array(data_buffer).copy()
                timestamps[role] = np.ctypeslib.as_array(time_buffer).copy()
            first_ts, last_ts = _validate_chunk_arrays(
                arrays=arrays,
                timestamps=timestamps,
                sample_rate_hz=metadata.sample_rate_hz,
            )
            if previous_last is not None:
                _validate_chunk_boundary(
                    previous_last_timestamp_s=previous_last,
                    current_first_timestamp_s=first_ts,
                    sample_rate_hz=metadata.sample_rate_hz,
                )
            previous_last = last_ts
            yield APEDXDChunk(
                source_file=source,
                start_index=start,
                stop_index=start + count,
                sample_rate_hz=metadata.sample_rate_hz,
                x_1f=arrays["x_1f"],
                y_1f=arrays["y_1f"],
                x_2f=arrays["x_2f"],
                y_2f=arrays["y_2f"],
                cavity_transmission=arrays["cavity_transmission"],
                first_timestamp_s=first_ts,
                last_timestamp_s=last_ts,
            )
    finally:
        _close_context(dll)


def validate_ape_dxd_stream(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_samples: int = 1_000_000,
) -> APEDXDMetadata:
    """Read the entire file chunk-by-chunk, validating without retaining it."""
    metadata = inspect_ape_dxd(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    )
    observed = 0
    for chunk in iter_ape_dxd_chunks(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
        chunk_samples=chunk_samples,
    ):
        if chunk.start_index != observed:
            raise DXDReaderError(
                f"stream index discontinuity: expected {observed}, got {chunk.start_index}"
            )
        observed = chunk.stop_index
    if observed != metadata.n_samples:
        raise DXDReaderError(
            f"stream sample-count mismatch: expected {metadata.n_samples}, observed {observed}"
        )
    return metadata


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
    """Compatibility helper used by controlled small-array tests."""
    missing = [key for key in _REQUIRED_KEYS if key not in arrays or key not in timestamps]
    if missing:
        raise DXDReaderError(f"missing required APE channel data: {missing}")
    values = {key: np.asarray(arrays[key], dtype=np.float64) for key in _REQUIRED_KEYS}
    times = {key: np.asarray(timestamps[key], dtype=np.float64) for key in _REQUIRED_KEYS}
    lengths = {key: values[key].size for key in _REQUIRED_KEYS}
    if len(set(lengths.values())) != 1:
        raise DXDReaderError(f"selected channel lengths disagree: {lengths}")
    n_samples = next(iter(lengths.values()))
    if n_samples <= 0:
        raise DXDReaderError("selected channels contain no samples")
    first_ts, _ = _validate_chunk_arrays(
        arrays=values,
        timestamps=times,
        sample_rate_hz=declared_sample_rate_hz or _regular_rate(times[_REQUIRED_KEYS[0]]) or -1.0,
    )
    del first_ts
    if declared_sample_rate_hz is not None:
        sample_rate_hz = float(declared_sample_rate_hz)
        regular_sampling = True
        relative_timestamps_s = None
    else:
        sample_rate_hz = _regular_rate(times[_REQUIRED_KEYS[0]])
        regular_sampling = sample_rate_hz is not None
        relative_timestamps_s = None if regular_sampling else (
            times[_REQUIRED_KEYS[0]] - times[_REQUIRED_KEYS[0]][0]
        )
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
    """Compatibility full-memory reader; do not use for production-size files."""
    metadata = inspect_ape_dxd(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    )
    arrays = {key: [] for key in _REQUIRED_KEYS}
    timestamps = []
    for chunk in iter_ape_dxd_chunks(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    ):
        arrays["x_1f"].append(chunk.x_1f)
        arrays["y_1f"].append(chunk.y_1f)
        arrays["x_2f"].append(chunk.x_2f)
        arrays["y_2f"].append(chunk.y_2f)
        arrays["cavity_transmission"].append(chunk.cavity_transmission)
        timestamps.append(
            chunk.first_timestamp_s
            + np.arange(chunk.stop_index - chunk.start_index, dtype=np.float64)
            / metadata.sample_rate_hz
        )
    merged = {key: np.concatenate(parts) for key, parts in arrays.items()}
    merged_t = np.concatenate(timestamps)
    time_map = {key: merged_t for key in _REQUIRED_KEYS}
    return _assemble_ape_data(
        source_file=metadata.source_file,
        arrays=merged,
        timestamps=time_map,
        units=metadata.channel_units,
        channel_names=metadata.channel_names,
        start_time=metadata.start_time,
        metadata=metadata.metadata,
        declared_sample_rate_hz=metadata.sample_rate_hz,
    )
