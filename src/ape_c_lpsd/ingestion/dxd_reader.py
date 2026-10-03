"""Canonical APE Dewesoft .dxd ingestion.

This module reads exactly the five raw APE channels needed by the first
analysis stage. It performs no calibration, filtering, cleaning, resampling,
interpolation, trimming, or gap filling.

Large files are processed in bounded chunks. The lightweight inspection path
validates the entire file without retaining the full recording in memory.
"""
from __future__ import annotations

from ctypes import byref, c_char_p, c_double, c_int, c_int64
from datetime import datetime, timedelta
import importlib.util
from pathlib import Path
from typing import Iterator, Mapping

import numpy as np

from .models import APEDXDChunk, APEDXDData, APEDXDInspection


class DXDReaderError(RuntimeError):
    """The requested .dxd data cannot be represented without hidden changes."""


_REQUIRED_KEYS = (
    "x_1f",
    "y_1f",
    "x_2f",
    "y_2f",
    "cavity_transmission",
)

_DEFAULT_CHUNK_SIZE_SAMPLES = 1_000_000


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
    # The SDK value itself does not encode a timezone, so we do not invent one.
    return datetime(1899, 12, 30) + timedelta(days=float(days))


def _timestamp_tolerance(sample_rate_hz: float) -> float:
    step = 1.0 / sample_rate_hz
    return max(1e-12, abs(step) * 1e-6)


def _validate_declared_rate(timestamps_s: np.ndarray, sample_rate_hz: float) -> None:
    """Validate timestamp steps against Dewesoft's declared sample rate."""
    if timestamps_s.size < 2:
        return
    dt = np.diff(timestamps_s)
    expected = 1.0 / sample_rate_hz
    tol = _timestamp_tolerance(sample_rate_hz)
    if (
        not np.all(np.isfinite(dt))
        or np.any(dt <= 0)
        or not np.allclose(dt, expected, rtol=0.0, atol=tol)
    ):
        raise DXDReaderError(
            "selected-channel timestamps are inconsistent with the Dewesoft "
            f"declared sample rate {sample_rate_hz} Hz"
        )


def _validate_chunk_boundary(
    previous_last_timestamp_s: float | None,
    current_first_timestamp_s: float,
    sample_rate_hz: float,
) -> None:
    if previous_last_timestamp_s is None:
        return
    expected = previous_last_timestamp_s + 1.0 / sample_rate_hz
    if not np.isclose(
        current_first_timestamp_s,
        expected,
        rtol=0.0,
        atol=_timestamp_tolerance(sample_rate_hz),
    ):
        raise DXDReaderError(
            "timestamp discontinuity between consecutive Dewesoft chunks"
        )


def _validate_mapping(channel_names: Mapping[str, str]) -> list[str]:
    keys = set(channel_names)
    missing = [key for key in _REQUIRED_KEYS if key not in keys]
    extra = sorted(keys.difference(_REQUIRED_KEYS))
    if missing or extra:
        raise DXDReaderError(
            f"channel mapping must contain exactly {_REQUIRED_KEYS}; "
            f"missing={missing}, extra={extra}"
        )
    requested = [str(channel_names[key]) for key in _REQUIRED_KEYS]
    if len(set(requested)) != len(requested):
        raise DXDReaderError(
            "the five APE roles must map to five distinct Dewesoft channels"
        )
    return requested


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
    """Synthetic/small-input validator retained for unit tests and compatibility."""
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

    reference = times["x_1f"]
    for role in _REQUIRED_KEYS[1:]:
        if not np.allclose(times[role], reference, rtol=0.0, atol=1e-9):
            raise DXDReaderError(
                f"selected channels are not synchronized: x_1f and {role} timestamps differ"
            )

    if declared_sample_rate_hz is not None:
        if not np.isfinite(declared_sample_rate_hz) or declared_sample_rate_hz <= 0:
            raise DXDReaderError(
                f"invalid Dewesoft declared sample rate: {declared_sample_rate_hz}"
            )
        _validate_declared_rate(reference, declared_sample_rate_hz)
        sample_rate_hz = float(declared_sample_rate_hz)
        regular_sampling = True
        relative_timestamps_s = None
    else:
        if reference.size < 2:
            sample_rate_hz = None
            regular_sampling = False
            relative_timestamps_s = reference - reference[0]
        else:
            dt = np.diff(reference)
            median_step = float(np.median(dt))
            tol = max(1e-12, abs(median_step) * 1e-6)
            regular_sampling = bool(
                np.all(np.isfinite(dt))
                and np.all(dt > 0)
                and np.allclose(dt, median_step, rtol=0.0, atol=tol)
            )
            if regular_sampling:
                sample_rate_hz = 1.0 / median_step
                relative_timestamps_s = None
            else:
                sample_rate_hz = None
                relative_timestamps_s = reference - reference[0]

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


class _DewesoftSession:
    def __init__(
        self,
        input_file: str | Path,
        *,
        channel_names: Mapping[str, str],
        wrapper_path: str | Path,
        library_path: str | Path,
    ):
        self.source = Path(input_file).expanduser().resolve()
        if not self.source.is_file():
            raise DXDReaderError(f"input .dxd file not found: {self.source}")
        if self.source.suffix.lower() != ".dxd":
            raise DXDReaderError(f"expected a .dxd input file, got: {self.source.name}")

        self.channel_names = dict(channel_names)
        self.requested_names = _validate_mapping(self.channel_names)
        self.header = _load_header(wrapper_path)
        self.library_path = Path(library_path).expanduser().resolve()
        if not self.library_path.is_file():
            raise DXDReaderError(f"Dewesoft native library not found: {self.library_path}")

        self.dll = None
        self.initialized = False
        self.opened = False
        self.file_info = None
        self.channels_by_role: dict[str, object] = {}
        self.units: dict[str, str] = {}
        self.sample_counts: dict[str, int] = {}
        self.channel_count = 0

    def __enter__(self):
        import ctypes

        self.dll = ctypes.cdll.LoadLibrary(str(self.library_path))
        if int(self.dll.DWInit()) != 0:
            raise DXDReaderError("Dewesoft DWInit() failed")
        self.initialized = True
        if int(self.dll.DWAddReader()) != 0:
            raise DXDReaderError("Dewesoft DWAddReader() failed")

        self.file_info = self.header.DWFileInfo(0, 0, 0)
        if int(
            self.dll.DWOpenDataFile(
                c_char_p(str(self.source).encode()), byref(self.file_info)
            )
        ) != 0:
            raise DXDReaderError(f"Dewesoft could not open: {self.source}")
        self.opened = True

        sample_rate = float(self.file_info.sample_rate)
        if not np.isfinite(sample_rate) or sample_rate <= 0:
            raise DXDReaderError(f"invalid Dewesoft declared sample rate: {sample_rate}")

        self.channel_count = int(self.dll.DWGetChannelListCount())
        if self.channel_count < 0:
            raise DXDReaderError("Dewesoft DWGetChannelListCount() failed")
        channel_list = (self.header.DWChannel * self.channel_count)()
        if int(self.dll.DWGetChannelList(byref(channel_list))) != 0:
            raise DXDReaderError("Dewesoft DWGetChannelList() failed")

        by_name: dict[str, object] = {}
        duplicates: set[str] = set()
        for channel in channel_list:
            name = _decode_c_string(channel.name)
            if name in by_name:
                duplicates.add(name)
            by_name[name] = channel

        requested_duplicates = sorted(set(self.requested_names).intersection(duplicates))
        if requested_duplicates:
            raise DXDReaderError(
                f"requested Dewesoft channel names are duplicated in the file: "
                f"{requested_duplicates}"
            )

        absent = [name for name in self.requested_names if name not in by_name]
        if absent:
            raise DXDReaderError(
                f"requested Dewesoft channels are absent: {absent}; "
                f"available channels={sorted(by_name)}"
            )

        for role in _REQUIRED_KEYS:
            recorded_name = str(self.channel_names[role])
            channel = by_name[recorded_name]
            if int(channel.array_size) != 1:
                raise DXDReaderError(
                    f"{recorded_name} has array_size={channel.array_size}; "
                    "only scalar channels are supported"
                )
            count = int(self.dll.DWGetScaledSamplesCount(c_int(int(channel.index))))
            if count <= 0:
                raise DXDReaderError(
                    f"invalid sample count for {recorded_name}: {count}"
                )
            self.channels_by_role[role] = channel
            self.units[role] = _decode_c_string(channel.unit)
            self.sample_counts[role] = count

        if len(set(self.sample_counts.values())) != 1:
            raise DXDReaderError(
                f"selected channel lengths disagree: {self.sample_counts}"
            )
        return self

    def __exit__(self, exc_type, exc, tb):
        close_error = None
        if self.opened and self.dll is not None:
            if int(self.dll.DWCloseDataFile()) != 0:
                close_error = DXDReaderError("Dewesoft DWCloseDataFile() failed")
        if self.initialized and self.dll is not None:
            if int(self.dll.DWDeInit()) != 0 and close_error is None:
                close_error = DXDReaderError("Dewesoft DWDeInit() failed")
        if close_error is not None and exc is None:
            raise close_error

    @property
    def n_samples(self) -> int:
        return next(iter(self.sample_counts.values()))

    @property
    def sample_rate_hz(self) -> float:
        return float(self.file_info.sample_rate)

    @property
    def start_time(self) -> datetime | None:
        return _ole_datetime(float(self.file_info.start_store_time))

    def file_metadata(self) -> dict:
        return {
            "dewesoft_file_sample_rate_hz": self.sample_rate_hz,
            "dewesoft_duration_s": float(self.file_info.duration),
            "dewesoft_start_store_time_ole_days": float(self.file_info.start_store_time),
            "dewesoft_channel_count": self.channel_count,
            "reader": "DWDataReader",
        }

    def read_role_chunk(
        self, role: str, raw_start: int, count: int
    ) -> tuple[np.ndarray, np.ndarray]:
        channel = self.channels_by_role[role]
        data_buffer = (c_double * count)()
        time_buffer = (c_double * count)()
        status = int(
            self.dll.DWGetScaledSamples(
                c_int(int(channel.index)),
                c_int64(raw_start),
                count,
                data_buffer,
                time_buffer,
            )
        )
        if status != 0:
            raise DXDReaderError(
                f"DWGetScaledSamples() failed for {self.channel_names[role]} "
                f"at raw_start={raw_start}, count={count}"
            )
        data = np.ctypeslib.as_array(data_buffer).copy()
        timestamps = np.ctypeslib.as_array(time_buffer).copy()
        return data, timestamps


def _assemble_chunk(
    session: _DewesoftSession,
    raw_start: int,
    count: int,
    previous_last_timestamp_s: float | None,
) -> tuple[APEDXDChunk, float]:
    arrays: dict[str, np.ndarray] = {}

    reference_data, reference_times = session.read_role_chunk("x_1f", raw_start, count)
    arrays["x_1f"] = reference_data

    if reference_data.shape != (count,) or reference_times.shape != (count,):
        raise DXDReaderError(f"x_1f chunk shape mismatch at raw_start={raw_start}")
    if not np.all(np.isfinite(reference_data)) or not np.all(np.isfinite(reference_times)):
        raise DXDReaderError(f"x_1f contains non-finite values at raw_start={raw_start}")

    _validate_declared_rate(reference_times, session.sample_rate_hz)
    _validate_chunk_boundary(
        previous_last_timestamp_s,
        float(reference_times[0]),
        session.sample_rate_hz,
    )

    for role in _REQUIRED_KEYS[1:]:
        values, times = session.read_role_chunk(role, raw_start, count)
        if values.shape != (count,) or times.shape != (count,):
            raise DXDReaderError(f"{role} chunk shape mismatch at raw_start={raw_start}")
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(times)):
            raise DXDReaderError(f"{role} contains non-finite values at raw_start={raw_start}")
        if not np.allclose(times, reference_times, rtol=0.0, atol=1e-9):
            raise DXDReaderError(
                f"selected channels are not synchronized at raw_start={raw_start}: "
                f"x_1f and {role} timestamps differ"
            )
        arrays[role] = values

    metadata = session.file_metadata()
    chunk = APEDXDChunk(
        source_file=session.source,
        start_index=raw_start,
        stop_index=raw_start + count,
        sample_rate_hz=session.sample_rate_hz,
        x_1f=arrays["x_1f"],
        y_1f=arrays["y_1f"],
        x_2f=arrays["x_2f"],
        y_2f=arrays["y_2f"],
        cavity_transmission=arrays["cavity_transmission"],
        channel_names={key: str(session.channel_names[key]) for key in _REQUIRED_KEYS},
        channel_units=dict(session.units),
        first_timestamp_s=float(reference_times[0]),
        last_timestamp_s=float(reference_times[-1]),
        start_time=session.start_time,
        metadata=metadata,
    )
    return chunk, float(reference_times[-1])


def iter_ape_dxd_chunks(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_size_samples: int = _DEFAULT_CHUNK_SIZE_SAMPLES,
) -> Iterator[APEDXDChunk]:
    """Yield synchronized raw chunks without retaining the full file in memory."""
    if chunk_size_samples <= 1:
        raise DXDReaderError("chunk_size_samples must be greater than 1")

    with _DewesoftSession(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    ) as session:
        previous_last = None
        for raw_start in range(0, session.n_samples, chunk_size_samples):
            count = min(chunk_size_samples, session.n_samples - raw_start)
            chunk, previous_last = _assemble_chunk(
                session, raw_start, count, previous_last
            )
            yield chunk


def inspect_dxd(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_size_samples: int = _DEFAULT_CHUNK_SIZE_SAMPLES,
) -> APEDXDInspection:
    """Validate the entire .dxd file in bounded memory and return metadata only."""
    if chunk_size_samples <= 1:
        raise DXDReaderError("chunk_size_samples must be greater than 1")

    with _DewesoftSession(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    ) as session:
        chunk_count = 0
        previous_last = None
        first_sdk_timestamp_s: float | None = None
        last_sdk_timestamp_s: float | None = None

        for raw_start in range(0, session.n_samples, chunk_size_samples):
            count = min(chunk_size_samples, session.n_samples - raw_start)
            chunk, previous_last = _assemble_chunk(
                session, raw_start, count, previous_last
            )
            if first_sdk_timestamp_s is None:
                first_sdk_timestamp_s = chunk.first_timestamp_s
            last_sdk_timestamp_s = chunk.last_timestamp_s
            chunk_count += 1

        metadata = session.file_metadata()
        metadata.update(
            {
                "first_selected_timestamp_s": first_sdk_timestamp_s,
                "last_selected_timestamp_s": last_sdk_timestamp_s,
                "validation_mode": "full_file_chunked",
                "full_record_retained_in_memory": False,
            }
        )

        return APEDXDInspection(
            source_file=session.source,
            sample_rate_hz=session.sample_rate_hz,
            n_samples=session.n_samples,
            duration_s=session.n_samples / session.sample_rate_hz,
            start_time=session.start_time,
            channel_names={key: str(session.channel_names[key]) for key in _REQUIRED_KEYS},
            channel_units=dict(session.units),
            chunk_size_samples=chunk_size_samples,
            chunk_count=chunk_count,
            metadata=metadata,
        )


def validate_ape_dxd_stream(*args, **kwargs) -> APEDXDInspection:
    """Compatibility alias for the bounded-memory full-file inspection."""
    return inspect_dxd(*args, **kwargs)


def inspect_ape_dxd(*args, **kwargs) -> APEDXDInspection:
    """Alias with explicit APE naming for callers that prefer it."""
    return inspect_dxd(*args, **kwargs)


def read_ape_dxd(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_size_samples: int = _DEFAULT_CHUNK_SIZE_SAMPLES,
) -> APEDXDData:
    """Materialize all chunks in memory.

    This compatibility helper is intentionally not used by the production
    wrapper for large APE records. Use inspect_dxd or iter_ape_dxd_chunks.
    """
    chunks = list(
        iter_ape_dxd_chunks(
            input_file,
            channel_names=channel_names,
            wrapper_path=wrapper_path,
            library_path=library_path,
            chunk_size_samples=chunk_size_samples,
        )
    )
    if not chunks:
        raise DXDReaderError("no data chunks returned")

    def join(name: str) -> np.ndarray:
        return np.concatenate([getattr(chunk, name) for chunk in chunks])

    first = chunks[0]
    return APEDXDData(
        source_file=first.source_file,
        sample_rate_hz=first.sample_rate_hz,
        n_samples=sum(chunk.n_samples for chunk in chunks),
        start_time=first.start_time,
        x_1f=join("x_1f"),
        y_1f=join("y_1f"),
        x_2f=join("x_2f"),
        y_2f=join("y_2f"),
        cavity_transmission=join("cavity_transmission"),
        channel_names=dict(first.channel_names),
        channel_units=dict(first.channel_units),
        regular_sampling=True,
        relative_timestamps_s=None,
        metadata={**first.metadata, "materialized_from_chunks": True},
    )
