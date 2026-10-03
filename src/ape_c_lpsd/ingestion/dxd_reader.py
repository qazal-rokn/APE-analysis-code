"""Canonical APE Dewesoft .dxd ingestion.

The production interface is chunked. Large files are read in bounded blocks
and validated without retaining the entire recording in memory.

No calibration, filtering, cleaning, resampling, interpolation, trimming, or
gap filling is performed here.
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

DEFAULT_CHUNK_SIZE_SAMPLES = 1_000_000
_TIMESTAMP_SYNC_ATOL_S = 1e-9


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


def _timestamp_tolerance(sample_rate_hz: float) -> float:
    step = 1.0 / sample_rate_hz
    return max(1e-12, abs(step) * 1e-6)


def _validate_declared_rate(timestamps_s: np.ndarray, sample_rate_hz: float) -> None:
    if timestamps_s.size < 2:
        return
    dt = np.diff(timestamps_s)
    expected = 1.0 / sample_rate_hz
    tolerance = _timestamp_tolerance(sample_rate_hz)
    if (
        not np.all(np.isfinite(dt))
        or np.any(dt <= 0)
        or not np.allclose(dt, expected, rtol=0.0, atol=tolerance)
    ):
        raise DXDReaderError(
            "selected-channel timestamps are inconsistent with the Dewesoft "
            f"declared sample rate {sample_rate_hz} Hz"
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

    def _cleanup(self) -> None:
        close_error = None
        if self.opened and self.dll is not None:
            if int(self.dll.DWCloseDataFile()) != 0:
                close_error = DXDReaderError("Dewesoft DWCloseDataFile() failed")
            self.opened = False
        if self.initialized and self.dll is not None:
            if int(self.dll.DWDeInit()) != 0 and close_error is None:
                close_error = DXDReaderError("Dewesoft DWDeInit() failed")
            self.initialized = False
        if close_error is not None:
            raise close_error

    def __enter__(self):
        import ctypes

        self.dll = ctypes.cdll.LoadLibrary(str(self.library_path))
        try:
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
                raise DXDReaderError(
                    f"invalid Dewesoft declared sample rate: {sample_rate}"
                )

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

            requested_duplicates = sorted(
                set(self.requested_names).intersection(duplicates)
            )
            if requested_duplicates:
                raise DXDReaderError(
                    "requested Dewesoft channel names are duplicated in the file: "
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
                count = int(
                    self.dll.DWGetScaledSamplesCount(c_int(int(channel.index)))
                )
                if count < 0:
                    raise DXDReaderError(
                        f"DWGetScaledSamplesCount() failed for {recorded_name}"
                    )
                if count == 0:
                    raise DXDReaderError(f"{recorded_name} contains no samples")

                self.channels_by_role[role] = channel
                self.units[role] = _decode_c_string(channel.unit)
                self.sample_counts[role] = count

            if len(set(self.sample_counts.values())) != 1:
                raise DXDReaderError(
                    f"selected channel lengths disagree: {self.sample_counts}"
                )
            return self
        except Exception:
            try:
                self._cleanup()
            finally:
                raise

    def __exit__(self, exc_type, exc, tb):
        try:
            self._cleanup()
        except Exception:
            if exc is None:
                raise

    @property
    def n_samples(self) -> int:
        return next(iter(self.sample_counts.values()))

    @property
    def sample_rate_hz(self) -> float:
        return float(self.file_info.sample_rate)

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
    reference_timestamps: np.ndarray | None = None

    for role in _REQUIRED_KEYS:
        data, timestamps = session.read_role_chunk(role, raw_start, count)

        if data.shape != (count,) or timestamps.shape != (count,):
            raise DXDReaderError(
                f"{role} chunk shape mismatch at raw_start={raw_start}"
            )
        if not np.all(np.isfinite(data)) or not np.all(np.isfinite(timestamps)):
            raise DXDReaderError(
                f"{role} contains non-finite values at raw_start={raw_start}"
            )

        arrays[role] = data
        if reference_timestamps is None:
            reference_timestamps = timestamps
        elif not np.allclose(
            timestamps,
            reference_timestamps,
            rtol=0.0,
            atol=_TIMESTAMP_SYNC_ATOL_S,
        ):
            raise DXDReaderError(
                f"selected channels are not synchronized at raw_start={raw_start}: "
                f"x_1f and {role} timestamps differ"
            )

    assert reference_timestamps is not None
    _validate_declared_rate(reference_timestamps, session.sample_rate_hz)

    if previous_last_timestamp_s is not None:
        expected = previous_last_timestamp_s + 1.0 / session.sample_rate_hz
        if not np.isclose(
            reference_timestamps[0],
            expected,
            rtol=0.0,
            atol=_timestamp_tolerance(session.sample_rate_hz),
        ):
            raise DXDReaderError(
                f"timestamp discontinuity between chunks before raw_start={raw_start}"
            )

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
        first_timestamp_s=float(reference_timestamps[0]),
        last_timestamp_s=float(reference_timestamps[-1]),
    )
    return chunk, float(reference_timestamps[-1])


def iter_ape_dxd_chunks(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_size_samples: int = DEFAULT_CHUNK_SIZE_SAMPLES,
) -> Iterator[APEDXDChunk]:
    """Yield synchronized raw chunks without retaining the full file in memory."""
    if int(chunk_size_samples) <= 1:
        raise DXDReaderError("chunk_size_samples must be greater than 1")
    chunk_size_samples = int(chunk_size_samples)

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
    chunk_size_samples: int = DEFAULT_CHUNK_SIZE_SAMPLES,
) -> APEDXDInspection:
    """Validate an entire .dxd file in bounded memory and return metadata only."""
    if int(chunk_size_samples) <= 1:
        raise DXDReaderError("chunk_size_samples must be greater than 1")
    chunk_size_samples = int(chunk_size_samples)

    with _DewesoftSession(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    ) as session:
        chunk_count = 0
        previous_last = None
        first_sdk_timestamp_s: float | None = None

        for raw_start in range(0, session.n_samples, chunk_size_samples):
            count = min(chunk_size_samples, session.n_samples - raw_start)
            chunk, previous_last = _assemble_chunk(
                session, raw_start, count, previous_last
            )
            if first_sdk_timestamp_s is None:
                first_sdk_timestamp_s = chunk.first_timestamp_s
            chunk_count += 1

        metadata = {
            "dewesoft_file_sample_rate_hz": session.sample_rate_hz,
            "dewesoft_duration_s": float(session.file_info.duration),
            "dewesoft_start_store_time_ole_days": float(
                session.file_info.start_store_time
            ),
            "dewesoft_first_sdk_timestamp_s": first_sdk_timestamp_s,
            "dewesoft_channel_count": session.channel_count,
            "reader": "DWDataReader",
            "validation_mode": "full_file_chunked",
        }

        return APEDXDInspection(
            source_file=session.source,
            sample_rate_hz=session.sample_rate_hz,
            n_samples=session.n_samples,
            duration_s=float(session.file_info.duration),
            start_time=_ole_datetime(float(session.file_info.start_store_time)),
            channel_names={
                key: str(session.channel_names[key]) for key in _REQUIRED_KEYS
            },
            channel_units=dict(session.units),
            chunk_size_samples=chunk_size_samples,
            chunk_count=chunk_count,
            metadata=metadata,
        )


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
    """Synthetic/small-input helper retained for controlled tests."""
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
        if not np.allclose(
            times[role], reference, rtol=0.0, atol=_TIMESTAMP_SYNC_ATOL_S
        ):
            raise DXDReaderError(
                f"selected channels are not synchronized: x_1f and {role} timestamps differ"
            )

    if declared_sample_rate_hz is not None:
        if not np.isfinite(declared_sample_rate_hz) or declared_sample_rate_hz <= 0:
            raise DXDReaderError(
                f"invalid Dewesoft declared sample rate: {declared_sample_rate_hz}"
            )
        _validate_declared_rate(reference, declared_sample_rate_hz)
        sample_rate_hz: float | None = float(declared_sample_rate_hz)
        regular_sampling = True
        relative_timestamps_s = None
    else:
        sample_rate_hz = None
        regular_sampling = False
        relative_timestamps_s = reference - reference[0]
        if reference.size >= 2:
            dt = np.diff(reference)
            if np.all(np.isfinite(dt)) and np.all(dt > 0):
                step = float(np.median(dt))
                tolerance = max(1e-12, abs(step) * 1e-6)
                if np.allclose(dt, step, rtol=0.0, atol=tolerance):
                    sample_rate_hz = 1.0 / step
                    regular_sampling = True
                    relative_timestamps_s = None

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
