"""Canonical bounded-memory APE Dewesoft .dxd ingestion.

Exactly five raw channels are read. This module performs no calibration,
filtering, cleaning, interpolation, resampling, trimming, or gap filling.
"""
from __future__ import annotations

from ctypes import byref, c_char_p, c_double, c_int, c_int64
from datetime import datetime, timedelta
import importlib.util
from pathlib import Path
from typing import Iterator, Mapping

import numpy as np

from .models import APECavityChunk, APEDXDChunk, APEDXDData, APEDXDInspection, APEDXDMetadata


class DXDReaderError(RuntimeError):
    pass


_REQUIRED_KEYS = (
    "x_1f",
    "y_1f",
    "x_2f",
    "y_2f",
    "cavity_transmission",
)


def _decode_c_string(value) -> str:
    return bytes(value).split(b"\x00", 1)[0].decode("utf-8", errors="replace")


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
    expected = 1.0 / declared_sample_rate_hz
    tolerance = max(1e-12, abs(expected) * 1e-6)
    return bool(np.allclose(dt, expected, rtol=0.0, atol=tolerance))


def _validate_mapping(channel_names: Mapping[str, str]) -> list[str]:
    keys = set(channel_names)
    missing = [k for k in _REQUIRED_KEYS if k not in keys]
    extra = sorted(keys.difference(_REQUIRED_KEYS))
    if missing or extra:
        raise DXDReaderError(
            f"channel mapping must contain exactly {_REQUIRED_KEYS}; missing={missing}, extra={extra}"
        )
    names = [str(channel_names[k]) for k in _REQUIRED_KEYS]
    if len(set(names)) != len(names):
        raise DXDReaderError("the five APE roles must map to five distinct Dewesoft channels")
    return names


def _validate_source(input_file: str | Path) -> Path:
    source = Path(input_file).expanduser().resolve()
    if not source.is_file():
        raise DXDReaderError(f"input .dxd file not found: {source}")
    if source.suffix.lower() != ".dxd":
        raise DXDReaderError(f"expected a .dxd input file, got: {source.name}")
    return source


class _DewesoftSession:
    def __init__(
        self,
        input_file: str | Path,
        *,
        channel_names: Mapping[str, str],
        wrapper_path: str | Path,
        library_path: str | Path,
    ):
        self.source = _validate_source(input_file)
        self.channel_names = dict(channel_names)
        _validate_mapping(self.channel_names)
        self.wrapper_path = Path(wrapper_path).expanduser().resolve()
        self.library_path = Path(library_path).expanduser().resolve()
        self.dll = None
        self.opened = False

    def __enter__(self):
        import ctypes

        header = _load_header(self.wrapper_path)
        if not self.library_path.is_file():
            raise DXDReaderError(f"Dewesoft native library not found: {self.library_path}")
        self.dll = ctypes.cdll.LoadLibrary(str(self.library_path))
        if int(self.dll.DWInit()) != 0:
            raise DXDReaderError("Dewesoft DWInit() failed")
        try:
            if int(self.dll.DWAddReader()) != 0:
                raise DXDReaderError("Dewesoft DWAddReader() failed")
            self.file_info = header.DWFileInfo(0, 0, 0)
            if int(
                self.dll.DWOpenDataFile(
                    c_char_p(str(self.source).encode()), byref(self.file_info)
                )
            ) != 0:
                raise DXDReaderError(f"Dewesoft could not open: {self.source}")
            self.opened = True

            count = int(self.dll.DWGetChannelListCount())
            if count < 0:
                raise DXDReaderError("Dewesoft DWGetChannelListCount() failed")
            channel_list = (header.DWChannel * count)()
            if int(self.dll.DWGetChannelList(byref(channel_list))) != 0:
                raise DXDReaderError("Dewesoft DWGetChannelList() failed")
            self.channel_count = count

            by_name = {}
            duplicates = set()
            for channel in channel_list:
                name = _decode_c_string(channel.name)
                if name in by_name:
                    duplicates.add(name)
                by_name[name] = channel

            requested = _validate_mapping(self.channel_names)
            duplicate_requested = sorted(set(requested).intersection(duplicates))
            if duplicate_requested:
                raise DXDReaderError(
                    f"requested Dewesoft channel names are duplicated: {duplicate_requested}"
                )
            absent = [name for name in requested if name not in by_name]
            if absent:
                raise DXDReaderError(
                    f"requested Dewesoft channels are absent: {absent}; "
                    f"available channels={sorted(by_name)}"
                )

            self.selected = {
                role: by_name[str(self.channel_names[role])] for role in _REQUIRED_KEYS
            }
            counts = {}
            self.units = {}
            for role, channel in self.selected.items():
                if int(channel.array_size) != 1:
                    raise DXDReaderError(
                        f"{self.channel_names[role]} has array_size={channel.array_size}; "
                        "only scalar channels are supported"
                    )
                n = int(self.dll.DWGetScaledSamplesCount(c_int(int(channel.index))))
                if n <= 0:
                    raise DXDReaderError(
                        f"DWGetScaledSamplesCount() failed or returned no samples for "
                        f"{self.channel_names[role]}"
                    )
                counts[role] = n
                self.units[role] = _decode_c_string(channel.unit)

            if len(set(counts.values())) != 1:
                raise DXDReaderError(f"selected channel lengths disagree: {counts}")
            self.n_samples = next(iter(counts.values()))
            self.sample_rate_hz = float(self.file_info.sample_rate)
            if not np.isfinite(self.sample_rate_hz) or self.sample_rate_hz <= 0:
                raise DXDReaderError(
                    f"invalid Dewesoft declared sample rate: {self.sample_rate_hz}"
                )
            return self
        except Exception:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, exc_type, exc, tb):
        if self.dll is None:
            return False
        close_error = None
        if self.opened:
            if int(self.dll.DWCloseDataFile()) != 0:
                close_error = DXDReaderError("Dewesoft DWCloseDataFile() failed")
        if int(self.dll.DWDeInit()) != 0 and close_error is None:
            close_error = DXDReaderError("Dewesoft DWDeInit() failed")
        self.opened = False
        self.dll = None
        if exc_type is None and close_error is not None:
            raise close_error
        return False

    def metadata(self) -> APEDXDMetadata:
        payload = {
            "dewesoft_file_sample_rate_hz": self.sample_rate_hz,
            "dewesoft_duration_s": float(self.file_info.duration),
            "dewesoft_start_store_time_ole_days": float(self.file_info.start_store_time),
            "dewesoft_channel_count": self.channel_count,
            "reader": "DWDataReader",
            "read_mode": "streaming",
        }
        return APEDXDMetadata(
            source_file=self.source,
            sample_rate_hz=self.sample_rate_hz,
            n_samples=self.n_samples,
            duration_s=float(self.file_info.duration),
            start_time=_ole_datetime(float(self.file_info.start_store_time)),
            channel_names={k: str(self.channel_names[k]) for k in _REQUIRED_KEYS},
            channel_units=dict(self.units),
            dewesoft_channel_count=self.channel_count,
            metadata=payload,
        )

    def read_role_chunk(self, role: str, raw_start: int, count: int):
        channel = self.selected[role]
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
                f"at start={raw_start}, count={count}"
            )
        return (
            np.ctypeslib.as_array(data_buffer).copy(),
            np.ctypeslib.as_array(time_buffer).copy(),
        )


def _validate_chunk_arrays(
    *,
    arrays: Mapping[str, np.ndarray],
    timestamps: Mapping[str, np.ndarray],
    sample_rate_hz: float,
) -> tuple[float, float]:
    lengths = {k: len(np.asarray(arrays[k])) for k in _REQUIRED_KEYS}
    if len(set(lengths.values())) != 1 or next(iter(lengths.values())) <= 0:
        raise DXDReaderError(f"chunk channel lengths disagree: {lengths}")

    reference = np.asarray(timestamps[_REQUIRED_KEYS[0]], dtype=np.float64)
    if not np.all(np.isfinite(reference)):
        raise DXDReaderError("chunk contains non-finite timestamps")

    for role in _REQUIRED_KEYS:
        values = np.asarray(arrays[role], dtype=np.float64)
        times = np.asarray(timestamps[role], dtype=np.float64)
        if values.ndim != 1 or times.ndim != 1:
            raise DXDReaderError(f"{role} chunk must be one-dimensional")
        if not np.all(np.isfinite(values)):
            raise DXDReaderError(f"{role} contains non-finite raw samples")
        if times.shape != reference.shape or not np.allclose(
            times, reference, rtol=0.0, atol=1e-9
        ):
            raise DXDReaderError(
                f"selected channels are not synchronized within chunk: "
                f"{_REQUIRED_KEYS[0]} vs {role}"
            )

    if not _timestamps_match_declared_rate(reference, sample_rate_hz):
        raise DXDReaderError(
            "chunk timestamps are inconsistent with Dewesoft declared sample rate "
            f"{sample_rate_hz} Hz"
        )
    return float(reference[0]), float(reference[-1])


def _validate_chunk_boundary(
    *,
    previous_last_timestamp_s: float,
    current_first_timestamp_s: float,
    sample_rate_hz: float,
) -> None:
    expected = previous_last_timestamp_s + 1.0 / sample_rate_hz
    tolerance = max(1e-12, (1.0 / sample_rate_hz) * 1e-6)
    if not np.isclose(current_first_timestamp_s, expected, rtol=0.0, atol=tolerance):
        raise DXDReaderError(
            "timestamp discontinuity between chunks: "
            f"expected {expected}, observed {current_first_timestamp_s}"
        )


def _assemble_chunk(
    session: _DewesoftSession,
    raw_start: int,
    count: int,
    previous_last: float | None,
) -> tuple[APEDXDChunk, float]:
    arrays = {}
    timestamps = {}
    for role in _REQUIRED_KEYS:
        arrays[role], timestamps[role] = session.read_role_chunk(role, raw_start, count)

    first_ts, last_ts = _validate_chunk_arrays(
        arrays=arrays,
        timestamps=timestamps,
        sample_rate_hz=session.sample_rate_hz,
    )
    if previous_last is not None:
        _validate_chunk_boundary(
            previous_last_timestamp_s=previous_last,
            current_first_timestamp_s=first_ts,
            sample_rate_hz=session.sample_rate_hz,
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
        first_timestamp_s=first_ts,
        last_timestamp_s=last_ts,
    )
    return chunk, last_ts


def iter_ape_dxd_cavity_chunks(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_size_samples: int = 1_000_000,
) -> Iterator[APECavityChunk]:
    """Yield validated cavity-only chunks without reading the other raw arrays."""
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
            values, timestamps = session.read_role_chunk(
                "cavity_transmission", raw_start, count
            )
            if not np.all(np.isfinite(values)):
                raise DXDReaderError("cavity_transmission contains non-finite raw samples")
            if not _timestamps_match_declared_rate(timestamps, session.sample_rate_hz):
                raise DXDReaderError(
                    "cavity chunk timestamps are inconsistent with Dewesoft declared "
                    f"sample rate {session.sample_rate_hz} Hz"
                )
            first_ts = float(timestamps[0])
            last_ts = float(timestamps[-1])
            if previous_last is not None:
                _validate_chunk_boundary(
                    previous_last_timestamp_s=previous_last,
                    current_first_timestamp_s=first_ts,
                    sample_rate_hz=session.sample_rate_hz,
                )
            previous_last = last_ts
            yield APECavityChunk(
                source_file=session.source,
                start_index=raw_start,
                stop_index=raw_start + count,
                sample_rate_hz=session.sample_rate_hz,
                cavity_transmission=values,
                first_timestamp_s=first_ts,
                last_timestamp_s=last_ts,
            )


def iter_ape_dxd_range_chunks(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    start_index: int,
    stop_index: int,
    chunk_size_samples: int = 1_000_000,
) -> Iterator[APEDXDChunk]:
    """Yield synchronized five-channel chunks for one half-open raw index range."""
    if chunk_size_samples <= 1:
        raise DXDReaderError("chunk_size_samples must be greater than 1")
    if start_index < 0 or stop_index <= start_index:
        raise DXDReaderError("invalid requested DXD sample range")

    with _DewesoftSession(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    ) as session:
        if stop_index > session.n_samples:
            raise DXDReaderError(
                f"requested stop_index={stop_index} exceeds file length {session.n_samples}"
            )
        previous_last = None
        for raw_start in range(start_index, stop_index, chunk_size_samples):
            count = min(chunk_size_samples, stop_index - raw_start)
            chunk, previous_last = _assemble_chunk(
                session, raw_start, count, previous_last
            )
            yield chunk


def iter_ape_dxd_chunks(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_size_samples: int = 1_000_000,
) -> Iterator[APEDXDChunk]:
    """Yield synchronized validated chunks without retaining the full file."""
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
    chunk_size_samples: int = 1_000_000,
) -> APEDXDInspection:
    """Validate the entire file chunk-by-chunk and return metadata only."""
    if chunk_size_samples <= 1:
        raise DXDReaderError("chunk_size_samples must be greater than 1")
    with _DewesoftSession(
        input_file,
        channel_names=channel_names,
        wrapper_path=wrapper_path,
        library_path=library_path,
    ) as session:
        metadata = session.metadata()
        previous_last = None
        chunk_count = 0
        first_timestamp_s = None

        for raw_start in range(0, session.n_samples, chunk_size_samples):
            count = min(chunk_size_samples, session.n_samples - raw_start)
            chunk, previous_last = _assemble_chunk(
                session, raw_start, count, previous_last
            )
            if first_timestamp_s is None:
                first_timestamp_s = chunk.first_timestamp_s
            chunk_count += 1

        payload = dict(metadata.metadata)
        payload.update(
            {
                "validation_mode": "full_file_chunked",
                "chunk_size_samples": chunk_size_samples,
                "chunk_count": chunk_count,
                "dewesoft_first_sdk_timestamp_s": first_timestamp_s,
            }
        )
        return APEDXDInspection(
            source_file=metadata.source_file,
            sample_rate_hz=metadata.sample_rate_hz,
            n_samples=metadata.n_samples,
            duration_s=metadata.duration_s,
            start_time=metadata.start_time,
            channel_names=metadata.channel_names,
            channel_units=metadata.channel_units,
            dewesoft_channel_count=metadata.dewesoft_channel_count,
            chunk_size_samples=chunk_size_samples,
            chunk_count=chunk_count,
            metadata=payload,
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
    """Compatibility helper for small synthetic/unit-test arrays."""
    missing = [k for k in _REQUIRED_KEYS if k not in arrays or k not in timestamps]
    if missing:
        raise DXDReaderError(f"missing required APE channel data: {missing}")
    values = {k: np.asarray(arrays[k], dtype=np.float64) for k in _REQUIRED_KEYS}
    times = {k: np.asarray(timestamps[k], dtype=np.float64) for k in _REQUIRED_KEYS}
    lengths = {k: values[k].size for k in _REQUIRED_KEYS}
    if len(set(lengths.values())) != 1:
        raise DXDReaderError(f"selected channel lengths disagree: {lengths}")
    n_samples = next(iter(lengths.values()))
    if n_samples <= 0:
        raise DXDReaderError("selected channels contain no samples")

    reference = times[_REQUIRED_KEYS[0]]
    for role in _REQUIRED_KEYS:
        if values[role].ndim != 1 or times[role].ndim != 1:
            raise DXDReaderError(f"{role} must be one-dimensional")
        if not np.all(np.isfinite(values[role])) or not np.all(np.isfinite(times[role])):
            raise DXDReaderError(f"{role} contains non-finite samples or timestamps")
        if times[role].shape != reference.shape or not np.allclose(
            times[role], reference, rtol=0.0, atol=1e-9
        ):
            raise DXDReaderError(
                f"selected channels are not synchronized: {_REQUIRED_KEYS[0]} and {role}"
            )

    if declared_sample_rate_hz is not None:
        if not _timestamps_match_declared_rate(reference, declared_sample_rate_hz):
            raise DXDReaderError(
                "selected-channel timestamps are inconsistent with the Dewesoft "
                f"declared sample rate {declared_sample_rate_hz} Hz"
            )
        fs = float(declared_sample_rate_hz)
        regular = True
        relative = None
    else:
        fs = _regular_rate(reference)
        regular = fs is not None
        relative = None if regular else reference - reference[0]

    return APEDXDData(
        source_file=source_file,
        sample_rate_hz=fs,
        n_samples=n_samples,
        start_time=start_time,
        x_1f=values["x_1f"],
        y_1f=values["y_1f"],
        x_2f=values["x_2f"],
        y_2f=values["y_2f"],
        cavity_transmission=values["cavity_transmission"],
        channel_names={k: str(channel_names[k]) for k in _REQUIRED_KEYS},
        channel_units={k: str(units.get(k, "")) for k in _REQUIRED_KEYS},
        regular_sampling=regular,
        relative_timestamps_s=relative,
        metadata=dict(metadata),
    )


def read_ape_dxd(
    input_file: str | Path,
    *,
    channel_names: Mapping[str, str],
    wrapper_path: str | Path,
    library_path: str | Path,
    chunk_size_samples: int = 1_000_000,
) -> APEDXDData:
    """Compatibility full-memory reader; avoid for production-size files."""
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
    n_samples = sum(chunk.n_samples for chunk in chunks)
    return APEDXDData(
        source_file=first.source_file,
        sample_rate_hz=first.sample_rate_hz,
        n_samples=n_samples,
        start_time=None,
        x_1f=join("x_1f"),
        y_1f=join("y_1f"),
        x_2f=join("x_2f"),
        y_2f=join("y_2f"),
        cavity_transmission=join("cavity_transmission"),
        channel_names={k: str(channel_names[k]) for k in _REQUIRED_KEYS},
        channel_units={k: "" for k in _REQUIRED_KEYS},
        regular_sampling=True,
        relative_timestamps_s=None,
        metadata={"reader": "DWDataReader", "materialized_from_chunks": True},
    )
