#!/usr/bin/env python3
"""Efficient cavity-first search for the longest final-good APE segment.

Workflow
--------
1. Run the approved cavity cleaner over the complete researcher-declared
   continuous file series.
2. Rank cavity-good segments by their maximum possible duration.
3. Starting from the longest candidate, calculate sample-level eta0 only inside
   that candidate. The protected eta0 reference intentionally resets at each
   candidate start and bootstraps from the first valid eta0 block mean.
4. Keep the longest eta0-good subsegment found so far.
5. Stop as soon as that surviving duration is at least as long as every
   unprocessed cavity-good candidate. At that point no unprocessed candidate can
   win, because eta0 cleaning can only shorten a cavity-good segment.

This tool is read-only and writes no data files.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ape_c_lpsd.cleaning import (
    build_cavity_reference_plan,
    build_eta0_segment_plan,
)
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _natural_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def _true_runs(
    mask: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
) -> list[tuple[int, int]]:
    """Return contiguous True runs as half-open raw-sample intervals."""
    mask = np.asarray(mask, dtype=bool)
    starts = np.asarray(starts, dtype=np.int64)
    stops = np.asarray(stops, dtype=np.int64)
    if mask.ndim != 1 or starts.shape != mask.shape or stops.shape != mask.shape:
        raise ValueError("run arrays must be aligned one-dimensional arrays")

    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []

    split = np.flatnonzero(np.diff(idx) > 1)
    first_pos = np.concatenate(([0], split + 1))
    last_pos = np.concatenate((split, [idx.size - 1]))
    first_blocks = idx[first_pos]
    last_blocks = idx[last_pos]
    return [
        (int(starts[first]), int(stops[last]))
        for first, last in zip(first_blocks, last_blocks)
    ]


def _summarize_runs(
    runs: list[tuple[int, int]],
    sample_rate_hz: float,
) -> dict[str, int | float | None]:
    if not runs:
        return {
            "count": 0,
            "total_duration_s": 0.0,
            "max_duration_s": None,
            "min_duration_s": None,
            "median_duration_s": None,
            "p90_duration_s": None,
            "p99_duration_s": None,
            "longest_start_index": None,
            "longest_stop_index_exclusive": None,
            "longest_start_time_s": None,
            "longest_stop_time_s": None,
        }

    durations = np.asarray(
        [(stop - start) / sample_rate_hz for start, stop in runs],
        dtype=np.float64,
    )
    longest_i = int(np.argmax(durations))
    longest_start, longest_stop = runs[longest_i]
    return {
        "count": len(runs),
        "total_duration_s": float(np.sum(durations)),
        "max_duration_s": float(np.max(durations)),
        "min_duration_s": float(np.min(durations)),
        "median_duration_s": float(np.median(durations)),
        "p90_duration_s": float(np.quantile(durations, 0.90)),
        "p99_duration_s": float(np.quantile(durations, 0.99)),
        "longest_start_index": longest_start,
        "longest_stop_index_exclusive": longest_stop,
        "longest_start_time_s": longest_start / sample_rate_hz,
        "longest_stop_time_s": longest_stop / sample_rate_hz,
    }


def _slice_chunk_to_global_interval(
    chunk,
    *,
    global_offset: int,
    start_index: int,
    stop_index: int,
):
    """Return the overlapping portion of one file-local chunk, or None."""
    chunk_global_start = global_offset + chunk.start_index
    chunk_global_stop = global_offset + chunk.stop_index
    overlap_start = max(start_index, chunk_global_start)
    overlap_stop = min(stop_index, chunk_global_stop)
    if overlap_start >= overlap_stop:
        return None

    lo = overlap_start - chunk_global_start
    hi = overlap_stop - chunk_global_start
    fs = chunk.sample_rate_hz
    first_ts = chunk.first_timestamp_s + lo / fs
    last_ts = chunk.first_timestamp_s + (hi - 1) / fs

    return replace(
        chunk,
        start_index=overlap_start,
        stop_index=overlap_stop,
        x_1f=chunk.x_1f[lo:hi],
        y_1f=chunk.y_1f[lo:hi],
        x_2f=chunk.x_2f[lo:hi],
        y_2f=chunk.y_2f[lo:hi],
        cavity_transmission=chunk.cavity_transmission[lo:hi],
        first_timestamp_s=first_ts,
        last_timestamp_s=last_ts,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Find the longest final-good APE segment using cavity-first candidate "
            "ranking and eta0 cleaning only where it can still affect the winner."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument("--input-dir", required=True, help="Directory with .dxd files")
    parser.add_argument("--pattern", default="*.dxd", help="Filename glob")
    parser.add_argument(
        "--chunk-size-samples",
        type=int,
        default=1_000_000,
        help="raw samples per Dewesoft chunk (default: 1000000)",
    )
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)
    input_dir = Path(args.input_dir).expanduser().resolve()
    files = sorted(input_dir.glob(args.pattern), key=_natural_key)
    if not files:
        raise SystemExit(f"no .dxd files matched {args.pattern!r} in {input_dir}")

    file_offsets: dict[Path, int] = {}
    file_sample_counts: dict[Path, int] = {}
    expected_sample_rate: float | None = None

    def cavity_chunk_factory():
        nonlocal expected_sample_rate
        global_offset = 0
        for file_number, source in enumerate(files, start=1):
            print(
                f"Cavity pass: file {file_number}/{len(files)} ({source.name})",
                file=sys.stderr,
                flush=True,
            )
            file_offsets[source] = global_offset
            local_expected_start = 0
            file_count = 0
            for chunk in iter_ape_dxd_chunks(
                source,
                channel_names=config.channels.as_mapping(),
                wrapper_path=dewesoft.wrapper_path,
                library_path=dewesoft.library_path,
                chunk_size_samples=args.chunk_size_samples,
            ):
                if expected_sample_rate is None:
                    expected_sample_rate = chunk.sample_rate_hz
                elif chunk.sample_rate_hz != expected_sample_rate:
                    raise RuntimeError(
                        f"sample-rate change at {source}: "
                        f"{chunk.sample_rate_hz} vs {expected_sample_rate} Hz"
                    )
                if chunk.start_index != local_expected_start:
                    raise RuntimeError(
                        f"non-contiguous chunks inside {source}: "
                        f"expected {local_expected_start}, got {chunk.start_index}"
                    )
                local_start = chunk.start_index
                local_stop = chunk.stop_index
                local_expected_start = local_stop
                file_count = local_stop
                yield replace(
                    chunk,
                    start_index=global_offset + local_start,
                    stop_index=global_offset + local_stop,
                )
            file_sample_counts[source] = file_count
            global_offset += file_count

    print(
        f"Found {len(files)} continuous files. Building cavity-only plan...",
        file=sys.stderr,
        flush=True,
    )
    cavity_plan = build_cavity_reference_plan(cavity_chunk_factory, config.cleaning)
    fs = cavity_plan.sample_rate_hz
    starts = cavity_plan.block_start_indices
    stops = cavity_plan.block_stop_indices

    cavity_runs = _true_runs(cavity_plan.cavity_good_blocks, starts, stops)
    if not cavity_runs:
        raise RuntimeError("cavity cleaner found no good segment")

    ranked = sorted(
        cavity_runs,
        key=lambda interval: (interval[1] - interval[0], -interval[0]),
        reverse=True,
    )

    print(
        f"Cavity pass complete: {len(ranked)} candidate segments. "
        "Testing eta0 only where a candidate can still win...",
        file=sys.stderr,
        flush=True,
    )

    def eta_candidate_factory(candidate_start: int, candidate_stop: int):
        for source in files:
            file_start = file_offsets[source]
            file_stop = file_start + file_sample_counts[source]
            if file_stop <= candidate_start or file_start >= candidate_stop:
                continue

            print(
                f"  eta0 candidate reads {source.name}",
                file=sys.stderr,
                flush=True,
            )
            for chunk in iter_ape_dxd_chunks(
                source,
                channel_names=config.channels.as_mapping(),
                wrapper_path=dewesoft.wrapper_path,
                library_path=dewesoft.library_path,
                chunk_size_samples=args.chunk_size_samples,
            ):
                sliced = _slice_chunk_to_global_interval(
                    chunk,
                    global_offset=file_start,
                    start_index=candidate_start,
                    stop_index=candidate_stop,
                )
                if sliced is not None:
                    yield sliced

    best_start: int | None = None
    best_stop: int | None = None
    best_duration_samples = -1
    evaluated: list[dict[str, int | float]] = []

    for rank, (candidate_start, candidate_stop) in enumerate(ranked, start=1):
        candidate_samples = candidate_stop - candidate_start
        if best_duration_samples >= candidate_samples:
            print(
                f"Stopping before candidate {rank}: current best "
                f"{best_duration_samples / fs:.2f} s >= remaining maximum "
                f"{candidate_samples / fs:.2f} s.",
                file=sys.stderr,
                flush=True,
            )
            break

        print(
            f"Eta0 candidate {rank}/{len(ranked)}: "
            f"{candidate_start / fs:.2f} to {candidate_stop / fs:.2f} s "
            f"({candidate_samples / fs / 3600.0:.3f} h)",
            file=sys.stderr,
            flush=True,
        )

        eta_plan = build_eta0_segment_plan(
            lambda cs=candidate_start, ce=candidate_stop: eta_candidate_factory(cs, ce),
            config.cleaning,
            config.calibration,
        )
        eta_runs = _true_runs(
            eta_plan.eta0_good_blocks,
            eta_plan.block_start_indices,
            eta_plan.block_stop_indices,
        )
        eta_summary = _summarize_runs(eta_runs, fs)
        surviving_s = eta_summary["max_duration_s"]
        surviving_samples = (
            0
            if surviving_s is None
            else int(round(float(surviving_s) * fs))
        )

        evaluated.append(
            {
                "rank": rank,
                "cavity_start_index": candidate_start,
                "cavity_stop_index_exclusive": candidate_stop,
                "cavity_duration_s": candidate_samples / fs,
                "eta_longest_surviving_duration_s": (
                    0.0 if surviving_s is None else float(surviving_s)
                ),
            }
        )

        if surviving_samples > best_duration_samples and eta_runs:
            longest = max(eta_runs, key=lambda interval: interval[1] - interval[0])
            best_start, best_stop = longest
            best_duration_samples = best_stop - best_start

    if best_start is None or best_stop is None:
        raise RuntimeError("eta0 cleaning left no final-good candidate")

    total = cavity_plan.total_samples

    def file_for_index(index: int) -> str | None:
        for source in files:
            start = file_offsets[source]
            stop = start + file_sample_counts[source]
            if start <= index < stop:
                return source.name
        if index == total:
            return files[-1].name
        return None

    next_unprocessed_bound_s = None
    if len(evaluated) < len(ranked):
        next_start, next_stop = ranked[len(evaluated)]
        next_unprocessed_bound_s = (next_stop - next_start) / fs

    summary = {
        "method": "cavity_first_ranked_eta0_candidate_search_v1",
        "continuity_assumption": (
            "researcher-declared continuous acquisition for cavity cleaning and "
            "cavity-good segment construction"
        ),
        "eta0_candidate_reference_policy": (
            "reset protected eta0 history at each cavity-good candidate start "
            "and bootstrap from the first valid eta0 block mean"
        ),
        "file_count": len(files),
        "file_order": [p.name for p in files],
        "n_samples": int(total),
        "duration_h": total / fs / 3600.0,
        "sample_rate_hz": fs,
        "cavity": {
            "candidate_count": len(ranked),
            "good_fraction": float(
                np.sum(
                    (stops - starts)[cavity_plan.cavity_good_blocks],
                    dtype=np.int64,
                )
                / total
            ),
            "good_segments": _summarize_runs(cavity_runs, fs),
        },
        "eta0_search": {
            "evaluated_candidate_count": len(evaluated),
            "total_candidate_count": len(ranked),
            "evaluated_candidates": evaluated,
            "next_unprocessed_candidate_upper_bound_s": next_unprocessed_bound_s,
            "stopped_by_upper_bound": len(evaluated) < len(ranked),
        },
        "longest_good_segment": {
            "duration_s": best_duration_samples / fs,
            "duration_h": best_duration_samples / fs / 3600.0,
            "global_start_index": best_start,
            "global_stop_index_exclusive": best_stop,
            "global_start_time_s": best_start / fs,
            "global_stop_time_s": best_stop / fs,
            "starts_in_file": file_for_index(best_start),
            "ends_in_file": file_for_index(best_stop - 1),
        },
        "cleaning": {
            "block_size_seconds": config.cleaning.block_size_seconds,
            "reference_window_seconds": config.cleaning.reference_window_seconds,
            "reference_quantile": config.cleaning.reference_quantile,
            "cavity_lower_fraction": config.cleaning.cavity_lower_fraction,
            "eta0_lower_fraction": config.cleaning.eta0_lower_fraction,
        },
        "memory_behavior": {
            "raw_full_record_retained": False,
            "writes_output_files": False,
        },
    }
    print(json.dumps(summary, indent=2, sort_keys=False))


if __name__ == "__main__":
    main()
