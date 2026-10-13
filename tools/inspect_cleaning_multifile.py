#!/usr/bin/env python3
"""Efficient longest-good-segment search across a continuous APE file series.

Scientific workflow:
1. read only cavity transmission across the complete acquisition;
2. build the approved full-record cavity mask;
3. rank cavity-good segments by their maximum possible final duration;
4. read only X2f/Y2f/cavity for candidate segments that can still win;
5. reset the protected eta0 reference at each cavity-good segment start,
   as explicitly approved by the researcher;
6. stop once no unprocessed cavity segment can exceed the best final-good run.

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
from ape_c_lpsd.ingestion import (
    iter_ape_dxd_cavity_chunks,
    iter_ape_dxd_eta_chunks,
    load_dxd_ingestion_config,
)
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _natural_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def _true_runs(
    mask: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
) -> list[tuple[int, int]]:
    """Return half-open sample-index runs for contiguous True reduced blocks."""
    mask = np.asarray(mask, dtype=bool)
    starts = np.asarray(starts, dtype=np.int64)
    stops = np.asarray(stops, dtype=np.int64)
    if mask.ndim != 1 or starts.shape != mask.shape or stops.shape != mask.shape:
        raise ValueError("run arrays must be aligned one-dimensional arrays")

    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []

    split = np.flatnonzero(np.diff(idx) > 1)
    first_positions = np.concatenate(([0], split + 1))
    last_positions = np.concatenate((split, [idx.size - 1]))

    return [
        (int(starts[idx[first]]), int(stops[idx[last]]))
        for first, last in zip(first_positions, last_positions)
    ]


def _run_summary(
    runs: list[tuple[int, int]],
    sample_rate_hz: float,
) -> dict[str, int | float | None]:
    if not runs:
        return {
            "count": 0,
            "total_duration_s": 0.0,
            "max_duration_s": None,
            "longest_start_index": None,
            "longest_stop_index_exclusive": None,
        }

    lengths = np.asarray([stop - start for start, stop in runs], dtype=np.int64)
    longest = int(np.argmax(lengths))
    return {
        "count": len(runs),
        "total_duration_s": float(np.sum(lengths) / sample_rate_hz),
        "max_duration_s": float(lengths[longest] / sample_rate_hz),
        "longest_start_index": int(runs[longest][0]),
        "longest_stop_index_exclusive": int(runs[longest][1]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Find the exact longest final-good APE segment with cavity-first "
            "candidate pruning."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument("--input-dir", required=True, help="continuous .dxd directory")
    parser.add_argument("--pattern", default="*.dxd", help="filename glob")
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

    def cavity_factory():
        nonlocal expected_sample_rate
        global_offset = 0

        for file_number, source in enumerate(files, start=1):
            print(
                f"Cavity pass: file {file_number}/{len(files)} ({source.name})",
                file=sys.stderr,
                flush=True,
            )
            file_offsets[source] = global_offset
            local_count = 0

            for chunk in iter_ape_dxd_cavity_chunks(
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

                local_count = chunk.stop_index
                yield replace(
                    chunk,
                    start_index=global_offset + chunk.start_index,
                    stop_index=global_offset + chunk.stop_index,
                )

            file_sample_counts[source] = local_count
            global_offset += local_count

    print(
        f"Found {len(files)} continuous files. Running cavity-only pass...",
        file=sys.stderr,
        flush=True,
    )
    cavity_plan = build_cavity_reference_plan(cavity_factory, config.cleaning)
    fs = cavity_plan.sample_rate_hz

    cavity_runs = _true_runs(
        cavity_plan.cavity_good_blocks,
        cavity_plan.block_start_indices,
        cavity_plan.block_stop_indices,
    )
    candidates = sorted(
        cavity_runs,
        key=lambda run: (-(run[1] - run[0]), run[0]),
    )

    if not candidates:
        raise RuntimeError("no cavity-good segment exists")

    print(
        f"Cavity pass complete: {len(candidates)} cavity-good segments. "
        "Testing only candidates that can still win...",
        file=sys.stderr,
        flush=True,
    )

    def file_for_index(index: int) -> str:
        total = cavity_plan.total_samples
        if index == total:
            return files[-1].name
        for source in files:
            start = file_offsets[source]
            stop = start + file_sample_counts[source]
            if start <= index < stop:
                return source.name
        raise RuntimeError(f"global index {index} does not map to a source file")

    def eta_candidate_factory(global_start: int, global_stop: int):
        for source in files:
            file_start = file_offsets[source]
            file_stop = file_start + file_sample_counts[source]
            overlap_start = max(global_start, file_start)
            overlap_stop = min(global_stop, file_stop)
            if overlap_stop <= overlap_start:
                continue

            local_start = overlap_start - file_start
            local_stop = overlap_stop - file_start

            for chunk in iter_ape_dxd_eta_chunks(
                source,
                channel_names=config.channels.as_mapping(),
                wrapper_path=dewesoft.wrapper_path,
                library_path=dewesoft.library_path,
                start_index=local_start,
                stop_index=local_stop,
                chunk_size_samples=args.chunk_size_samples,
            ):
                yield replace(
                    chunk,
                    start_index=file_start + chunk.start_index,
                    stop_index=file_start + chunk.stop_index,
                )

    best_start: int | None = None
    best_stop: int | None = None
    best_samples = 0
    processed: list[dict[str, object]] = []

    for rank, (candidate_start, candidate_stop) in enumerate(candidates, start=1):
        candidate_samples = candidate_stop - candidate_start

        # Since candidates are sorted by cavity-good duration, neither this nor
        # any later candidate can beat the current best once this bound is met.
        if candidate_samples <= best_samples:
            break

        print(
            f"Eta0 candidate {rank}: "
            f"{candidate_samples / fs / 3600.0:.4f} h "
            f"[{candidate_start}, {candidate_stop})",
            file=sys.stderr,
            flush=True,
        )

        eta_plan = build_eta0_segment_plan(
            lambda s=candidate_start, e=candidate_stop: eta_candidate_factory(s, e),
            config.cleaning,
            config.calibration,
        )

        eta_runs = _true_runs(
            eta_plan.eta0_good_blocks,
            eta_plan.block_start_indices,
            eta_plan.block_stop_indices,
        )
        eta_summary = _run_summary(eta_runs, fs)

        candidate_best_samples = 0
        candidate_best_start = None
        candidate_best_stop = None
        for start, stop in eta_runs:
            length = stop - start
            if length > candidate_best_samples:
                candidate_best_samples = length
                candidate_best_start = start
                candidate_best_stop = stop

        processed.append(
            {
                "rank": rank,
                "cavity_start_index": candidate_start,
                "cavity_stop_index_exclusive": candidate_stop,
                "cavity_duration_s": candidate_samples / fs,
                "longest_eta_good_duration_s": candidate_best_samples / fs,
                "longest_eta_good_start_index": candidate_best_start,
                "longest_eta_good_stop_index_exclusive": candidate_best_stop,
                "eta_good_run_count": eta_summary["count"],
            }
        )

        if candidate_best_samples > best_samples:
            best_samples = candidate_best_samples
            best_start = candidate_best_start
            best_stop = candidate_best_stop

    if best_start is None or best_stop is None:
        raise RuntimeError("no final-good segment survived eta0 cleaning")

    next_unprocessed_bound_s = None
    if len(processed) < len(candidates):
        next_run = candidates[len(processed)]
        next_unprocessed_bound_s = (next_run[1] - next_run[0]) / fs

    summary = {
        "method": "cavity_first_pruned_longest_segment_search_v1",
        "scientific_policy": {
            "cavity_cleaning": "unchanged full-record centered q90 method",
            "eta0_candidate_start": (
                "reset protected eta0 reference at each cavity-good segment start "
                "and bootstrap from first valid eta0 block mean"
            ),
            "eta0_reference": "past-only accepted-history q90",
            "eta0_lower_fraction": config.cleaning.eta0_lower_fraction,
        },
        "file_order": [p.name for p in files],
        "file_count": len(files),
        "sample_rate_hz": fs,
        "n_samples": cavity_plan.total_samples,
        "duration_h": cavity_plan.total_samples / fs / 3600.0,
        "cavity": {
            "good_segment_count": len(cavity_runs),
            "summary": _run_summary(cavity_runs, fs),
        },
        "candidate_search": {
            "candidate_count_total": len(candidates),
            "candidate_count_eta0_processed": len(processed),
            "processed_candidates": processed,
            "next_unprocessed_max_possible_duration_s": next_unprocessed_bound_s,
            "stopping_rule": (
                "stop when the longest unprocessed cavity segment cannot exceed "
                "the best final-good segment already found"
            ),
        },
        "longest_good_segment": {
            "duration_s": best_samples / fs,
            "duration_h": best_samples / fs / 3600.0,
            "global_start_index": best_start,
            "global_stop_index_exclusive": best_stop,
            "global_start_time_s": best_start / fs,
            "global_stop_time_s": best_stop / fs,
            "starts_in_file": file_for_index(best_start),
            "ends_in_file": file_for_index(best_stop - 1),
        },
        "performance": {
            "full_record_channels_read_first_pass": ["cavity_transmission"],
            "candidate_channels_read": ["x_2f", "y_2f", "cavity_transmission"],
            "x_1f_y_1f_read_during_selection": False,
            "eta0_computed_for_entire_42h": False,
            "writes_output_files": False,
        },
    }

    print(json.dumps(summary, indent=2, sort_keys=False))


if __name__ == "__main__":
    main()
