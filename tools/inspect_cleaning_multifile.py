#!/usr/bin/env python3
"""Find the longest final-good APE segment with cavity-first pruning.

The researcher has declared the input files to be consecutive parts of one
continuous acquisition. The expensive sample-level eta0 calculation is applied
only to cavity-good candidate segments that can still beat the current best.

Approved efficient-search rule: eta0 protected-reference history resets at the
start of each cavity-good candidate segment and bootstraps from that segment's
first valid eta0 block mean.
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


def _true_runs(mask: np.ndarray, starts: np.ndarray, stops: np.ndarray):
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    split = np.flatnonzero(np.diff(idx) != 1) + 1
    groups = np.split(idx, split)
    return [
        (int(starts[g[0]]), int(stops[g[-1]]))
        for g in groups
    ]


def _summarize_runs(runs, fs: float):
    if not runs:
        return {
            "count": 0,
            "total_duration_s": 0.0,
            "max_duration_s": None,
            "longest_start_index": None,
            "longest_stop_index_exclusive": None,
            "longest_start_time_s": None,
            "longest_stop_time_s": None,
        }
    durations = np.asarray([stop - start for start, stop in runs], dtype=np.int64)
    i = int(np.argmax(durations))
    start, stop = runs[i]
    return {
        "count": len(runs),
        "total_duration_s": float(np.sum(durations) / fs),
        "max_duration_s": float(durations[i] / fs),
        "longest_start_index": start,
        "longest_stop_index_exclusive": stop,
        "longest_start_time_s": float(start / fs),
        "longest_stop_time_s": float(stop / fs),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find longest final-good segment with cavity-first eta0 pruning."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--pattern", default="*.dxd")
    parser.add_argument("--chunk-size-samples", type=int, default=1_000_000)
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)
    input_dir = Path(args.input_dir).expanduser().resolve()
    files = sorted(input_dir.glob(args.pattern), key=_natural_key)
    if not files:
        raise SystemExit(f"no .dxd files matched {args.pattern!r} in {input_dir}")

    file_offsets: dict[Path, int] = {}
    file_sample_counts: dict[Path, int] = {}
    expected_fs: float | None = None

    def cavity_chunk_factory():
        nonlocal expected_fs
        global_offset = 0
        for file_number, source in enumerate(files, start=1):
            print(
                f"Cavity pass: file {file_number}/{len(files)} ({source.name})",
                file=sys.stderr,
                flush=True,
            )
            file_offsets[source] = global_offset
            local_expected = 0
            file_count = 0
            for chunk in iter_ape_dxd_chunks(
                source,
                channel_names=config.channels.as_mapping(),
                wrapper_path=dewesoft.wrapper_path,
                library_path=dewesoft.library_path,
                chunk_size_samples=args.chunk_size_samples,
            ):
                if expected_fs is None:
                    expected_fs = chunk.sample_rate_hz
                elif chunk.sample_rate_hz != expected_fs:
                    raise RuntimeError("sample-rate change across input series")
                if chunk.start_index != local_expected:
                    raise RuntimeError(f"non-contiguous chunks inside {source}")
                local_expected = chunk.stop_index
                file_count = chunk.stop_index
                yield replace(
                    chunk,
                    start_index=global_offset + chunk.start_index,
                    stop_index=global_offset + chunk.stop_index,
                )
            file_sample_counts[source] = file_count
            global_offset += file_count

    print(
        f"Found {len(files)} files. Running cavity-only pass...",
        file=sys.stderr,
        flush=True,
    )
    cavity_plan = build_cavity_reference_plan(cavity_chunk_factory, config.cleaning)
    fs = cavity_plan.sample_rate_hz

    cavity_runs = _true_runs(
        cavity_plan.cavity_good_blocks,
        cavity_plan.block_start_indices,
        cavity_plan.block_stop_indices,
    )
    candidates = sorted(
        cavity_runs,
        key=lambda run: run[1] - run[0],
        reverse=True,
    )

    print(
        f"Cavity pass complete: {len(candidates)} cavity-good candidate segments.",
        file=sys.stderr,
        flush=True,
    )

    def candidate_chunk_factory(candidate_start: int, candidate_stop: int):
        for source in files:
            file_start = file_offsets[source]
            file_stop = file_start + file_sample_counts[source]
            if file_stop <= candidate_start or file_start >= candidate_stop:
                continue

            print(
                f"  eta0 reads {source.name}",
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
                global_chunk_start = file_start + chunk.start_index
                global_chunk_stop = file_start + chunk.stop_index

                if global_chunk_stop <= candidate_start:
                    continue
                if global_chunk_start >= candidate_stop:
                    break

                overlap_start = max(candidate_start, global_chunk_start)
                overlap_stop = min(candidate_stop, global_chunk_stop)
                lo = overlap_start - global_chunk_start
                hi = overlap_stop - global_chunk_start
                if hi <= lo:
                    continue

                first_ts = chunk.first_timestamp_s + lo / chunk.sample_rate_hz
                last_ts = chunk.first_timestamp_s + (hi - 1) / chunk.sample_rate_hz

                yield replace(
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

    best_start = None
    best_stop = None
    best_length = 0
    evaluated = []

    for rank, (candidate_start, candidate_stop) in enumerate(candidates, start=1):
        candidate_length = candidate_stop - candidate_start
        if candidate_length <= best_length:
            break

        print(
            f"Eta0 candidate {rank}: {candidate_length/fs/3600.0:.6f} h "
            f"[{candidate_start/fs:.2f}, {candidate_stop/fs:.2f}) s",
            file=sys.stderr,
            flush=True,
        )

        eta_plan = build_eta0_segment_plan(
            lambda s=candidate_start, e=candidate_stop: candidate_chunk_factory(s, e),
            config.cleaning,
            config.calibration,
        )
        eta_runs = _true_runs(
            eta_plan.eta0_good_blocks,
            eta_plan.block_start_indices,
            eta_plan.block_stop_indices,
        )
        eta_summary = _summarize_runs(eta_runs, fs)
        evaluated.append(
            {
                "rank": rank,
                "candidate_start_index": candidate_start,
                "candidate_stop_index_exclusive": candidate_stop,
                "candidate_duration_s": candidate_length / fs,
                "longest_eta_good_duration_s": eta_summary["max_duration_s"],
            }
        )

        if eta_summary["max_duration_s"] is not None:
            eta_best_length = int(
                eta_summary["longest_stop_index_exclusive"]
                - eta_summary["longest_start_index"]
            )
            if eta_best_length > best_length:
                best_length = eta_best_length
                best_start = int(eta_summary["longest_start_index"])
                best_stop = int(eta_summary["longest_stop_index_exclusive"])

    total = int(cavity_plan.total_samples)

    def file_for_index(index: int | None):
        if index is None:
            return None
        for source in files:
            start = file_offsets[source]
            stop = start + file_sample_counts[source]
            if start <= index < stop:
                return source.name
        if index == total:
            return files[-1].name
        return None

    cavity_bad_runs = _true_runs(
        ~cavity_plan.cavity_good_blocks,
        cavity_plan.block_start_indices,
        cavity_plan.block_stop_indices,
    )

    summary = {
        "method": "cavity_first_ranked_eta0_candidate_search_v1",
        "continuity_assumption": (
            "researcher-declared continuous acquisition across file boundaries"
        ),
        "eta0_candidate_start_policy": (
            "reset protected eta0 history at each cavity-good candidate start "
            "and bootstrap from first valid eta0 block mean"
        ),
        "file_count": len(files),
        "file_order": [p.name for p in files],
        "n_samples": total,
        "duration_h": total / fs / 3600.0,
        "sample_rate_hz": fs,
        "cavity": {
            "good_segment_count": len(candidates),
            "bad_intervals": _summarize_runs(cavity_bad_runs, fs),
            "longest_good_segment_duration_h": (
                None if not candidates else (candidates[0][1] - candidates[0][0]) / fs / 3600.0
            ),
        },
        "eta0_search": {
            "candidate_count_total": len(candidates),
            "candidate_count_evaluated": len(evaluated),
            "candidate_count_skipped_by_bound": len(candidates) - len(evaluated),
            "evaluated_candidates": evaluated,
        },
        "longest_good_segment": {
            "duration_s": None if best_start is None else best_length / fs,
            "duration_h": None if best_start is None else best_length / fs / 3600.0,
            "global_start_index": best_start,
            "global_stop_index_exclusive": best_stop,
            "global_start_time_s": None if best_start is None else best_start / fs,
            "global_stop_time_s": None if best_stop is None else best_stop / fs,
            "starts_in_file": file_for_index(best_start),
            "ends_in_file": (
                None if best_stop is None else file_for_index(max(0, best_stop - 1))
            ),
        },
        "performance": {
            "full_record_sample_level_eta0": False,
            "cavity_only_full_record_first": True,
            "eta0_only_for_candidates_that_can_still_win": True,
            "early_stop_bound": "stop when next cavity candidate length <= current best final-good length",
        },
        "writes_output_files": False,
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
