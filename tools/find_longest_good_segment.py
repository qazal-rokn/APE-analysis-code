#!/usr/bin/env python3
"""Efficient cavity-first search for the longest final-good APE segment.

Scientific method:
1. run the unchanged cavity-quality cleaner over the complete continuous record;
2. rank cavity-good segments by duration;
3. for each candidate that can still beat the current best, re-read only that
   candidate and apply the approved protected eta0 cleaner;
4. eta0 protected-reference state is intentionally reset at the beginning of
   each cavity-good candidate segment, per researcher approval;
5. stop when no unprocessed cavity-good segment can exceed the best final-good
   segment already found.

This tool is read-only and writes no data products.
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
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0:
        return []

    padded = np.concatenate(([False], mask, [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    runs = []
    for left, right in edges.reshape(-1, 2):
        runs.append(
            {
                "start_index": int(starts[left]),
                "stop_index_exclusive": int(stops[right - 1]),
                "n_blocks": int(right - left),
            }
        )
    return runs


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Find the longest final-good segment efficiently: full-record cavity "
            "cleaning first, then protected eta0 cleaning only for candidates "
            "that can still become the longest."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument("--input-dir", required=True, help="Directory of continuous DXD files")
    parser.add_argument("--pattern", default="*.dxd", help="Filename glob")
    parser.add_argument(
        "--chunk-size-samples",
        type=int,
        default=1_000_000,
        help="raw samples per DXD chunk (default: 1000000)",
    )
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)
    input_dir = Path(args.input_dir).expanduser().resolve()
    files = sorted(input_dir.glob(args.pattern), key=_natural_key)
    if not files:
        raise SystemExit(f"no DXD files matched {args.pattern!r} in {input_dir}")

    file_offsets: dict[Path, int] = {}
    file_sample_counts: dict[Path, int] = {}
    expected_fs: float | None = None

    def full_chunk_factory():
        nonlocal expected_fs
        global_offset = 0

        for number, source in enumerate(files, start=1):
            print(
                f"Cavity pass: file {number}/{len(files)} ({source.name})",
                file=sys.stderr,
                flush=True,
            )
            file_offsets.setdefault(source, global_offset)
            local_stop = 0

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
                    raise RuntimeError(
                        f"sample-rate change: {chunk.sample_rate_hz} vs {expected_fs}"
                    )

                if chunk.start_index != local_stop:
                    raise RuntimeError(
                        f"non-contiguous chunks in {source.name}: "
                        f"expected {local_stop}, got {chunk.start_index}"
                    )

                local_start = chunk.start_index
                local_stop = chunk.stop_index
                yield replace(
                    chunk,
                    start_index=global_offset + local_start,
                    stop_index=global_offset + local_stop,
                )

            file_sample_counts[source] = local_stop
            global_offset += local_stop

    print(
        f"Found {len(files)} continuous files. Running cavity-only cleaning...",
        file=sys.stderr,
        flush=True,
    )

    cavity_plan = build_cavity_reference_plan(
        full_chunk_factory,
        config.cleaning,
    )

    fs = cavity_plan.sample_rate_hz
    cavity_runs = _true_runs(
        cavity_plan.cavity_good_blocks,
        cavity_plan.block_start_indices,
        cavity_plan.block_stop_indices,
    )

    for run in cavity_runs:
        run["duration_samples"] = (
            run["stop_index_exclusive"] - run["start_index"]
        )
        run["duration_s"] = run["duration_samples"] / fs

    cavity_runs.sort(
        key=lambda item: item["duration_samples"],
        reverse=True,
    )

    print(
        f"Cavity pass complete: {len(cavity_runs)} cavity-good segments.",
        file=sys.stderr,
        flush=True,
    )

    def candidate_chunk_factory(global_start: int, global_stop: int):
        for source in files:
            file_start = file_offsets[source]
            file_stop = file_start + file_sample_counts[source]

            overlap_start = max(global_start, file_start)
            overlap_stop = min(global_stop, file_stop)
            if overlap_start >= overlap_stop:
                continue

            local_start = overlap_start - file_start
            local_stop = overlap_stop - file_start

            for chunk in iter_ape_dxd_chunks(
                source,
                channel_names=config.channels.as_mapping(),
                wrapper_path=dewesoft.wrapper_path,
                library_path=dewesoft.library_path,
                chunk_size_samples=args.chunk_size_samples,
                start_index=local_start,
                stop_index=local_stop,
            ):
                yield replace(
                    chunk,
                    start_index=file_start + chunk.start_index,
                    stop_index=file_start + chunk.stop_index,
                )

    best = None
    processed = []

    for rank, candidate in enumerate(cavity_runs, start=1):
        candidate_samples = candidate["duration_samples"]

        if best is not None and candidate_samples <= best["duration_samples"]:
            print(
                "Stopping: every remaining cavity-good candidate is no longer "
                "than the current best final-good segment.",
                file=sys.stderr,
                flush=True,
            )
            break

        print(
            f"Eta0 candidate {rank}/{len(cavity_runs)}: "
            f"{candidate['duration_s'] / 3600.0:.4f} h "
            f"[{candidate['start_index']}, {candidate['stop_index_exclusive']})",
            file=sys.stderr,
            flush=True,
        )

        eta_plan = build_eta0_segment_plan(
            lambda c=candidate: candidate_chunk_factory(
                c["start_index"],
                c["stop_index_exclusive"],
            ),
            config.cleaning,
            config.calibration,
        )

        eta_runs = _true_runs(
            eta_plan.eta0_good_blocks,
            eta_plan.block_start_indices,
            eta_plan.block_stop_indices,
        )

        if eta_runs:
            longest = max(
                eta_runs,
                key=lambda item: item["stop_index_exclusive"] - item["start_index"],
            )
            longest["duration_samples"] = (
                longest["stop_index_exclusive"] - longest["start_index"]
            )
            longest["duration_s"] = longest["duration_samples"] / fs
        else:
            longest = {
                "start_index": None,
                "stop_index_exclusive": None,
                "duration_samples": 0,
                "duration_s": 0.0,
            }

        processed.append(
            {
                "rank": rank,
                "cavity_start_index": candidate["start_index"],
                "cavity_stop_index_exclusive": candidate["stop_index_exclusive"],
                "cavity_duration_s": candidate["duration_s"],
                "longest_final_start_index": longest["start_index"],
                "longest_final_stop_index_exclusive": longest["stop_index_exclusive"],
                "longest_final_duration_s": longest["duration_s"],
            }
        )

        if best is None or longest["duration_samples"] > best["duration_samples"]:
            best = longest

    if best is None:
        best = {
            "start_index": None,
            "stop_index_exclusive": None,
            "duration_samples": 0,
            "duration_s": 0.0,
        }

    total = cavity_plan.total_samples

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

    best_start = best["start_index"]
    best_stop = best["stop_index_exclusive"]

    summary = {
        "method": "cavity_first_candidate_eta0_search_v1",
        "continuity_assumption": (
            "researcher-declared continuous acquisition across input files"
        ),
        "eta0_candidate_initialization": (
            "protected eta0 reference resets at the beginning of each cavity-good "
            "candidate segment, per researcher approval"
        ),
        "file_count": len(files),
        "file_order": [p.name for p in files],
        "sample_rate_hz": fs,
        "n_samples": total,
        "duration_h": total / fs / 3600.0,
        "cleaning": {
            "block_size_seconds": config.cleaning.block_size_seconds,
            "reference_window_seconds": config.cleaning.reference_window_seconds,
            "reference_quantile": config.cleaning.reference_quantile,
            "cavity_lower_fraction": config.cleaning.cavity_lower_fraction,
            "eta0_lower_fraction": config.cleaning.eta0_lower_fraction,
        },
        "cavity_good_segment_count": len(cavity_runs),
        "largest_cavity_good_segments": [
            {
                "rank": i + 1,
                "start_index": run["start_index"],
                "stop_index_exclusive": run["stop_index_exclusive"],
                "duration_s": run["duration_s"],
                "duration_h": run["duration_s"] / 3600.0,
            }
            for i, run in enumerate(cavity_runs[:10])
        ],
        "eta0_candidates_processed": len(processed),
        "eta0_candidates_skipped_by_bound": len(cavity_runs) - len(processed),
        "processed_candidates": processed,
        "longest_final_good_segment": {
            "global_start_index": best_start,
            "global_stop_index_exclusive": best_stop,
            "global_start_time_s": None if best_start is None else best_start / fs,
            "global_stop_time_s": None if best_stop is None else best_stop / fs,
            "duration_s": best["duration_s"],
            "duration_h": best["duration_s"] / 3600.0,
            "starts_in_file": file_for_index(best_start),
            "ends_in_file": (
                None
                if best_stop is None
                else file_for_index(max(0, best_stop - 1))
            ),
        },
        "memory_behavior": {
            "raw_full_record_retained": False,
            "eta0_computed_only_for_processed_candidates": True,
            "writes_output_files": False,
        },
    }

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
