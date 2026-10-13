#!/usr/bin/env python3
"""Read-only cleaning inspection across one continuous APE file series.

The researcher has declared the supplied files to be consecutive parts of one
continuous acquisition. File boundaries therefore do not break good-data runs.

The expensive raw data are read exactly once. The cleaning plan already stores
all final quality masks at reduced-block resolution, so diagnostics and the
longest-good-segment search are computed directly from those blocks instead of
re-reading every DXD file merely to expand masks back to raw samples.
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

from ape_c_lpsd.cleaning import build_cleaning_reference_plan
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _natural_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def _summarize_block_runs(
    mask: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
    sample_rate_hz: float,
) -> dict[str, int | float | None]:
    """Summarize True runs directly from contiguous reduced blocks."""
    mask = np.asarray(mask, dtype=bool)
    starts = np.asarray(starts, dtype=np.int64)
    stops = np.asarray(stops, dtype=np.int64)
    if mask.ndim != 1 or starts.shape != mask.shape or stops.shape != mask.shape:
        raise ValueError("block run arrays must be aligned one-dimensional arrays")

    idx = np.flatnonzero(mask)
    if idx.size == 0:
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

    split = np.flatnonzero(np.diff(idx) > 1)
    run_first_positions = np.concatenate(([0], split + 1))
    run_last_positions = np.concatenate((split, [idx.size - 1]))
    first_blocks = idx[run_first_positions]
    last_blocks = idx[run_last_positions]

    run_starts = starts[first_blocks]
    run_stops = stops[last_blocks]
    durations_s = (run_stops - run_starts) / float(sample_rate_hz)
    longest_i = int(np.argmax(durations_s))

    return {
        "count": int(durations_s.size),
        "total_duration_s": float(np.sum(durations_s)),
        "max_duration_s": float(np.max(durations_s)),
        "min_duration_s": float(np.min(durations_s)),
        "median_duration_s": float(np.median(durations_s)),
        "p90_duration_s": float(np.quantile(durations_s, 0.90)),
        "p99_duration_s": float(np.quantile(durations_s, 0.99)),
        "longest_start_index": int(run_starts[longest_i]),
        "longest_stop_index_exclusive": int(run_stops[longest_i]),
        "longest_start_time_s": float(run_starts[longest_i] / sample_rate_hz),
        "longest_stop_time_s": float(run_stops[longest_i] / sample_rate_hz),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect approved cleaning across a researcher-declared continuous "
            "series of Dewesoft .dxd files. Read-only; writes no data."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing the continuous .dxd file series",
    )
    parser.add_argument(
        "--pattern",
        default="*.dxd",
        help="Filename glob inside input-dir (default: *.dxd)",
    )
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

    def chunk_factory():
        nonlocal expected_sample_rate
        global_offset = 0
        for file_number, source in enumerate(files, start=1):
            print(
                f"Cleaning pass: file {file_number}/{len(files)} ({source.name})",
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
                file_count = local_stop
                local_expected_start = local_stop
                yield replace(
                    chunk,
                    start_index=global_offset + local_start,
                    stop_index=global_offset + local_stop,
                )
            file_sample_counts[source] = file_count
            global_offset += file_count

    print(
        f"Found {len(files)} continuous files. Running one raw-data cleaning pass...",
        file=sys.stderr,
        flush=True,
    )
    plan = build_cleaning_reference_plan(
        chunk_factory,
        config.cleaning,
        config.calibration,
    )
    print(
        "Raw-data pass complete. Computing diagnostics from reduced blocks...",
        file=sys.stderr,
        flush=True,
    )

    starts = plan.block_start_indices
    stops = plan.block_stop_indices
    lengths = stops - starts
    total = int(plan.total_samples)

    cavity_good_blocks = plan.cavity_good_blocks
    eta_good_blocks = plan.eta0_good_blocks
    final_good_blocks = plan.final_good_blocks
    eta_only_bad_blocks = cavity_good_blocks & ~eta_good_blocks

    cavity_good_count = int(np.sum(lengths[cavity_good_blocks], dtype=np.int64))
    eta_good_count = int(np.sum(lengths[eta_good_blocks], dtype=np.int64))
    final_good_count = int(np.sum(lengths[final_good_blocks], dtype=np.int64))

    fs = plan.sample_rate_hz
    cavity_bad_summary = _summarize_block_runs(
        ~cavity_good_blocks, starts, stops, fs
    )
    eta_only_bad_summary = _summarize_block_runs(
        eta_only_bad_blocks, starts, stops, fs
    )
    final_bad_summary = _summarize_block_runs(
        ~final_good_blocks, starts, stops, fs
    )
    final_good_summary = _summarize_block_runs(
        final_good_blocks, starts, stops, fs
    )

    per_file: dict[str, dict[str, int | float]] = {}
    for source in files:
        file_start = file_offsets[source]
        file_stop = file_start + file_sample_counts[source]
        b0 = int(np.searchsorted(starts, file_start, side="left"))
        b1 = int(np.searchsorted(starts, file_stop, side="left"))

        if b0 >= b1:
            raise RuntimeError(f"no cleaning blocks found for {source}")
        if int(starts[b0]) != file_start or int(stops[b1 - 1]) != file_stop:
            raise RuntimeError(
                f"file boundary is not aligned to cleaning blocks for {source}"
            )

        local_lengths = lengths[b0:b1]
        local_cavity = cavity_good_blocks[b0:b1]
        local_eta = eta_good_blocks[b0:b1]
        local_final = final_good_blocks[b0:b1]
        n = int(np.sum(local_lengths, dtype=np.int64))

        cavity_n = int(np.sum(local_lengths[local_cavity], dtype=np.int64))
        eta_n = int(np.sum(local_lengths[local_eta], dtype=np.int64))
        final_n = int(np.sum(local_lengths[local_final], dtype=np.int64))

        per_file[source.name] = {
            "n_samples": n,
            "cavity_good_samples": cavity_n,
            "eta0_good_samples": eta_n,
            "final_good_samples": final_n,
            "cavity_good_fraction": cavity_n / n,
            "eta0_good_fraction": eta_n / n,
            "final_good_fraction": final_n / n,
        }

    longest_start = final_good_summary["longest_start_index"]
    longest_stop = final_good_summary["longest_stop_index_exclusive"]

    def file_for_index(index: int | None) -> str | None:
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

    summary = {
        "continuity_assumption": (
            "researcher-declared continuous acquisition; file boundaries do not "
            "break good-data segments"
        ),
        "file_order": [p.name for p in files],
        "file_count": len(files),
        "n_samples": total,
        "duration_s": total / fs,
        "duration_h": total / fs / 3600.0,
        "sample_rate_hz": fs,
        "cleaning": {
            "block_size_seconds": config.cleaning.block_size_seconds,
            "reference_window_seconds": config.cleaning.reference_window_seconds,
            "reference_quantile": config.cleaning.reference_quantile,
            "refinement_iterations": config.cleaning.refinement_iterations,
            "cavity_lower_fraction": config.cleaning.cavity_lower_fraction,
            "eta0_lower_fraction": config.cleaning.eta0_lower_fraction,
            "one_sided_lower_cuts_only": True,
            "reduced_block_decisions": True,
            "raw_sample_eta0_calculation": True,
        },
        "overall": {
            "cavity_good_fraction": cavity_good_count / total,
            "eta0_good_fraction": eta_good_count / total,
            "eta0_good_given_cavity_good": (
                None if cavity_good_count == 0 else eta_good_count / cavity_good_count
            ),
            "final_good_fraction": final_good_count / total,
        },
        "interval_diagnostics": {
            "cavity_bad": cavity_bad_summary,
            "eta0_only_bad": eta_only_bad_summary,
            "final_bad": final_bad_summary,
            "final_good_segments": final_good_summary,
        },
        "longest_good_segment": {
            "duration_s": final_good_summary["max_duration_s"],
            "duration_h": (
                None
                if final_good_summary["max_duration_s"] is None
                else final_good_summary["max_duration_s"] / 3600.0
            ),
            "global_start_index": longest_start,
            "global_stop_index_exclusive": longest_stop,
            "global_start_time_s": final_good_summary["longest_start_time_s"],
            "global_stop_time_s": final_good_summary["longest_stop_time_s"],
            "starts_in_file": file_for_index(longest_start),
            "ends_in_file": (
                None
                if longest_stop is None
                else file_for_index(max(0, int(longest_stop) - 1))
            ),
        },
        "performance": {
            "raw_dxd_passes": 1,
            "diagnostics_from_reduced_blocks": True,
            "python_loop_over_full_0p05s_blocks": False,
        },
        "memory_behavior": {
            "raw_full_record_retained": False,
            "full_record_references_retained_at_block_resolution": True,
            "writes_output_files": False,
        },
    }
    print(json.dumps(summary, indent=2, sort_keys=False))


if __name__ == "__main__":
    main()
