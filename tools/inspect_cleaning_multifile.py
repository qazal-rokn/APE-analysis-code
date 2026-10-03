#!/usr/bin/env python3
"""Read-only staged-cleaning inspection across one continuous APE file series.

The researcher has declared the supplied files to be consecutive parts of one
continuous acquisition. File boundaries therefore do not break good-data runs.
The tool still refuses sample-rate changes or malformed/non-contiguous chunks
within any individual file.
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

from ape_c_lpsd.cleaning import build_cleaning_reference_plan, clean_chunk_with_plan
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths

# Reuse the already-tested streaming interval accumulator.
from inspect_cleaning import _IntervalAccumulator


def _natural_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


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
        raise SystemExit(
            f"no .dxd files matched {args.pattern!r} in {input_dir}"
        )

    file_offsets: dict[Path, int] = {}
    file_sample_counts: dict[Path, int] = {}
    expected_sample_rate: float | None = None
    progress_phase = "Reading"

    def chunk_factory():
        nonlocal expected_sample_rate
        global_offset = 0
        for file_number, source in enumerate(files, start=1):
            print(
                f"{progress_phase}: file {file_number}/{len(files)} "
                f"({source.name})",
                file=sys.stderr,
                flush=True,
            )
            file_offsets.setdefault(source, global_offset)
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
            previous = file_sample_counts.setdefault(source, file_count)
            if previous != file_count:
                raise RuntimeError(
                    f"sample count changed between passes for {source}: "
                    f"{previous} vs {file_count}"
                )
            global_offset += file_count

    print(
        f"Found {len(files)} continuous files. Building reduced-block cleaning references...",
        file=sys.stderr,
        flush=True,
    )
    progress_phase = "Reference pass"
    plan = build_cleaning_reference_plan(
        chunk_factory,
        config.cleaning,
        config.calibration,
    )
    print(
        "Reference pass complete. Inspecting final masks and intervals...",
        file=sys.stderr,
        flush=True,
    )

    total = 0
    cavity_good_count = 0
    eta_good_count = 0
    final_good_count = 0
    cavity_bad_runs = _IntervalAccumulator(target_value=False)
    eta_only_bad_runs = _IntervalAccumulator(target_value=True)
    final_bad_runs = _IntervalAccumulator(target_value=False)
    final_good_runs = _IntervalAccumulator(target_value=True)
    per_file: dict[str, dict[str, int | float]] = {}

    progress_phase = "Inspection pass"
    for chunk in chunk_factory():
        result = clean_chunk_with_plan(
            chunk,
            plan,
            config.cleaning,
            config.calibration,
        )
        total += chunk.n_samples
        cavity_good = int(np.count_nonzero(result.cavity_good))
        eta_good = int(np.count_nonzero(result.eta_good))
        final_good = int(np.count_nonzero(result.final_good))
        cavity_good_count += cavity_good
        eta_good_count += eta_good
        final_good_count += final_good

        cavity_bad_runs.feed(result.cavity_good, chunk.start_index)
        eta_only_bad_runs.feed(
            result.cavity_good & ~result.eta_good,
            chunk.start_index,
        )
        final_bad_runs.feed(result.final_good, chunk.start_index)
        final_good_runs.feed(result.final_good, chunk.start_index)

        name = chunk.source_file.name
        stats = per_file.setdefault(
            name,
            {
                "n_samples": 0,
                "cavity_good_samples": 0,
                "eta0_good_samples": 0,
                "final_good_samples": 0,
            },
        )
        stats["n_samples"] += chunk.n_samples
        stats["cavity_good_samples"] += cavity_good
        stats["eta0_good_samples"] += eta_good
        stats["final_good_samples"] += final_good

    if total != plan.total_samples:
        raise RuntimeError(
            f"inspection sample count changed: plan={plan.total_samples}, observed={total}"
        )

    cavity_bad_runs.finish()
    eta_only_bad_runs.finish()
    final_bad_runs.finish()
    final_good_runs.finish()

    for stats in per_file.values():
        n = int(stats["n_samples"])
        stats["cavity_good_fraction"] = stats["cavity_good_samples"] / n
        stats["eta0_good_fraction"] = stats["eta0_good_samples"] / n
        stats["final_good_fraction"] = stats["final_good_samples"] / n

    fs = plan.sample_rate_hz
    longest = final_good_runs.summary(fs)
    longest_start = longest["longest_start_index"]
    longest_stop = longest["longest_stop_index_exclusive"]

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

    print(
        "Inspection complete. Writing JSON summary...",
        file=sys.stderr,
        flush=True,
    )

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
            "raw_sample_decisions": False,
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
            "cavity_bad": cavity_bad_runs.summary(fs),
            "eta0_only_bad": eta_only_bad_runs.summary(fs),
            "final_bad": final_bad_runs.summary(fs),
            "final_good_segments": longest,
        },
        "longest_good_segment": {
            "duration_s": longest["max_duration_s"],
            "duration_h": (
                None
                if longest["max_duration_s"] is None
                else longest["max_duration_s"] / 3600.0
            ),
            "global_start_index": longest_start,
            "global_stop_index_exclusive": longest_stop,
            "global_start_time_s": longest["longest_start_time_s"],
            "global_stop_time_s": longest["longest_stop_time_s"],
            "starts_in_file": file_for_index(longest_start),
            "ends_in_file": (
                None
                if longest_stop is None
                else file_for_index(max(0, longest_stop - 1))
            ),
        },
        "per_file": per_file,
        "memory_behavior": {
            "raw_full_record_retained": False,
            "full_record_references_retained_at_block_resolution": True,
            "writes_output_files": False,
        },
    }
    print(json.dumps(summary, indent=2, sort_keys=False))


if __name__ == "__main__":
    main()
