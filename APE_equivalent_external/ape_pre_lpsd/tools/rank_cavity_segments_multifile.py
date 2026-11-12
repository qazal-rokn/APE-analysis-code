#!/usr/bin/env python3
"""Rank cavity-good segments across a continuous APE file series.

This is the first stage of the efficient longest-good-segment search. It applies
only the approved cavity cleaner, deliberately skipping eta0 calibration.
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

from ape_c_lpsd.cleaning import build_cavity_reference_plan
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _natural_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def _true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    values = np.asarray(mask, dtype=bool)
    padded = np.concatenate(([False], values, [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    return [(int(a), int(b)) for a, b in zip(edges[::2], edges[1::2])]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rank cavity-good segments before any eta0 calculation."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--pattern", default="*.dxd")
    parser.add_argument("--chunk-size-samples", type=int, default=1_000_000)
    parser.add_argument("--top", type=int, default=20)
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
                f"Cavity pass: file {file_number}/{len(files)} ({source.name})",
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

            file_sample_counts[source] = file_count
            global_offset += file_count

    plan = build_cavity_reference_plan(chunk_factory, config.cleaning)
    fs = plan.sample_rate_hz

    runs = []
    for block_start, block_stop in _true_runs(plan.cavity_good_blocks):
        sample_start = int(plan.block_start_indices[block_start])
        sample_stop = int(plan.block_stop_indices[block_stop - 1])
        n_samples = sample_stop - sample_start
        runs.append(
            {
                "global_start_index": sample_start,
                "global_stop_index_exclusive": sample_stop,
                "duration_s": n_samples / fs,
                "duration_h": n_samples / fs / 3600.0,
            }
        )

    runs.sort(key=lambda item: item["duration_s"], reverse=True)

    def file_for_index(index: int) -> str:
        for source in files:
            start = file_offsets[source]
            stop = start + file_sample_counts[source]
            if start <= index < stop:
                return source.name
        if index == plan.total_samples:
            return files[-1].name
        return "UNKNOWN"

    for rank, run in enumerate(runs, start=1):
        run["rank"] = rank
        run["starts_in_file"] = file_for_index(run["global_start_index"])
        run["ends_in_file"] = file_for_index(
            max(run["global_start_index"], run["global_stop_index_exclusive"] - 1)
        )
        run["global_start_time_s"] = run["global_start_index"] / fs
        run["global_stop_time_s"] = run["global_stop_index_exclusive"] / fs

    summary = {
        "file_count": len(files),
        "sample_rate_hz": fs,
        "total_samples": int(plan.total_samples),
        "total_duration_h": plan.total_samples / fs / 3600.0,
        "cavity_good_fraction": float(
            np.count_nonzero(plan.cavity_good_blocks) / plan.cavity_good_blocks.size
        ),
        "cavity_good_segment_count": len(runs),
        "top_cavity_good_segments": runs[: max(1, args.top)],
        "note": (
            "eta0 has not been calculated. These are candidate upper bounds for "
            "the final longest-good segment."
        ),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
