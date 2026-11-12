#!/usr/bin/env python3
"""Rank cavity-good segments across a continuous APE DXD file series.

This tool intentionally performs only the approved cavity-cleaning stage.
It does not calculate eta0. The purpose is to identify and rank candidate
cavity-good intervals before spending sample-level calibration work on eta0.
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
    _refined_mask,
    _samples_per_block,
)
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _natural_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def _true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 1:
        raise ValueError("mask must be one-dimensional")
    padded = np.concatenate(([False], mask, [False]))
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    return [(int(a), int(b)) for a, b in changes.reshape(-1, 2)]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run only the approved cavity-quality stage across a continuous "
            "series of APE Dewesoft files and rank cavity-good candidates."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument("--input-dir", required=True, help="directory containing DXD files")
    parser.add_argument("--pattern", default="*.dxd", help="filename glob")
    parser.add_argument("--chunk-size-samples", type=int, default=1_000_000)
    parser.add_argument("--top", type=int, default=20, help="number of ranked candidates")
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)
    input_dir = Path(args.input_dir).expanduser().resolve()
    files = sorted(input_dir.glob(args.pattern), key=_natural_key)
    if not files:
        raise SystemExit(f"no .dxd files matched {args.pattern!r} in {input_dir}")

    fs = None
    block_samples = None
    global_offset = 0
    file_offsets: dict[Path, int] = {}
    file_counts: dict[Path, int] = {}

    starts_parts = []
    stops_parts = []
    mids_parts = []
    avg_parts = []
    min_parts = []

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
            if fs is None:
                fs = float(chunk.sample_rate_hz)
                block_samples = _samples_per_block(
                    fs, config.cleaning.block_size_seconds
                )
            elif chunk.sample_rate_hz != fs:
                raise RuntimeError(
                    f"sample-rate change at {source}: "
                    f"{chunk.sample_rate_hz} vs {fs} Hz"
                )

            if chunk.start_index != local_expected:
                raise RuntimeError(
                    f"non-contiguous chunks inside {source}: "
                    f"expected {local_expected}, got {chunk.start_index}"
                )
            if chunk.start_index % block_samples != 0:
                raise RuntimeError("chunk boundary is not cleaning-block aligned")
            if chunk.n_samples % block_samples != 0:
                raise RuntimeError(
                    "this fast cavity-only tool requires block-aligned chunk lengths"
                )

            v = np.asarray(chunk.cavity_transmission, dtype=np.float64)
            if not np.all(np.isfinite(v)):
                raise RuntimeError(f"non-finite cavity samples in {source}")

            n_blocks = chunk.n_samples // block_samples
            vb = v.reshape(n_blocks, block_samples)

            local_starts = (
                chunk.start_index
                + np.arange(n_blocks, dtype=np.int64) * block_samples
            )
            starts = global_offset + local_starts
            stops = starts + block_samples
            mids = (starts + stops) / (2.0 * fs)

            starts_parts.append(starts)
            stops_parts.append(stops)
            mids_parts.append(mids)
            avg_parts.append(np.mean(vb, axis=1))
            min_parts.append(np.min(vb, axis=1))

            local_expected = chunk.stop_index
            file_count = chunk.stop_index

        file_counts[source] = file_count
        global_offset += file_count

    if fs is None or block_samples is None:
        raise RuntimeError("no samples were read")

    starts = np.concatenate(starts_parts)
    stops = np.concatenate(stops_parts)
    mids = np.concatenate(mids_parts)
    cavity_avg = np.concatenate(avg_parts)
    cavity_min = np.concatenate(min_parts)

    cavity_ref, cavity_good = _refined_mask(
        cavity_avg,
        cavity_min,
        mids,
        np.isfinite(cavity_avg),
        window_seconds=config.cleaning.reference_window_seconds,
        quantile=config.cleaning.reference_quantile,
        lower_fraction=config.cleaning.cavity_lower_fraction,
        refinement_iterations=config.cleaning.refinement_iterations,
    )

    runs = _true_runs(cavity_good)
    candidates = []
    for block_start, block_stop in runs:
        sample_start = int(starts[block_start])
        sample_stop = int(stops[block_stop - 1])
        candidates.append(
            {
                "global_start_index": sample_start,
                "global_stop_index_exclusive": sample_stop,
                "global_start_time_s": sample_start / fs,
                "global_stop_time_s": sample_stop / fs,
                "duration_s": (sample_stop - sample_start) / fs,
                "duration_h": (sample_stop - sample_start) / fs / 3600.0,
            }
        )

    candidates.sort(key=lambda item: item["duration_s"], reverse=True)

    def file_for_index(index: int) -> str | None:
        for source in files:
            lo = file_offsets[source]
            hi = lo + file_counts[source]
            if lo <= index < hi:
                return source.name
        if index == global_offset:
            return files[-1].name
        return None

    for item in candidates:
        item["starts_in_file"] = file_for_index(item["global_start_index"])
        stop_probe = max(
            item["global_start_index"],
            item["global_stop_index_exclusive"] - 1,
        )
        item["ends_in_file"] = file_for_index(stop_probe)

    top = candidates[: max(1, args.top)]
    summary = {
        "method": "approved cavity cleaning only; no eta0 calculation",
        "file_count": len(files),
        "duration_h": global_offset / fs / 3600.0,
        "sample_rate_hz": fs,
        "block_size_seconds": config.cleaning.block_size_seconds,
        "cavity_good_fraction": float(np.mean(cavity_good)),
        "candidate_count": len(candidates),
        "top_candidates": top,
        "selection_rule": (
            "Rank cavity-good segments by duration. Eta0 evaluation may stop "
            "once the best surviving final-good duration exceeds the duration "
            "of every unprocessed cavity candidate."
        ),
        "note": (
            "This tool does not reset or redefine the approved eta0 reference; "
            "it only identifies cavity-good candidate intervals."
        ),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
