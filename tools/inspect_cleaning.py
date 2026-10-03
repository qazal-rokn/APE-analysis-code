#!/usr/bin/env python3
"""Read-only real-data inspection of the approved staged APE cleaning."""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ape_c_lpsd.cleaning import build_cleaning_reference_plan, clean_chunk_with_plan
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


@dataclass
class _IntervalAccumulator:
    """Streaming contiguous-interval statistics for a boolean target mask."""

    target_value: bool
    lengths_samples: list[int] = field(default_factory=list)
    longest_start_index: int | None = None
    longest_stop_index: int | None = None
    _open_start_index: int | None = None
    _open_length: int = 0
    _expected_next_index: int = 0

    def _close_open(self) -> None:
        if self._open_start_index is None or self._open_length <= 0:
            return
        start = self._open_start_index
        stop = start + self._open_length
        self.lengths_samples.append(self._open_length)
        longest = (
            0
            if self.longest_start_index is None or self.longest_stop_index is None
            else self.longest_stop_index - self.longest_start_index
        )
        if self._open_length > longest:
            self.longest_start_index = start
            self.longest_stop_index = stop
        self._open_start_index = None
        self._open_length = 0

    def feed(self, mask: np.ndarray, start_index: int) -> None:
        values = np.asarray(mask, dtype=bool)
        if values.ndim != 1:
            raise ValueError("interval mask must be one-dimensional")
        if start_index != self._expected_next_index:
            raise ValueError("interval chunks must be contiguous and ordered")

        target = values == self.target_value
        if target.size == 0:
            return

        padded = np.concatenate(([False], target, [False]))
        edges = np.flatnonzero(padded[1:] != padded[:-1])
        local_runs = list(zip(edges[::2], edges[1::2]))

        if self._open_start_index is not None:
            if local_runs and local_runs[0][0] == 0:
                first_start, first_stop = local_runs.pop(0)
                self._open_length += int(first_stop - first_start)
                if first_stop < target.size:
                    self._close_open()
            else:
                self._close_open()

        for run_start, run_stop in local_runs:
            absolute_start = start_index + int(run_start)
            run_length = int(run_stop - run_start)
            if run_stop == target.size:
                self._open_start_index = absolute_start
                self._open_length = run_length
            else:
                self.lengths_samples.append(run_length)
                longest = (
                    0
                    if self.longest_start_index is None or self.longest_stop_index is None
                    else self.longest_stop_index - self.longest_start_index
                )
                if run_length > longest:
                    self.longest_start_index = absolute_start
                    self.longest_stop_index = absolute_start + run_length

        self._expected_next_index = start_index + values.size

    def finish(self) -> None:
        self._close_open()

    def summary(self, sample_rate_hz: float) -> dict[str, object]:
        self.finish()
        lengths = np.asarray(self.lengths_samples, dtype=np.int64)
        if lengths.size == 0:
            duration_stats = {
                "count": 0,
                "total_duration_s": 0.0,
                "min_duration_s": None,
                "median_duration_s": None,
                "p90_duration_s": None,
                "p99_duration_s": None,
                "max_duration_s": None,
            }
        else:
            seconds = lengths.astype(np.float64) / sample_rate_hz
            duration_stats = {
                "count": int(lengths.size),
                "total_duration_s": float(seconds.sum()),
                "min_duration_s": float(seconds.min()),
                "median_duration_s": float(np.median(seconds)),
                "p90_duration_s": float(np.quantile(seconds, 0.90)),
                "p99_duration_s": float(np.quantile(seconds, 0.99)),
                "max_duration_s": float(seconds.max()),
            }

        if self.longest_start_index is not None and self.longest_stop_index is not None:
            duration_stats["longest_start_index"] = int(self.longest_start_index)
            duration_stats["longest_stop_index_exclusive"] = int(self.longest_stop_index)
            duration_stats["longest_start_time_s"] = (
                self.longest_start_index / sample_rate_hz
            )
            duration_stats["longest_stop_time_s"] = (
                self.longest_stop_index / sample_rate_hz
            )
        else:
            duration_stats["longest_start_index"] = None
            duration_stats["longest_stop_index_exclusive"] = None
            duration_stats["longest_start_time_s"] = None
            duration_stats["longest_stop_time_s"] = None
        return duration_stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect the approved cavity-first, eta0-second cleaning on one real "
            "Dewesoft file. Read-only: prints summary only and writes no data."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument(
        "--chunk-size-samples",
        type=int,
        default=1_000_000,
        help="raw samples per Dewesoft chunk (default: 1000000)",
    )
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)

    def chunk_factory():
        return iter_ape_dxd_chunks(
            config.input_file,
            channel_names=config.channels.as_mapping(),
            wrapper_path=dewesoft.wrapper_path,
            library_path=dewesoft.library_path,
            chunk_size_samples=args.chunk_size_samples,
        )

    plan = build_cleaning_reference_plan(
        chunk_factory,
        config.cleaning,
        config.calibration,
    )

    total = 0
    cavity_good_count = 0
    eta_good_count = 0
    final_good_count = 0
    cavity_bad_runs = _IntervalAccumulator(target_value=False)
    eta_only_bad_runs = _IntervalAccumulator(target_value=True)
    final_bad_runs = _IntervalAccumulator(target_value=False)
    final_good_runs = _IntervalAccumulator(target_value=True)
    cavity_block_medians: list[float] = []
    eta_block_medians: list[float] = []

    block = plan.block_samples
    for chunk in chunk_factory():
        result = clean_chunk_with_plan(
            chunk,
            plan,
            config.cleaning,
            config.calibration,
        )
        total += chunk.n_samples
        cavity_good_count += int(np.count_nonzero(result.cavity_good))
        eta_good_count += int(np.count_nonzero(result.eta_good))
        final_good_count += int(np.count_nonzero(result.final_good))

        cavity_bad_runs.feed(result.cavity_good, chunk.start_index)
        eta_only_bad = result.cavity_good & ~result.eta_good
        eta_only_bad_runs.feed(eta_only_bad, chunk.start_index)
        final_bad_runs.feed(result.final_good, chunk.start_index)
        final_good_runs.feed(result.final_good, chunk.start_index)

        v = np.asarray(chunk.cavity_transmission, dtype=np.float64)
        e = np.asarray(result.eta0_preliminary, dtype=np.float64)
        for start in range(0, chunk.n_samples, block):
            stop = min(chunk.n_samples, start + block)
            vb = v[start:stop]
            vb = vb[np.isfinite(vb)]
            if vb.size:
                cavity_block_medians.append(float(np.median(vb)))
            eb = e[start:stop]
            eb = eb[np.isfinite(eb)]
            if eb.size:
                eta_block_medians.append(float(np.median(eb)))

    if total != plan.total_samples:
        raise RuntimeError(
            f"inspection sample count changed: plan={plan.total_samples}, observed={total}"
        )

    cavity_bad_runs.finish()
    eta_only_bad_runs.finish()
    final_bad_runs.finish()
    final_good_runs.finish()

    summary = {
        "source_file": str(config.input_file),
        "n_samples": total,
        "sample_rate_hz": plan.sample_rate_hz,
        "cleaning": {
            "block_size_seconds": config.cleaning.block_size_seconds,
            "reference_window_seconds": config.cleaning.reference_window_seconds,
            "cavity_lower_fraction": config.cleaning.cavity_lower_fraction,
            "eta0_lower_fraction": config.cleaning.eta0_lower_fraction,
            "one_sided_lower_cuts_only": True,
        },
        "fractions": {
            "cavity_good_of_total": cavity_good_count / total,
            "eta0_good_of_total": eta_good_count / total,
            "eta0_good_given_cavity_good": (
                None if cavity_good_count == 0 else eta_good_count / cavity_good_count
            ),
            "final_good_of_total": final_good_count / total,
        },
        "counts": {
            "cavity_good_samples": cavity_good_count,
            "eta0_good_samples": eta_good_count,
            "final_good_samples": final_good_count,
            "cavity_bad_intervals": len(cavity_bad_runs.lengths_samples),
            "eta0_only_bad_intervals": len(eta_only_bad_runs.lengths_samples),
            "final_bad_intervals": len(final_bad_runs.lengths_samples),
        },
        "diagnostics": {
            "median_cavity_transmission_v_from_0p05s_block_medians": (
                None if not cavity_block_medians else float(np.median(cavity_block_medians))
            ),
            "median_eta0_from_cavity_good_0p05s_block_medians": (
                None if not eta_block_medians else float(np.median(eta_block_medians))
            ),
        },
        "interval_diagnostics": {
            "cavity_bad": cavity_bad_runs.summary(plan.sample_rate_hz),
            "eta0_only_bad": eta_only_bad_runs.summary(plan.sample_rate_hz),
            "final_bad": final_bad_runs.summary(plan.sample_rate_hz),
            "final_good_segments": final_good_runs.summary(plan.sample_rate_hz),
        },
        "memory_behavior": {
            "raw_full_record_retained": False,
            "full_record_references_retained_at_block_resolution": True,
            "writes_output_files": False,
        },
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
