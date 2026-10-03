#!/usr/bin/env python3
"""Read-only real-data inspection of the approved staged APE cleaning."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ape_c_lpsd.cleaning import build_cleaning_reference_plan, clean_chunk_with_plan
from ape_c_lpsd.ingestion import iter_ape_dxd_chunks, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def _count_bad_intervals(mask: np.ndarray, previous_good: bool) -> tuple[int, bool]:
    bad = ~np.asarray(mask, dtype=bool)
    if bad.size == 0:
        return 0, previous_good
    starts = int(np.count_nonzero(bad & np.concatenate(([previous_good], ~bad[:-1]))))
    return starts, bool(mask[-1])


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
    cavity_bad_intervals = 0
    eta_bad_intervals = 0
    final_bad_intervals = 0
    previous_cavity_good = True
    previous_eta_good = True
    previous_final_good = True
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

        added, previous_cavity_good = _count_bad_intervals(
            result.cavity_good, previous_cavity_good
        )
        cavity_bad_intervals += added
        added, previous_eta_good = _count_bad_intervals(
            result.eta_good, previous_eta_good
        )
        eta_bad_intervals += added
        added, previous_final_good = _count_bad_intervals(
            result.final_good, previous_final_good
        )
        final_bad_intervals += added

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
            "cavity_bad_intervals": cavity_bad_intervals,
            "eta0_bad_intervals": eta_bad_intervals,
            "final_bad_intervals": final_bad_intervals,
        },
        "diagnostics": {
            "median_cavity_transmission_v_from_0p05s_block_medians": (
                None if not cavity_block_medians else float(np.median(cavity_block_medians))
            ),
            "median_eta0_from_cavity_good_0p05s_block_medians": (
                None if not eta_block_medians else float(np.median(eta_block_medians))
            ),
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
