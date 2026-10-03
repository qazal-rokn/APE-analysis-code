#!/usr/bin/env python3
"""Validate the five raw APE channels from one Dewesoft .dxd file."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ape_c_lpsd.ingestion import inspect_dxd, load_dxd_ingestion_config
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "APE stage-1 ingestion: validate X1f, Y1f, X2f, Y2f, and cavity "
            "transmission from a Dewesoft .dxd file without calibration."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    parser.add_argument(
        "--chunk-size-samples",
        type=int,
        default=1_000_000,
        help="samples per bounded-memory validation chunk (default: 1000000)",
    )
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)
    inspection = inspect_dxd(
        config.input_file,
        channel_names=config.channels.as_mapping(),
        wrapper_path=dewesoft.wrapper_path,
        library_path=dewesoft.library_path,
        chunk_size_samples=args.chunk_size_samples,
    )

    summary = {
        "source_file": str(inspection.source_file),
        "n_samples": inspection.n_samples,
        "duration_s": inspection.duration_s,
        "sample_rate_hz": inspection.sample_rate_hz,
        "start_time": (
            None if inspection.start_time is None else inspection.start_time.isoformat()
        ),
        "channel_names": inspection.channel_names,
        "channel_units": inspection.channel_units,
        "validation": {
            "mode": "full_file_chunked",
            "chunk_size_samples": inspection.chunk_size_samples,
            "chunk_count": inspection.chunk_count,
            "full_record_retained_in_memory": False,
        },
        "calibration_config_recorded_not_applied": {
            "cavity_transmission_w_per_v": config.calibration.cavity_transmission_w_per_v,
            "responsivity_a_per_w": config.calibration.responsivity_a_per_w,
            "transimpedance_gain_v_per_a": config.calibration.transimpedance_gain_v_per_a,
            "second_harmonic_pd_response": config.calibration.second_harmonic_pd_response,
        },
        "reader_metadata": inspection.metadata,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
