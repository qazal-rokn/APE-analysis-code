#!/usr/bin/env python3
"""Read and validate the five raw APE channels from one Dewesoft .dxd file."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ape_c_lpsd.ingestion import load_dxd_ingestion_config, validate_ape_dxd_stream
from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "APE stage-1 ingestion: read exactly X1f, Y1f, X2f, Y2f, and "
            "cavity transmission from a Dewesoft .dxd file without calibration."
        )
    )
    parser.add_argument("--config", required=True, help="DXD ingestion YAML")
    args = parser.parse_args()

    config = load_dxd_ingestion_config(args.config)
    dewesoft = resolve_reference_dewesoft_paths(config.reference_project)
    data = validate_ape_dxd_stream(
        config.input_file,
        channel_names=config.channels.as_mapping(),
        wrapper_path=dewesoft.wrapper_path,
        library_path=dewesoft.library_path,
        chunk_samples=1_000_000,
    )

    summary = {
        "source_file": str(data.source_file),
        "n_samples": data.n_samples,
        "sample_rate_hz": data.sample_rate_hz,
        "duration_s": data.duration_s,
        "start_time": None if data.start_time is None else data.start_time.isoformat(),
        "channel_names": data.channel_names,
        "channel_units": data.channel_units,
        "stream_validation": {
            "validated_full_file": True,
            "chunk_samples": 1000000,
            "full_arrays_retained_in_memory": False,
        },
        "calibration_config_recorded_not_applied": {
            "cavity_transmission_w_per_v": config.calibration.cavity_transmission_w_per_v,
            "responsivity_a_per_w": config.calibration.responsivity_a_per_w,
            "transimpedance_gain_v_per_a": config.calibration.transimpedance_gain_v_per_a,
            "second_harmonic_pd_response": config.calibration.second_harmonic_pd_response,
        },
        "reader_metadata": data.metadata,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
