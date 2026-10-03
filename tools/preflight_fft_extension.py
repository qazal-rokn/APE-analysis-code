#!/usr/bin/env python3
"""Report no-data FFT resource/range feasibility; run with --help."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ape_c_lpsd.fft_extension import load_fft_extension_config, preflight_full_record_fft

def main() -> None:
    parser = argparse.ArgumentParser(description="Preflight a future full-record FFT without reading a signal sequence.")
    parser.add_argument("--config", type=Path, default=ROOT / "configs/fft_extension.yaml")
    parser.add_argument("--sample-count", type=int, required=True, help="Validated contiguous sequence length supplied by the caller.")
    parser.add_argument("--sample-rate-hz", type=float, required=True, help="Validated sequence sampling rate supplied by the caller.")
    args = parser.parse_args()
    config = load_fft_extension_config(args.config)
    plan = preflight_full_record_fft(sample_count=args.sample_count, sample_rate_hz=args.sample_rate_hz, options=config.analysis)
    print(json.dumps({"method_identity": "fft_full_record", "expected_memory_bytes": plan.expected_memory_bytes, "memory_limit_bytes": plan.memory_limit_bytes, "bin_width_hz": plan.bin_width_hz, "fft_calculable_range_hz": [plan.fft_calculable_min_hz, plan.fft_calculable_max_hz], "configured_handoff_hz": plan.c_lpsd_handoff_hz, "signal_data_read": False}, indent=2))
if __name__ == "__main__":
    main()
