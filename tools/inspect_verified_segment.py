#!/usr/bin/env python3
"""Validate the approved segment without exporting it or invoking C-LPSD."""
from __future__ import annotations
import sys
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = ["/Users/user/Desktop/Code/Projects/APE-data-analysis/src", str(ROOT / "src")]
from ape_c_lpsd.verified_segment import Part, validate, validate_reference_equivalence

config = yaml.safe_load((ROOT / "configs/verified_segment.yaml").read_text())
parts = tuple(Part(Path(row["path"]), row["raw_start"], row["raw_stop"]) for row in config["segment"]["parts"])
validate_reference_equivalence(parts[1])
report = validate(parts)
print({"sample_count": report.sample_count, "chunks": report.chunks, "elapsed_seconds": report.elapsed_seconds, "part_counts": report.part_counts})
