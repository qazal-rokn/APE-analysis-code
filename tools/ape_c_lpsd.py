#!/usr/bin/env python3
"""Run `python tools/ape_c_lpsd.py --help` for all commands."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from ape_c_lpsd.cli import main
main()

