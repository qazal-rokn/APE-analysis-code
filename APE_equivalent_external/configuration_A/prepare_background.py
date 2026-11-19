"""Configuration A: sum X/Y linear PSD, then run UNCHANGED vendor background fit.

Uses the existing sum adapter and original
scalar_dark_matter_lpsd.lvk.fit_background main (subprocess).
Never alters vendor source or overwrites previously fitted X/Y backgrounds.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

from configuration_A.combine_quadrature_psd import combine_psd_hdf5


def run_combined_background(
    x_h5, y_h5, output_dir, *,
    ana_fmin, ana_fmax, n_processes=1, segment_size=1000,
):
    out = Path(output_dir).expanduser().resolve()
    if not (0 < ana_fmin < ana_fmax):
        raise ValueError("Require positive analysis range with fmin < fmax")
    if n_processes < 1 or segment_size < 5:
        raise ValueError("Invalid original vendor background fitter settings")
    if (out / "background_Q.json").exists():
        raise FileExistsError("Refusing to overwrite existing Q background")
    summary = combine_psd_hdf5(x_h5, y_h5, out)
    spectrum = Path(summary["output_file"])
    background = out / "background_Q.json"
    command = [
        sys.executable, "-m", "scalar_dark_matter_lpsd.lvk.fit_background",
        "--data-path", str(spectrum),
        "--output-path", str(background),
        "--ana-fmin", str(ana_fmin),
        "--ana-fmax", str(ana_fmax),
        "--n-processes", str(n_processes),
        "--segment-size", str(segment_size),
    ]
    # Original vendor program performs ALL spline/background distribution fits.
    # The combined PSD is retained if the fit fails, so recovery is explicit.
    subprocess.run(command, check=True)
    if not background.is_file():
        raise RuntimeError("Original vendor background fit did not create its JSON")
    return {"combined_psd": str(spectrum), "fitted_background": str(background),
            "vendor_background_command": command, "configuration": "A_Q"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--x-h5", required=True)
    p.add_argument("--y-h5", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--ana-fmin", type=float, required=True)
    p.add_argument("--ana-fmax", type=float, required=True)
    p.add_argument("--n-processes", type=int, default=1)
    p.add_argument("--segment-size", type=int, default=1000)
    a = p.parse_args()
    print(json.dumps(run_combined_background(
        a.x_h5, a.y_h5, a.output_dir, ana_fmin=a.ana_fmin,
        ana_fmax=a.ana_fmax, n_processes=a.n_processes,
        segment_size=a.segment_size), indent=2))


if __name__ == "__main__":
    main()
