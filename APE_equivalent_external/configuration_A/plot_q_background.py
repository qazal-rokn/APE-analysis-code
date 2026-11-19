"""Plot the actual 4.2875-hour Q PSD, original background fit and diagnostics.

Visualization only. Uses unchanged vendor LPSDJSONIO and bkg_model.
Never changes PSD, background JSON, inference results or vendor source.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO
from scalar_dark_matter_lpsd.lvk.fit_background import bkg_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", required=True)
    parser.add_argument("--background-json", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--fmin", type=float, default=65.0)
    parser.add_argument("--fmax", type=float, default=8000.0)
    parser.add_argument("--max-points", type=int, default=180000)
    args = parser.parse_args()

    source = Path(args.spectrum).expanduser().resolve()
    background = Path(args.background_json).expanduser().resolve()
    output = Path(args.output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.max_points < 1000:
        parser.error("--max-points must be at least 1000")

    with h5py.File(source, "r") as h:
        f_all = h["frequency"]
        # Constant-Q PSD is already calibrated and combined, in rad^2/Hz.
        i0 = int(np.searchsorted(f_all[:], args.fmin))
        i1 = int(np.searchsorted(f_all[:], args.fmax, side="right"))
        if i1 <= i0:
            raise ValueError("No frequencies in requested interval")
        indices = np.unique(np.linspace(i0, i1 - 1,
                                        min(i1-i0, args.max_points), dtype=np.int64))
        frequency = f_all[indices]
        psd = h["PSD"][indices]
        n_total = i1-i0

    io = LPSDJSONIO(str(background))
    df = io.get_df(io.get_label(filepath=str(source)))
    if len(df) == 0:
        raise RuntimeError("No fitted segments found")
    # Assign each sampled frequency to the first matching original segment.
    fitted_log = np.full(frequency.shape, np.nan)
    log_f = np.log(frequency)
    for _, segment in df.iterrows():
        lo, hi = float(segment["fmin"]), float(segment["fmax"])
        inside = ((frequency >= lo) & (frequency <= hi)
                  & ~np.isfinite(fitted_log))
        if np.any(inside):
            knots = np.r_[np.asarray(segment["x_knots"], dtype=float),
                          np.asarray(segment["y_knots"], dtype=float)]
            fitted_log[inside] = bkg_model(log_f[inside], knots)
    usable = np.isfinite(fitted_log) & np.isfinite(psd) & (psd > 0)
    if not np.all(usable):
        print(f"WARNING: {np.count_nonzero(~usable)} sampled points are not covered by valid fits")
    if not np.any(usable):
        raise RuntimeError("No valid data/background overlap")
    fitted_psd = np.exp(fitted_log[usable])
    residual = np.log(psd[usable]) - fitted_log[usable]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(frequency[usable], np.sqrt(psd[usable]), lw=0.5, alpha=.42,
            label="Measured summed Q ASD (4.2875 h)")
    ax.plot(frequency[usable], np.sqrt(fitted_psd), lw=1.0,
            label="Original vendor spline background")
    ax.set(xlabel="Frequency (Hz)", ylabel="Rotation ASD (rad/sqrt(Hz))",
           title="APE 4.2875-hour summed-quadrature spectrum and background")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.grid(True, which="both", alpha=.25)
    fig.tight_layout()
    fig.savefig(output/"Q_background_overlay.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(frequency[usable], residual, s=.3, alpha=.35, rasterized=True)
    ax.axhline(0, color="black", ls="--", lw=.8)
    ax.set(xlabel="Frequency (Hz)", ylabel="ln(PSD_Q) - ln(background)",
           title="Background-subtracted log-PSD residuals")
    ax.set_xscale("log")
    ax.grid(True, alpha=.2)
    fig.tight_layout()
    fig.savefig(output/"Q_log_residuals.png", dpi=180)
    plt.close(fig)

    chi = df["chi_sqr"].to_numpy(dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(residual[np.isfinite(residual)], bins=140, density=True)
    axes[0].set(xlabel="Log-PSD residual", ylabel="Density",
                title="Sampled residual distribution")
    centers = (df["fmin"].to_numpy(dtype=float) + df["fmax"].to_numpy(dtype=float))/2
    axes[1].scatter(centers, chi, s=2, alpha=.7)
    axes[1].axhline(10, ls="--", lw=1, color="black", label="Likelihood cutoff")
    axes[1].set_xscale("log")
    axes[1].set(xlabel="Frequency (Hz)", ylabel="Original fit chi-square",
                title="Background-fit quality")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(output/"Q_fit_diagnostics.png", dpi=180)
    plt.close(fig)

    print(f"Original background segments: {len(df)}")
    print(f"Segments passing chi2 < 10: {np.count_nonzero(np.isfinite(chi)&(chi<10))}")
    print(f"Analyzed native-frequency bins: {n_total}")
    print(f"Sampled points for plot: {len(indices)}")
    print(f"Median sampled log-residual: {np.median(residual):.4g}")
    print(f"Images saved to {output}")


if __name__ == "__main__":
    main()
