"""PORTABLE accelerated, pilot-gated single-bin likelihood diagnostic for APE Q.

For peak_shape=[1], maximize exactly the SAME fitted skew-normal log-likelihood
as original lvk.stats.log_likelihood.  The residual skew-normal MODE is
independent of the observed PSD amplitude for each background segment.
This replaces original vendor numeric optimizer/Fisher by a mathematically
derived mode and analytic curvature; therefore NOT bitwise-original outputs.
Original source and original result files are never changed.

Mandatory gate compares all 20 real-data native bins against untouched
calc_max_lkl pilot before computing full 65--8000 Hz native-grid result.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
import numpy as np
import h5py
from scipy.optimize import brentq
from scipy.special import log_ndtr
from scipy.stats import skewnorm

# Portable exact-model equivalents of original vendor JSON reader and cubic spline.
# Only numerical search/Fisher implementation differs; original vendor source untouched.
import pandas as pd
from scipy.interpolate import CubicSpline


class LPSDJSONIO:
    def __init__(self, filename):
        self.data = json.loads(Path(filename).read_text())
    @staticmethod
    def get_label(filepath):
        parts = Path(filepath).stem.split("_")
        return "bkginfo_" + "_".join(parts[-3:])
    def get_df(self, name):
        return pd.DataFrame(json.loads(self.data[name]))


def bkg_model(x, knots):
    k = len(knots) // 2
    return CubicSpline(knots[:k], knots[k:], extrapolate=True)(x)


LOG_SQRT_2PI = 0.5 * np.log(2.0 * np.pi)


def skew_mode_and_curvature(loc, scale, alpha):
    """Mode and d²lnpdf/dr² for the fitted skew-normal, in log-PSD residuals."""
    a = float(alpha)
    if not (np.isfinite(scale) and scale > 0 and np.isfinite(a)):
        raise ValueError("Invalid original-vendor skew-normal parameters")
    def inverse_mills(t):
        return np.exp(-0.5*t*t - LOG_SQRT_2PI - log_ndtr(t))
    def score(z):
        return -z + a * inverse_mills(a*z)
    zm = 0.0 if a == 0 else brentq(score, -15.0, 15.0, xtol=1e-13)
    t = a*zm
    lam = inverse_mills(t)
    d2logpdf_dr2 = (-1.0 - a*a*lam*(t+lam)) / (scale*scale)
    return float(loc+scale*zm), float(d2logpdf_dr2)


def compute(freq, logp, beta, table):
    """All native bins; first matching background segment at shared endpoints."""
    starts = table["fmin"].to_numpy(float)
    ends = table["fmax"].to_numpy(float)
    if not np.all(np.diff(starts) > 0) or np.any(ends[:-1] > starts[1:]):
        raise RuntimeError("Background segments overlap or are unordered")
    out = np.full((len(freq), 5), np.nan, dtype=np.float64)
    out[:, 0] = freq
    used = np.zeros(len(freq), dtype=bool)
    # Fast contiguous ranges, per-segment scientific model unchanged.
    for i, row in enumerate(table.itertuples(index=False)):
        lo = int(np.searchsorted(freq, row.fmin, side="left" if i == 0 else "right"))
        hi = int(np.searchsorted(freq, row.fmax, side="right"))
        if hi <= lo:
            continue
        ff, yy, bbeta = freq[lo:hi], logp[lo:hi], beta[lo:hi]
        if not np.isfinite(row.chi_sqr) or row.chi_sqr >= 10:
            continue
        if not (np.isfinite(ff).all() and np.isfinite(yy).all()
                and np.isfinite(bbeta).all() and np.all(bbeta > 0)):
            raise RuntimeError(f"Invalid input in segment {i}")
        knots = np.r_[np.asarray(row.x_knots, float), np.asarray(row.y_knots, float)]
        b = np.asarray(bkg_model(np.log(ff), knots), float)
        mode, curvature = skew_mode_and_curvature(row.loc_skew, row.sigma_skew, row.alpha_skew)
        # y - ln(B+beta*mu_hat) = mode:
        # mu_hat = (exp(y-mode)-B)/beta; use expm1 for cancellation stability.
        mu = np.exp(b) * np.expm1(yy - b - mode) / bbeta
        l0 = skewnorm.logpdf(yy-b, row.alpha_skew, loc=row.loc_skew, scale=row.sigma_skew)
        lm = skewnorm.logpdf(mode, row.alpha_skew, loc=row.loc_skew, scale=row.sigma_skew)
        # Original Fisher definition: sigma_mu = sqrt(-1/L''(mu_hat))
        # L'' = logpdf''(mode) * [beta/exp(y-mode)]^2
        sigma = np.exp(yy-mode) / bbeta / np.sqrt(-curvature)
        out[lo:hi, 1] = l0
        out[lo:hi, 2] = lm
        out[lo:hi, 3] = mu
        out[lo:hi, 4] = sigma
        used[lo:hi] = True
    return out, used


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spectrum", required=True)
    p.add_argument("--background-json", required=True)
    p.add_argument("--response-npz", required=True)
    p.add_argument("--pilot-npy", required=True)
    p.add_argument("--output-npy", required=True)
    p.add_argument("--fmin", type=float, default=65)
    p.add_argument("--fmax", type=float, default=8000)
    p.add_argument("--max-mu-z", type=float, default=.07,
                   help="Max |mu_analytic-mu_vendor|/sigma_vendor allowed for every pilot bin")
    p.add_argument("--max-sigma-rel", type=float, default=.10,
                   help="Max relative sigma disagreement allowed for every pilot bin")
    args = p.parse_args()
    start_time = time.monotonic()
    output = Path(args.output_npy)
    if output.exists() or output.with_suffix(".json").exists():
        raise RuntimeError("Refusing overwrite of existing diagnostic output")
    if not output.parent.is_dir():
        raise RuntimeError("Output directory does not exist")

    with h5py.File(args.spectrum, "r") as h:
        f = np.asarray(h["frequency"])
        y = np.asarray(h["logPSD"])
    with np.load(args.response_npz, allow_pickle=False) as npz:
        response_freq = np.asarray(npz["frequency_hz"])
        beta = np.asarray(npz["peak_norm_q"])
        response_metadata = json.loads(str(npz["metadata_json"].item()))
    if not np.array_equal(f, response_freq) or len(y) != len(f) or not np.all(beta > 0):
        raise RuntimeError("Response and original Q spectrum do not match exactly")
    io = LPSDJSONIO(args.background_json)
    table = io.get_df(io.get_label(filepath=args.spectrum)).sort_values("fmin").reset_index(drop=True)

    # In the pilot, original vendor's background first-match convention is
    # intentionally reproduced at shared boundaries; pilot is around 5 kHz.
    pilot = np.load(args.pilot_npy, allow_pickle=False)
    if pilot.shape != (20, 5):
        raise RuntimeError("Expected original 20-row vendor pilot")
    ids = np.searchsorted(f, pilot[:, 0])
    if np.any(ids >= len(f)) or not np.array_equal(f[ids], pilot[:, 0]):
        raise RuntimeError("Original pilot frequency grid mismatch")
    check, ok = compute(f[ids], y[ids], beta[ids], table)
    if not np.all(ok):
        raise RuntimeError("Original pilot bins not all covered")
    mu_z = np.abs(check[:, 3] - pilot[:, 3]) / pilot[:, 4]
    sigma_rel = np.abs(check[:, 4]/pilot[:, 4] - 1.0)
    l0_diff = np.abs(check[:, 1] - pilot[:, 1])
    lm_diff = np.abs(check[:, 2] - pilot[:, 2])
    report = dict(pilot_bins=len(pilot),
                  max_mu_difference_in_vendor_sigma=float(np.nanmax(mu_z)),
                  max_sigma_relative_difference=float(np.nanmax(sigma_rel)),
                  max_abs_zero_logl_difference=float(np.nanmax(l0_diff)),
                  max_abs_max_logl_difference=float(np.nanmax(lm_diff)))
    print("MANDATORY original-vendor pilot comparison:", json.dumps(report, indent=2), flush=True)
    if (not np.isfinite(mu_z).all() or not np.isfinite(sigma_rel).all()
            or np.max(mu_z) > args.max_mu_z
            or np.max(sigma_rel) > args.max_sigma_rel
            or np.max(l0_diff) > 1e-6 or np.max(lm_diff) > 1e-6):
        raise RuntimeError("FAILED vendor-pilot comparison. Refusing accelerated full scan.")

    i0 = int(np.searchsorted(f, args.fmin, side="left"))
    i1 = int(np.searchsorted(f, args.fmax, side="left"))
    fast, usable = compute(f[i0:i1], y[i0:i1], beta[i0:i1], table)
    valid = np.isfinite(fast[:, 1:]).all(axis=1) & usable
    summary = dict(
        designation="ACCELERATED SINGLE-BIN DIAGNOSTIC, NOT ORIGINAL VENDOR OPTIMIZER",
        full_duration_seconds=15435,
        input_hdf5=args.spectrum,
        input_background=args.background_json,
        input_response=args.response_npz,
        pilot_reference=args.pilot_npy,
        signal_shape="peak_shape=[1], unchanged original skew-normal likelihood",
        method="Analytic skew-normal mode and exact log-likelihood Hessian",
        physical_response_status=response_metadata.get("response_status"),
        pilot_comparison=report,
        frequencies_total=int(len(fast)),
        frequencies_finite=int(np.count_nonzero(valid)),
        unfilled_frequency_bins=int(np.count_nonzero(~valid)),
        analytic_mu_abs_gt_1e8=int(np.count_nonzero(valid & (np.abs(fast[:, 3]) > 1e8))),
        wall_seconds=float(time.monotonic()-start_time),
        note="Not original numerical optimizer/Fisher. Only pilot-verified numerical equivalent where optimization is interior. Physical UL not certified.",
    )
    # Pilot-gated output uses separate name; never touches original chunks.
    np.save(output, fast)
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2)+"\n")
    print("Accelerated result:", json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
