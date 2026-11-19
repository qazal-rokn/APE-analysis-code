"""One read-only Q-background QA job: local fits, peak absorption, injection recovery.

Uses original vendor bkg_model, pdf_skewnorm, stats.log_likelihood and
process_iteration.  The existing full fit is never overwritten.
Injection tests are synthetic single-bin power excesses; they are a software
diagnostic, not a physical-halo/coverage study.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO
from scalar_dark_matter_lpsd.lvk.fit_background import bkg_model, process_iteration
from scalar_dark_matter_lpsd.lvk.stats import pdf_skewnorm, log_likelihood


def eval_model(freq, row):
    knots = np.r_[np.asarray(row.x_knots, float), np.asarray(row.y_knots, float)]
    return np.asarray(bkg_model(np.log(freq), knots), float)


def scan_one(freq, logp, row):
    b = eval_model(freq, row)
    return b, logp - b


def single_bin_profile(y, b, model_args, *, beta, mu_grid):
    # Original vendor likelihood: no substitute residual PDF or statistic.
    Y = np.array([[float(y)]])
    B = np.array([[float(b)]])
    M = np.array([[float(beta)]])
    A = np.array([list(model_args)], dtype=float)
    shape = np.ones(1)
    values = np.array([log_likelihood(float(mu), Y, B, M, shape, A)
                       for mu in mu_grid])
    return float(mu_grid[np.nanargmax(values)])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spectrum", required=True)
    ap.add_argument("--background-json", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--refit", action="store_true",
                    help="Run ONE original-vendor background refit on an injected segment")
    args = ap.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    q = Path(args.spectrum)
    io = LPSDJSONIO(args.background_json)
    table = io.get_df(io.get_label(filepath=str(q))).sort_values("fmin").reset_index(drop=True)
    if len(table) == 0:
        raise ValueError("No fitted segments")
    with h5py.File(q, "r") as h:
        freq = np.asarray(h["frequency"])
        logp = np.asarray(h["logPSD"])
    good = np.isfinite(freq) & np.isfinite(logp) & (freq >= 65) & (freq <= 8000)
    if not np.any(good):
        raise RuntimeError("No valid original Q data")
    starts = table.fmin.to_numpy(float)
    ends = table.fmax.to_numpy(float)
    if np.any(ends[:-1] > starts[1:]):
        raise RuntimeError("Overlapping intervals: QA refuses ambiguous assignment")
    # Representative fitted windows across the entire frequency band.
    targets = [80, 120, 350, 1000, 3000, 5000, 7500]
    chosen = np.unique([int(np.argmin(abs((starts + ends)/2 - t))) for t in targets])
    def indices_for(row):
        i = np.searchsorted(freq, float(row.fmin), side="left")
        j = np.searchsorted(freq, float(row.fmax), side="right")
        return i, j

    fig, axes = plt.subplots(len(chosen), 2, figsize=(13, 2.75*len(chosen)))
    diagnostic = []
    for axrow, ix in zip(axes, chosen):
        row = table.iloc[int(ix)]
        i, j = indices_for(row)
        ff, yy = freq[i:j], logp[i:j]
        ok = np.isfinite(yy) & np.isfinite(ff) & (ff > 0)
        ff, yy = ff[ok], yy[ok]
        bg, residual = scan_one(ff, yy, row)
        a, b = axrow
        a.plot(ff, yy, ".", markersize=1.7, alpha=.32, label="log PSD")
        a.plot(ff, bg, "-", lw=1.7, label="original spline")
        a.set(title=f"{ff[0]:.1f}–{ff[-1]:.1f} Hz; chi2={float(row.chi_sqr):.2f}",
              xlabel="Frequency (Hz)", ylabel="ln(PSD)")
        a.legend(fontsize=7)
        b.hist(residual, bins=35, density=True, alpha=.55, label="residuals")
        xx = np.linspace(np.quantile(residual, .005), np.quantile(residual, .995), 250)
        b.plot(xx, pdf_skewnorm(xx, 1, float(row.loc_skew),
                                float(row.sigma_skew), float(row.alpha_skew)),
               lw=1.7, label="original skew-normal PDF")
        b.set(xlabel="ln(PSD/background)", ylabel="Density")
        b.legend(fontsize=7)
        diagnostic.append(dict(f_hz=float(np.median(ff)), chi2=float(row.chi_sqr),
                               n_bins=len(ff), median_log_residual=float(np.median(residual))))
    fig.tight_layout()
    fig.savefig(out/"01_segment_models_and_residual_pdfs.png", dpi=160)
    plt.close(fig)
    print("CHECK 1: plotted representative PSD fits and original skew-normal PDFs", flush=True)

    # Peak absorption: inspect ALL fitted segments, sampling full local windows.
    # Compare the spline level at strongest data peak with local background levels.
    # This is a flagging heuristic, NOT proof a real DM signal was absorbed.
    peak_records = []
    for ix, row in table.iterrows():
        i, j = indices_for(row)
        if j-i < 25:
            continue
        ff, yy = freq[i:j], logp[i:j]
        ok = np.isfinite(ff) & np.isfinite(yy) & (ff > 0)
        if np.count_nonzero(ok) < 25:
            continue
        ff, yy = ff[ok], yy[ok]
        bg = eval_model(ff, row)
        peak = int(np.argmax(yy))
        # Fitted-background excursion across local region.
        amp = float(bg[peak] - np.median(bg))
        peak_records.append(dict(segment=int(ix), frequency_hz=float(ff[peak]),
                                 log_excess_above_fit=float(yy[peak]-bg[peak]),
                                 spline_excursion_log=float(amp),
                                 peak_to_fit_power_ratio=float(np.exp(np.clip(yy[peak]-bg[peak], -700, 700))),
                                 chi2=float(row.chi_sqr)))
    peak_records.sort(key=lambda v: v["spline_excursion_log"], reverse=True)
    fig, ax = plt.subplots(figsize=(11, 4))
    p = peak_records[:30]
    ax.bar(np.arange(len(p)), [v["spline_excursion_log"] for v in p])
    ax.set(xlabel="Ranked fitted segments", ylabel="Spline peak excursion (natural log PSD)",
           title="Strongest local spline responses: review for potential peak absorption")
    fig.tight_layout()
    fig.savefig(out/"02_peak_absorption_flags.png", dpi=160)
    plt.close(fig)
    print("CHECK 2: scanned all fitted segments for peak-absorption flags", flush=True)

    # Controlled injections. Use original vendor log_likelihood, on the SAME
    # measured single-bin residual before and after adding known power.
    # Use beta = 0.5*B to keep a dimensionless test without claiming g units.
    # A paired delta-mu comparison cancels most pre-existing noise fluctuation.
    inj_results = []
    for target in [120, 1000, 5000]:
        ix = int(np.argmin(abs((starts+ends)/2 - target)))
        row = table.iloc[ix]
        i, j = indices_for(row)
        ff, yy = freq[i:j], logp[i:j]
        ok = np.isfinite(ff) & np.isfinite(yy) & (ff > 0)
        ff, yy = ff[ok], yy[ok]
        bg = eval_model(ff, row)
        # Choose typical residual near median, to avoid selecting a large line.
        kk = int(np.argmin(abs((yy-bg) - np.median(yy-bg))))
        beta = float(.5*np.exp(bg[kk]))
        observed = float(np.exp(yy[kk]))
        A = [float(row.loc_skew), float(row.sigma_skew), float(row.alpha_skew)]
        # Dense deterministic grid; original vendor likelihood values at every
        # point. This diagnostic is NOT the production vendor optimizer.
        grid = np.linspace(-1.9, 5.0, 850)
        null_mu = single_bin_profile(yy[kk], bg[kk], A, beta=beta, mu_grid=grid)
        injected_mu = 1.0
        new_y = float(np.log(observed + injected_mu*beta))
        found_mu = single_bin_profile(new_y, bg[kk], A, beta=beta, mu_grid=grid)
        inj_results.append(dict(frequency_hz=float(ff[kk]), injected_mu=1.,
                                baseline_best_mu=null_mu, injected_best_mu=found_mu,
                                recovered_delta_mu=found_mu-null_mu,
                                recovery_error=(found_mu-null_mu)-1.,
                                assumption="Frozen original background, single-bin additive power"))
    print("CHECK 3: three paired injections evaluated with original vendor likelihood", flush=True)

    refit = {"performed": False, "reason": "Optional original-vendor refit requires --refit"}
    if args.refit:
        # One *local* true vendor refit: how much of an injected narrow peak
        # is swallowed when the background is fit again?
        ix = int(np.argmin(abs((starts+ends)/2 - 5000)))
        row = table.iloc[ix]
        i, j = indices_for(row)
        i, j = max(0, i-50), min(len(freq), j+50)
        ff, yy = freq[i:j].copy(), logp[i:j].copy()
        ok = np.isfinite(ff) & np.isfinite(yy) & (ff > 0)
        ff, yy = ff[ok], yy[ok]
        orig_bg = eval_model(ff, row)
        idx = int(np.argmin(abs(ff-5000)))
        injection_power = np.exp(orig_bg[idx])
        injected = yy.copy()
        injected[idx] = np.log(np.exp(yy[idx])+injection_power)
        # Same process_iteration and default fit hyperparameters as production.
        kwargs = dict(kernel_size=30, plot_mean=False, verbose=False,
                      plot_path=None, plot_prefix="QA_only",
                      k_min=4, k_pruning=1, k_plateau=3, bin_scaling=.15,
                      cfd_fraction=.1, disable_tqdm=True)
        _, fit_knots, fit_distribution, chi = process_iteration(
            (0, np.log(ff), injected, kwargs))
        if fit_knots is not None:
            new_bg = float(bkg_model(np.log([ff[idx]]), np.asarray(fit_knots))[0])
            old_bg = float(orig_bg[idx])
            refit = dict(performed=True, frequency_hz=float(ff[idx]),
                         injected_power_over_old_background=1.,
                         fitted_background_power_change_fraction=float(np.exp(new_bg-old_bg)-1),
                         chi2_new=float(chi),
                         note="One selected segment only, original vendor fit; not general detection efficiency")
        else:
            refit = dict(performed=True, success=False, reason="Vendor local refit returned no model")
        print("CHECK 3b: original vendor fit recomputed once on one injected segment", flush=True)

    summary = dict(status="DIAGNOSTIC, not coverage validation",
                   source_spectrum=str(q), background_json=args.background_json,
                   original_segments=int(len(table)), local_segments=diagnostic,
                   top_peak_flags=peak_records[:30],
                   injection_recovery=inj_results,
                   injection_refit=refit,
                   interpretation="A passing software injection check does not establish false-alarm rates, halo-model recovery or broadband signal absorption")
    (out/"QA_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    fig, ax = plt.subplots(figsize=(7, 4))
    xs = [x["frequency_hz"] for x in inj_results]
    ys = [x["recovered_delta_mu"] for x in inj_results]
    ax.scatter(xs, ys, label="Paired vendor-likelihood recovery")
    ax.axhline(1, color="black", ls="--", label="Injected delta mu = 1")
    ax.set_xscale("log")
    ax.set(xlabel="Frequency (Hz)", ylabel="Recovered delta mu",
           title="Single-bin frozen-background injection test")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out/"03_injection_recovery.png", dpi=160)
    plt.close(fig)
    print(json.dumps({"plots": str(out), "injection_recovery": inj_results,
                      "single_refit": refit}, indent=2), flush=True)


if __name__ == "__main__":
    main()
