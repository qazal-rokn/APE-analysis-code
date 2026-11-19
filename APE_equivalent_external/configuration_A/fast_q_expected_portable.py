"""Portable fast PRELIMINARY median-expected 95% UL from analytic Q likelihood.

Original vendor get_upper_limits_approx numerical formula, and its
1000-point sigma lookup/interpolation strategy. No observed UL.
Uses APE g_ref=1e-10 only as provisional response metadata convention.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Transcribed unchanged statistical equations from original vendor
# scalar_dark_matter_lpsd/lvk/stats.py + lvk/calc_upper_limit.py.
def cumf_q_mu_scalar(q_mu, mu, mu_prime, sigma):
    if q_mu > (mu / sigma) ** 2:
        return norm.cdf((q_mu - (mu ** 2 - 2 * mu * mu_prime) / sigma ** 2) / (2 * mu / sigma))
    return norm.cdf(np.sqrt(q_mu) - (mu - mu_prime) / sigma)


def get_upper_limits_approx(sigma, alpha_CL=.95, tol=1e-10):
    if np.isnan(sigma):
        return np.nan, np.nan
    def opt_pval(param, _sigma, n=1000):
        delta_q_mu = max(25, (param / _sigma) ** 2 + 1) / n
        q_mu = (param / _sigma) ** 2
        F_q_mu_0 = cumf_q_mu_scalar(q_mu, param, 0, _sigma)
        if F_q_mu_0 < 0.5:
            while F_q_mu_0 < 0.5:
                q_mu += delta_q_mu
                F_q_mu_0 = cumf_q_mu_scalar(q_mu, param, 0, _sigma)
        else:
            while F_q_mu_0 > 0.5:
                q_mu -= delta_q_mu
                F_q_mu_0 = cumf_q_mu_scalar(q_mu, param, 0, _sigma)
        p_val = 1. - cumf_q_mu_scalar(q_mu, param, param, _sigma)
        return (p_val - (1 - alpha_CL))**2
    popt = minimize(opt_pval, sigma, args=(sigma,), bounds=[(0, None)],
                    method="Nelder-Mead", tol=tol)
    if not popt.success:
        raise RuntimeError("Original upper-limit solver did not converge")
    upper_Lambda = np.sqrt(popt.x)
    sigma_Lambda = sigma / (2. * upper_Lambda)
    return upper_Lambda, sigma_Lambda


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--likelihood", required=True)
    p.add_argument("--g-ref", type=float, default=1e-10)
    args = p.parse_args()
    source = Path(args.likelihood)
    metadata = json.loads(source.with_suffix(".json").read_text())
    if metadata.get("designation") != "ACCELERATED SINGLE-BIN DIAGNOSTIC, NOT ORIGINAL VENDOR OPTIMIZER":
        raise RuntimeError("Only for pilot-validated accelerated likelihood diagnostic")
    data = np.load(source, mmap_mode="r")
    valid = np.isfinite(data[:,4]) & (data[:,4] > 0)
    sigma = np.asarray(data[valid,4])
    if not len(sigma):
        raise ValueError("No valid positive sigma")
    if np.max(sigma) == np.min(sigma):
        ul,unc = get_upper_limits_approx(float(sigma[0]))
        y = np.full(len(sigma), float(np.asarray(ul).item()))
        err = np.full(len(sigma), float(np.asarray(unc).item()))
    else:
        x = np.geomspace(float(np.min(sigma)),float(np.max(sigma)),1000)
        ul_vec = np.empty(len(x))
        unc_vec = np.empty(len(x))
        for i, s in enumerate(x):
            u, e = get_upper_limits_approx(float(s))
            ul_vec[i] = float(np.asarray(u).item())
            unc_vec[i] = float(np.asarray(e).item())
        y = interp1d(x, ul_vec, bounds_error=True)(sigma)
        err = interp1d(x, unc_vec, bounds_error=True)(sigma)
    result = np.column_stack([data[valid,0],y,err,args.g_ref*y,args.g_ref*err])
    out = source.with_name("Q_expected_UL_diagnostic.npy")
    if out.exists():
        raise FileExistsError(out)
    np.save(out, result)
    idx = np.unique(np.linspace(0,len(result)-1,min(5000,len(result)),dtype=int))
    np.savetxt(source.with_name("Q_expected_UL_plot_sample.csv"), result[idx],
               delimiter=",",
               header="frequency_hz,vendor_sqrt_mu_UL,vendor_sqrt_mu_uncertainty,provisional_g_UL_GeV_minus1,provisional_g_uncertainty",
               comments="")
    fig,ax = plt.subplots(figsize=(10,5))
    ax.loglog(result[idx,0],result[idx,3])
    ax.set(title="APE Q provisional expected 95% limit (diagnostic, NOT validated)",
           xlabel="Frequency (Hz)",ylabel="g_UL (GeV^-1)")
    ax.grid(alpha=.25,which="both")
    fig.tight_layout()
    fig.savefig(source.with_name("Q_expected_UL_diagnostic.png"),dpi=180)
    print("Saved",out, "and Q_expected_UL_diagnostic.png")


if __name__ == "__main__":
    main()
