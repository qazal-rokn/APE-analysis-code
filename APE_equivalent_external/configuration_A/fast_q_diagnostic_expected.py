"""Write preliminary full-band EXPECTED UL and candidate screen from analytic QA result.

Uses original vendor upper-limit solver and original candidate screening adapter.
Input must be clearly labeled ANALYTIC DIAGNOSTIC, not vendor optimizer output.
Not a physically validated APE limit, not an observed UL.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ape_post_lpsd.ape_finalize_inference import original_expected_ul, original_candidate_screen


def main():
    root = Path("/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119/combined_Q/likelihood/fast_diagnostic")
    source = root/"Q_full_analytic_diagnostic.npy"
    report = json.loads(source.with_suffix(".json").read_text())
    if report.get("designation") != "ACCELERATED SINGLE-BIN DIAGNOSTIC, NOT ORIGINAL VENDOR OPTIMIZER":
        raise ValueError("Expected analytic diagnostic input only")
    data = np.load(source, allow_pickle=False)
    good = np.isfinite(data).all(axis=1) & (data[:, 4] > 0)
    if not np.any(good):
        raise ValueError("No finite sigma values")
    expected = original_expected_ul(data[good], g_ref=1e-10)
    out_npy = root/"Q_expected_UL_diagnostic.npy"
    if out_npy.exists():
        raise FileExistsError(out_npy)
    np.save(out_npy, expected)
    # Explicitly approximate visualization via fixed logarithmic samples.
    target = np.geomspace(expected[0,0], expected[-1,0], 2500)
    selected = np.unique(np.minimum(np.searchsorted(expected[:,0], target), len(expected)-1))
    np.savetxt(root/"Q_expected_UL_plot_sample.csv", expected[selected],
               delimiter=",", header="frequency_hz,vendor_sqrt_mu_UL,sqrt_mu_uncertainty,provisional_g_UL_GeV_minus1,g_uncertainty",
               comments="")
    candidates, screening = original_candidate_screen(data[good])
    np.savetxt(root/"Q_candidate_screen_diagnostic.csv", candidates,
               delimiter=",", header="frequency_hz,q0,signed_mu_hat,sqrt_mu_cluster",
               comments="")
    fig, ax = plt.subplots(figsize=(10,5))
    ax.loglog(expected[selected,0],expected[selected,3],lw=1)
    ax.set(xlabel="Frequency (Hz)", ylabel="Provisional expected 95% coupling UL (GeV$^{-1}$)",
           title="APE Q: single-bin analytic diagnostic (NOT physical validated limit)")
    ax.grid(alpha=.25, which="both")
    fig.tight_layout()
    fig.savefig(root/"Q_expected_UL_diagnostic.png",dpi=180)
    plt.close(fig)
    summary=dict(frequency_count=int(len(data)),
                 valid_frequencies=int(good.sum()),
                 expected_upper_limit_type="Original vendor median expected, 95% CL",
                 candidate_screening=screening,
                 source_fast_diagnostic=str(source),
                 physical_upper_limit_validated=False,
                 original_vendor_optimization_performed=False,
                 note="Diagnostic only: analytic acceleration must pass independent statistical/physical validation")
    (root/"Q_expected_UL_diagnostic_report.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2), flush=True)


if __name__=="__main__":
    main()
