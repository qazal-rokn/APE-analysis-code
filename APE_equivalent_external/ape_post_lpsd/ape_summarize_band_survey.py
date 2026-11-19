#!/usr/bin/env python3
"""Consolidate the resumable APE 64-mass survey without re-running its fits.

Produces a diagnostic CSV, JSON audit and PNG plots. No physical upper limit,
discovery significance, trial-factor correction, or valid coverage is claimed.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import numpy as np


def read_mass(directory: Path, index: int, reference_frequency: float, config_hash: str):
    completed = directory / "completed.json"
    partial = directory / "partial.json"
    if completed.exists() and partial.exists():
        raise ValueError(f"Conflicting complete and partial markers: {directory}")
    marker = completed if completed.exists() else partial
    if not marker.exists():
        raise FileNotFoundError(f"Mass index {index} has no complete/partial marker")
    obj = json.loads(marker.read_text())
    if (obj.get("config_hash") != config_hash or obj.get("index") != index
            or not np.isclose(obj.get("rest_hz", np.nan), reference_frequency,
                              rtol=1e-13, atol=0)):
        raise ValueError(f"Mass index {index} provenance mismatch")
    status = "complete" if marker == completed else "partial_sigma"
    rows = []
    for channel in ("X", "Y"):
        fit = obj.get(channel)
        if not isinstance(fit, dict):
            raise ValueError(f"{directory}: {channel} fit missing")
        ll0, llmax, mu = (fit.get("lnL0"), fit.get("lnLmax"), fit.get("mu_hat"))
        if any(x is None or not np.isfinite(x) for x in (ll0, llmax, mu)):
            raise ValueError(f"{directory}: invalid {channel} likelihood")
        sigma = fit.get("sigma_mu")
        if sigma is not None and (not np.isfinite(sigma) or sigma <= 0):
            raise ValueError(f"{directory}: invalid {channel} sigma")
        if status == "complete" and sigma is None:
            raise ValueError(f"{directory}: complete checkpoint without {channel} sigma")
        # Strictly a nominal diagnostic. In the vendor adapter, the unconstrained
        # signed best fit is set to zero when negative for one-sided q0.
        q0 = max(0., 2. * (float(llmax) - float(ll0))) if mu > 0 else 0.
        rows.append(dict(index=index,rest_hz=float(reference_frequency),
            channel=channel,status=status,lnL0=float(ll0),lnLmax=float(llmax),
            mu_hat_vendor_signed=float(mu),
            sigma_mu_vendor_uncalibrated=float(sigma) if sigma is not None else None,
            q0_nominal_uncalibrated=float(q0)))
    return rows


def collect(root: Path):
    manifest = json.loads((root / "scan_manifest.json").read_text())
    frequencies = manifest["rest_frequencies_hz"]
    provenance = manifest["config_hash"]
    if (manifest.get("status") != "PROVISIONAL_SURVEY_NOT_PHYSICAL_LIMIT"
            or not isinstance(frequencies,list) or not frequencies):
        raise ValueError("Unexpected survey manifest")
    rows = []
    for i,freq in enumerate(frequencies):
        rows.extend(read_mass(root / f"mass_{i:07d}", i, float(freq), provenance))
    statuses = {i:set(row["status"] for row in rows if row["index"] == i)
                for i in range(len(frequencies))}
    complete = sum(v == {"complete"} for v in statuses.values())
    partial = sum(v == {"partial_sigma"} for v in statuses.values())
    if complete+partial!=len(frequencies):
        raise ValueError("Incomplete or contradictory survey statuses")
    return manifest, rows, complete, partial


def make_plot(rows, output_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs=plt.subplots(2,1,figsize=(11,8),sharex=True)
    for channel in ("X","Y"):
        data=[r for r in rows if r["channel"]==channel]
        f=np.array([r["rest_hz"] for r in data],float)
        mu=np.array([r["mu_hat_vendor_signed"] for r in data],float)
        sig=np.array([np.nan if r["sigma_mu_vendor_uncalibrated"] is None
                     else r["sigma_mu_vendor_uncalibrated"] for r in data],float)
        ax=axs[0]
        label=f"{channel}, signed vendor mu"
        ax.plot(f,mu,marker=".",linewidth=1,label=label)
        good=np.isfinite(sig)
        if np.any(good):
            ax.fill_between(f,mu-sig,mu+sig,where=good,alpha=.12,
                            label=f"{channel} ± nominal sigma")
        bad=~good
        if np.any(bad):
            ax.scatter(f[bad],mu[bad],s=38,facecolors="none",
                       edgecolors="black",label=f"{channel}: missing sigma")
        q0=np.array([r["q0_nominal_uncalibrated"] for r in data],float)
        axs[1].plot(f,q0,marker=".",linewidth=1,label=channel)
    for ax in axs:
        ax.set_xscale("log")
        ax.grid(alpha=.25)
        ax.legend(fontsize=8,ncol=2)
    axs[0].axhline(0,color="gray",linewidth=.8)
    axs[0].set_yscale("symlog",linthresh=1000)
    axs[0].set_ylabel("Signed mu (vendor, uncalibrated)")
    axs[1].set_ylabel("Nominal q0 (NOT significance)")
    axs[1].set_xlabel("Axion rest frequency (Hz)")
    fig.suptitle("APE 4.28 h: diagnostic composite likelihood survey\n"
                 "Uncalibrated correlated bins, provisional optical response")
    fig.tight_layout()
    fig.savefig(output_path,dpi=170)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--survey-dir",required=True,type=Path)
    p.add_argument("--output-dir",required=True,type=Path,
                   help="New directory; refuses existing output files")
    args=p.parse_args()
    root=args.survey_dir.expanduser().resolve()
    out=args.output_dir.expanduser().resolve()
    manifest, rows, complete, partial=collect(root)
    if out.exists():
        p.error("Output directory exists. Use a new name; do not overwrite an audit.")
    out.mkdir(parents=True,exist_ok=False)
    fields=list(rows[0])
    with (out/"ape_diagnostic_likelihood_64.csv").open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader();writer.writerows(rows)
    summary={"status":"PROVISIONAL_COMPOSITE_LIKELIHOOD_NOT_PHYSICAL_LIMIT",
             "frequencies_total":len(manifest["rest_frequencies_hz"]),
             "channels":["X","Y"],"complete_frequencies":complete,
             "partial_sigma_frequencies":partial,"fit_rows":len(rows),
             "source_survey":str(root),"config_hash":manifest["config_hash"],
             "warnings":["Vendor sigma and nominal q0 assume independent C bins",
                         "No look-elsewhere or coverage calibration",
                         "Unknown X/Y projection coefficients: do not combine",
                         "Optical response, halo and production PSLL provisional",
                         "64 frequencies do not constitute continuous mass coverage"]}
    (out/"ape_diagnostic_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    make_plot(rows,out/"ape_diagnostic_likelihood_64.png")
    print(json.dumps({**summary,"output_dir":str(out)},indent=2))


if __name__=="__main__":
    main()
