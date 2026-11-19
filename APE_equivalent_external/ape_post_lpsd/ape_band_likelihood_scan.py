#!/usr/bin/env python3
"""Resumable, selected-mass APE band-wide diagnostic likelihood survey.

Uses the committed halo, QWP and unchanged vendor-likelihood adapters. This
does NOT compute a physically calibrated full-native-bin upper limit.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

import numpy as np

HERE = Path(__file__).resolve().parent


def selected_masses(count, fmin, fmax):
    if count < 2 or not (65 < fmin < fmax < 10000):
        raise ValueError("Need count >= 2 and 65 < fmin < fmax < 10000 Hz")
    return np.geomspace(fmin, fmax, count)


def config_hash(masses, settings):
    raw = json.dumps({"frequencies": masses.tolist(), "settings": settings},
                     sort_keys=True)
    return hashlib.sha256(raw.encode()).hexdigest()


def classify_vendor_results(cases):
    """Classify complete, partial-sigma, and unusable vendor likelihood results.

    sigma_mu may be null when unchanged vendor FinDiff errors; preserve a
    numerically valid signed-mu fit instead of discarding it.
    """
    if (len(cases) != 2 or {row.get("channel") for row in cases} != {"X", "Y"}):
        return "INVALID"
    for row in cases:
        result = row.get("vendor_results", {})
        for field in ("lnL0", "lnLmax", "mu_hat"):
            value = result.get(field)
            if value is None or not np.isfinite(value):
                return "INVALID"
    return ("COMPLETE" if all(
        row.get("vendor_results", {}).get("sigma_mu") is not None and
        np.isfinite(row["vendor_results"]["sigma_mu"]) and
        row["vendor_results"]["sigma_mu"] > 0 for row in cases
    ) else "PARTIAL_SIGMA")


def execute(script, args, env):
    subprocess.run([sys.executable, str(HERE / script), *map(str, args)],
                   check=True, cwd=str(HERE.parent), env=env)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--count", type=int, default=64)
    p.add_argument("--fmin", type=float, default=65.1)
    p.add_argument("--fmax", type=float, default=8000.)
    p.add_argument("--start-index", type=int, default=0)
    p.add_argument("--stop-index", type=int, default=None)
    p.add_argument("--seeds", type=int, default=25)
    p.add_argument("--half-width-bins", type=int, default=40)
    p.add_argument("--halo-samples", type=int, default=120000)
    p.add_argument("--halo-modes", type=int, default=64)
    p.add_argument("--psll", type=float, default=238.13)
    p.add_argument("--finesse-s", type=float, default=1500.)
    p.add_argument("--tau-s", type=float, default=2.5e-9)
    p.add_argument("--tau1-s", type=float, default=1e-8)
    p.add_argument("--rho-dm-gev4", type=float, default=3.1e-42)
    p.add_argument("--unit-conversion", type=float, default=1.519e24)
    p.add_argument("--g-ref-gev-inv", type=float, default=1e-10)
    a = p.parse_args()
    try:
        masses = selected_masses(a.count, a.fmin, a.fmax)
    except ValueError as e:
        p.error(str(e))
    stop = len(masses) if a.stop_index is None else a.stop_index
    if not 0 <= a.start_index < stop <= len(masses):
        p.error("Invalid start/stop index")
    if not (3 <= a.seeds <= 1000 and 5 <= a.half_width_bins <= 5000
            and a.halo_samples >= 1000 and 8 <= a.halo_modes <= 512
            and 25 <= a.psll <= 250):
        p.error("Invalid numerical settings")
    optical = dict(finesse_s=a.finesse_s, tau_s=a.tau_s, tau1_s=a.tau1_s,
                   rho_dm_gev4=a.rho_dm_gev4,
                   unit_conversion=a.unit_conversion, g_ref_gev_inv=a.g_ref_gev_inv)
    if any(not np.isfinite(x) or x <= 0 for x in optical.values()):
        p.error("Optical parameters must be finite and positive")
    root = a.root.expanduser().resolve()
    out = a.output_dir.expanduser().resolve()
    x = root / "handoff" / "APE_0_15435_X.h5"
    y = root / "handoff" / "APE_0_15435_Y.h5"
    backgrounds = [root/"background"/f"background_{ch}_FULL.json" for ch in ("X","Y")]
    if any(not item.is_file() for item in [x,y,*backgrounds]):
        p.error("Missing X/Y spectra or original fitted backgrounds")
    settings = dict(root=str(root), seeds=a.seeds, width=a.half_width_bins,
        halo_samples=a.halo_samples, halo_modes=a.halo_modes,
        pilot_psll=a.psll, optical=optical)
    sha = config_hash(masses, settings)
    out.mkdir(parents=True, exist_ok=True)
    manifest = out/"scan_manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text()).get("config_hash") != sha:
            p.error("Output directory belongs to a different configuration")
    else:
        with manifest.open("x") as f:
            json.dump(dict(status="PROVISIONAL_SURVEY_NOT_PHYSICAL_LIMIT",
                config_hash=sha, settings=settings,
                rest_frequencies_hz=masses.tolist(),
                caveats=["unknown quadrature projections",
                  "uncalibrated overlapping-bin likelihood",
                  "unverified production Kaiser PSLL",
                  "provisional fixed-velocity halo and QWP model"]),f,indent=2)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(HERE.parent) + os.pathsep + env.get("PYTHONPATH","")
    passed, partial, failed, skipped = 0, [], [], 0
    for index in range(a.start_index, stop):
        rest = float(masses[index])
        location = out/f"mass_{index:07d}"
        location.mkdir(exist_ok=True)
        marker=location/"completed.json"
        if marker.exists():
            info=json.loads(marker.read_text())
            if info.get("config_hash")!=sha or info.get("rest_hz")!=rest:
                raise RuntimeError(f"Conflicting completed marker at {index}")
            skipped+=1
            continue
        halo=location/"halo.npz"
        qwp=location/"qwp.npz"
        likelihood=location/"likelihood.json"
        try:
            if not halo.exists():
                execute("ape_axion_halo_kernel_diagnostic.py",
                    ["--spectrum",x,"--native-kernel-py",HERE/"ape_native_kaiser_kernel.py",
                     "--rest-hz",rest,"--half-width-bins",a.half_width_bins,
                     "--halo-samples",a.halo_samples,"--halo-modes",a.halo_modes,
                     "--psll",a.psll,"--seed",20261119+index,
                     "--output",halo],env)
            if not qwp.exists():
                execute("ape_qwp_halo_psd_preview.py",
                    ["--halo-npz",halo,"--spectrum-x",x,"--spectrum-y",y,
                     "--output",qwp,"--finesse-s",a.finesse_s,
                     "--tau-s",a.tau_s,"--tau1-s",a.tau1_s,
                     "--rho-dm-gev4",a.rho_dm_gev4,
                     "--unit-conversion",a.unit_conversion,
                     "--g-ref-gev-inv",a.g_ref_gev_inv],env)
            if not likelihood.exists():
                execute("ape_provisional_qwp_multibin_vendor.py",
                    ["--root",root,"--qwp-npz",qwp,
                     "--seeds",a.seeds,"--output-json",likelihood],env)
            obj=json.loads(likelihood.read_text())
            cases=obj.get("diagnostic_cases",[])
            quality=classify_vendor_results(cases)
            if quality=="INVALID":
                raise RuntimeError("Vendor likelihood or mu_hat missing/nonfinite")
            results={row["channel"]:row["vendor_results"] for row in cases}
            if quality=="PARTIAL_SIGMA":
                # Do not silently replace the original vendor Fisher uncertainty.
                # Keep the finite likelihood fit as a checkpoint for later
                # uncertainty/coverage work, without labeling it complete.
                info=dict(status="PROVISIONAL_PARTIAL_SIGMA_UNAVAILABLE",
                    config_hash=sha,index=index,rest_hz=rest,
                    X=results["X"],Y=results["Y"])
                partial_file=location/"partial.json"
                if partial_file.exists():
                    old=json.loads(partial_file.read_text())
                    if old.get("config_hash")!=sha or old.get("rest_hz")!=rest:
                        raise RuntimeError("Existing partial checkpoint belongs to different scan")
                else:
                    with partial_file.open("x") as f:json.dump(info,f,indent=2)
                (location/"error.txt").unlink(missing_ok=True)
                partial.append(dict(index=index,f_hz=rest,
                                    missing_sigma_channels=[row["channel"] for row in cases
                                        if row["vendor_results"].get("sigma_mu") is None]))
                print(f"[{index+1}/{len(masses)}] {rest:.6f} Hz PARTIAL: finite vendor likelihood, sigma missing",flush=True)
                continue
            with marker.open("x") as f:
                json.dump(dict(status="PROVISIONAL_NUMERICAL_ONLY",
                    config_hash=sha, index=index, rest_hz=rest,
                    X=results["X"],Y=results["Y"]),f,indent=2)
            (location/"error.txt").unlink(missing_ok=True)
            passed+=1
            print(f"[{index+1}/{len(masses)}] {rest:.6f} Hz PASS X/Y",flush=True)
        except Exception as err:
            (location/"error.txt").write_text(traceback.format_exc())
            failed.append(dict(index=index,f_hz=rest,error=str(err)))
            print(f"[{index+1}/{len(masses)}] {rest:.6f} Hz FAILED: {err}",
                  file=sys.stderr,flush=True)
    report=dict(status="PROVISIONAL_SURVEY_NOT_PHYSICAL_LIMIT",
                selected_frequencies=len(masses), processed_index_range=[a.start_index,stop],
                passed=passed,partial=partial,skipped=skipped,failed=failed,output_dir=str(out))
    print(json.dumps(report,indent=2))
    if failed:
        raise SystemExit(1)

if __name__=="__main__":
    main()
