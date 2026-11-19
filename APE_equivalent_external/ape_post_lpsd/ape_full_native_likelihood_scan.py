#!/usr/bin/env python3
"""Full native-grid *observed* APE likelihood scan, using the saved 4.28 h PSD.

No C-LPSD recomputation; no vendor background refitting; original vendor
stats.log_likelihood is called directly at every mass/amplitude trial.

This is a COMPOSITE-LIKELIHOOD engineering scan, not a statistically calibrated
physical coupling limit. X/Y each use a counterfactual unit-power projection;
they MUST NOT be combined until the measured optical projection is supplied.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import time

import h5py
import numpy as np
from scipy.optimize import minimize_scalar

from ape_post_lpsd.ape_axion_halo_kernel_diagnostic import halo_modes
from ape_post_lpsd.ape_native_kaiser_kernel import model_kernel
from ape_post_lpsd.physics.qwp_signal import qwp_amplitude_per_coupling
from scalar_dark_matter_lpsd.LPSDIO import LPSDJSONIO
from scalar_dark_matter_lpsd.lvk.fit_background import bkg_model
from scalar_dark_matter_lpsd.lvk.stats import log_likelihood

COLUMNS=("native_index","rest_hz","channel","status","lnL0",
         "lnLmax","mu_hat_signed","q0_nominal_uncalibrated",
         "max_T_over_B","n_bins","n_avg_at_center","background_segments")


def eligible_native_indices(freq, fmin, fmax, half_width):
    f=np.asarray(freq,float)
    valid=np.flatnonzero((f >= fmin) & (f <= fmax))
    if not len(valid):
        raise ValueError("No native frequencies in requested interval")
    return valid[(valid >= half_width) & (valid + half_width < len(f))]


def halo_benchmark(halo_samples, halo_modes_count, seed):
    modes,weights,report=halo_modes(rest=1.,n_samples=halo_samples,
          n_modes=halo_modes_count,seed=seed,
          sigma=229000*np.sqrt(3./2.),escape=528000.,lab_speed=229000.)
    if not np.isclose(weights.sum(),1.,atol=1e-12):
        raise ValueError("Halo weights do not sum to one")
    return modes, weights, report


def fitted_background_on_native_grid(table, freq, chi_cut=10.):
    """Map original fitted spline/skew-normal segments once per task window.

    Each native frequency receives parameters from the segment with latest
    fmin that contains it. Uncovered/overlap/failed-quality frequencies are
    marked invalid and never silently substituted or interpolated.
    """
    f=np.asarray(freq,float)
    required={"fmin","fmax","chi_sqr","x_knots","y_knots",
              "alpha_skew","loc_skew","sigma_skew"}
    if not required.issubset(table.columns):
        raise ValueError("Missing original vendor background fit fields")
    t=table.sort_values(["fmin","fmax"]).reset_index(drop=True)
    lower=t["fmin"].to_numpy(float)
    upper=t["fmax"].to_numpy(float)
    owners=np.searchsorted(lower,f,side="right")-1
    own=np.clip(owners,0,len(t)-1)
    valid=(owners>=0)&(f>=lower[own])&(f<=upper[own])
    # Reject overlapping interiors but allow exact shared endpoints. Check
    # neighboring segments; the original fit uses mostly adjacent ranges.
    prev=own-1
    overlapping=(prev>=0)&(f>lower[own])&(f<upper[np.maximum(prev,0)])
    valid &= ~overlapping
    logB=np.full(len(f),np.nan,float)
    pars=np.full((len(f),3),np.nan,float)
    for k in np.unique(own[valid]):
        fit=t.iloc[int(k)]
        chi=float(fit["chi_sqr"])
        indices=np.flatnonzero(valid & (own==k))
        if not (np.isfinite(chi) and chi<chi_cut):
            valid[indices]=False
            continue
        knots=np.r_[np.asarray(fit["x_knots"],float),
                    np.asarray(fit["y_knots"],float)]
        values=np.asarray(bkg_model(np.log(f[indices]),knots),float).reshape(-1)
        logB[indices]=values
        pars[indices]=[fit["loc_skew"],fit["sigma_skew"],fit["alpha_skew"]]
    valid &= np.isfinite(logB) & (np.exp(logB)>0)
    valid &= np.isfinite(pars).all(axis=1)&(pars[:,1]>0)
    return logB,pars,own,valid


def fit_signed_vendor_likelihood(observed_log_psd, logB, pars, T):
    """Bounded signed-mu optimization of unchanged vendor likelihood.

    Parameter s=mu*max(T/B) conditions the optimizer near unity. This is
    NOT the vendor Fisher sigma algorithm. The vendor likelihood *function*
    remains unchanged. Report a boundary hit rather than claiming success.
    """
    y=np.asarray(observed_log_psd,float).reshape(-1)
    b=np.asarray(logB,float).reshape(-1)
    pars=np.asarray(pars,float)
    T=np.asarray(T,float).reshape(-1)
    n=len(T)
    if (y.shape!=(n,) or b.shape!=(n,) or pars.shape!=(n,3)
            or not np.isfinite(y).all() or not np.isfinite(b).all()
            or not np.isfinite(pars).all() or (T<0).any()
            or not np.isfinite(T).all()):
        raise ValueError("Bad vendor likelihood inputs")
    ratio=float(np.max(T/np.exp(b)))
    if not (np.isfinite(ratio) and ratio>0):
        raise ValueError("Invalid signal/background scale")
    def ll(s):
        val=float(log_likelihood(float(s/ratio),y[:,None],b[:,None],
              T[:,None],np.ones(1,float),pars))
        return val if np.isfinite(val) else -np.inf
    l0=ll(0.)
    if not np.isfinite(l0):
        return dict(status="INVALID_L0",lnL0=None,lnLmax=None,
                    mu_hat_signed=None,q0_nominal_uncalibrated=None,
                    max_T_over_B=ratio)
    upper=10.
    status="FINITE"
    result=None
    for _ in range(3):
        # SciPy's bounded minimizer uses the ORIGINAL vendor likelihood.
        result=minimize_scalar(lambda s:-ll(s),method="bounded",
                  bounds=(-0.98,upper),
                  options={"xatol":2e-4,"maxiter":100})
        if not result.success or not np.isfinite(result.fun):
            status="OPTIMIZER_FAILED"
            break
        if result.x>=0.99*upper:
            upper*=10.
            status="UPPER_RANGE_EXTENDED"
            continue
        break
    if result is None or not result.success or not np.isfinite(result.fun):
        return dict(status="OPTIMIZER_FAILED",lnL0=l0,lnLmax=None,
                    mu_hat_signed=None,q0_nominal_uncalibrated=None,
                    max_T_over_B=ratio)
    if result.x>=0.99*upper:
        status="UPPER_RANGE_LIMIT"
    if result.x<=-0.975:
        status="LOWER_RANGE_LIMIT"
    lmax=max(l0,-float(result.fun))
    s_hat=float(result.x) if lmax>l0 else 0.
    return dict(status=status,lnL0=l0,lnLmax=lmax,
         mu_hat_signed=s_hat/ratio,
         q0_nominal_uncalibrated=(max(0.,2*(lmax-l0)) if s_hat>0 else 0.),
         max_T_over_B=ratio)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    p.add_argument("--task-id",type=int,default=None,
                   help="Default: SLURM_ARRAY_TASK_ID")
    p.add_argument("--task-bins",type=int,default=10000,
                   help="Native central frequencies per task")
    p.add_argument("--fmin",type=float,default=65.1)
    p.add_argument("--fmax",type=float,default=8000.)
    p.add_argument("--half-width-bins",type=int,default=40)
    p.add_argument("--halo-samples",type=int,default=120000)
    p.add_argument("--halo-modes",type=int,default=64)
    p.add_argument("--seed",type=int,default=20261119)
    p.add_argument("--psll",type=float,default=238.13)
    p.add_argument("--chi-cut",type=float,default=10.)
    p.add_argument("--finesse-s",type=float,default=1500.)
    p.add_argument("--tau-s",type=float,default=2.5e-9)
    p.add_argument("--tau1-s",type=float,default=1e-8)
    p.add_argument("--rho-dm-gev4",type=float,default=3.1e-42)
    p.add_argument("--unit-conversion",type=float,default=1.519e24)
    p.add_argument("--g-ref-gev-inv",type=float,default=1e-10)
    a=p.parse_args()
    if a.task_id is None:
        try:a.task_id=int(os.environ["SLURM_ARRAY_TASK_ID"])
        except (KeyError,ValueError):p.error("Provide --task-id or submit as Slurm array")
    if (a.task_id<0 or a.task_bins<1 or a.half_width_bins<5
          or a.half_width_bins>2000 or not 65<a.fmin<a.fmax<10000
          or a.halo_samples<1000 or not 8<=a.halo_modes<=512
          or not 25<=a.psll<=250 or a.chi_cut<=0):
        p.error("Invalid grid/halo/fit settings")
    opt=dict(finesse_s=a.finesse_s,tau_s=a.tau_s,tau1_s=a.tau1_s,
             rho_dm_gev4=a.rho_dm_gev4,
             unit_conversion=a.unit_conversion)
    if any(v<=0 or not np.isfinite(v) for v in [*opt.values(),a.g_ref_gev_inv]):
        p.error("Invalid provisional optical parameters")
    root=a.root.expanduser().resolve()
    out=a.outdir.expanduser().resolve()
    spectra={q:root/"handoff"/f"APE_0_15435_{q}.h5" for q in ("X","Y")}
    backgrounds={q:root/"background"/f"background_{q}_FULL.json"
                 for q in ("X","Y")}
    for f in [*spectra.values(),*backgrounds.values()]:
        if not f.is_file():p.error(f"Missing original four-hour input: {f}")
    modes_unit, weights, halo_report=halo_benchmark(
        a.halo_samples,a.halo_modes,a.seed)
    out.mkdir(parents=True,exist_ok=True)
    file=out/f"native_{a.task_id:06d}.npz"
    settings=dict(fmin=a.fmin,fmax=a.fmax,task_bins=a.task_bins,
        half_width_bins=a.half_width_bins,halo_samples=a.halo_samples,
        halo_modes=a.halo_modes,halo_seed=a.seed,psll=a.psll,
        chi_cut=a.chi_cut,optical=opt,g_ref_gev_inv=a.g_ref_gev_inv)
    if file.exists():
        with np.load(file,allow_pickle=False) as saved:
            meta=json.loads(str(saved["metadata_json"].item()))
            if (meta.get("root")!=str(root)
                    or meta.get("settings")!=settings
                    or meta.get("task_id")!=a.task_id
                    or meta.get("status")!="PROVISIONAL_FULL_NATIVE_GRID_COMPOSITE_LIKELIHOOD"):
                raise RuntimeError(f"Existing task output belongs to different scan: {file}")
        print(json.dumps(dict(status="SKIPPED_VERIFIED_EXISTING_TASK",
                              task_id=a.task_id,path=str(file))))
        return
    begin=time.monotonic()
    with h5py.File(spectra["X"],"r") as hx, h5py.File(spectra["Y"],"r") as hy:
        n=len(hx["frequency"])
        attrs=hx["logPSD"].attrs
        fs=float(attrs["fs"]);fmin0=float(attrs["fmin"])
        fmax0=float(attrs["fmax"]);jdes=int(attrs["Jdes"])
        if n!=jdes or len(hy["frequency"])!=n:
            raise ValueError("Native X/Y grid size mismatch")
        # Compute eligible central-bin indices without allocating full HDF5
        # frequency arrays. Verify the exact native values on each task.
        logstep=math.log(fmax0/fmin0)/(jdes-1)
        first=max(a.half_width_bins,
            int(math.ceil(math.log(a.fmin/fmin0)/logstep-1e-9)))
        last=min(n-1-a.half_width_bins,
            int(math.floor(math.log(a.fmax/fmin0)/logstep+1e-9)))
        total=max(0,last-first+1)
        n_tasks=(total+a.task_bins-1)//a.task_bins
        if a.task_id>=n_tasks:
            p.error(f"task-id {a.task_id} outside [0,{n_tasks-1}]")
        start=first+a.task_id*a.task_bins
        stop=min(first+(a.task_id+1)*a.task_bins,last+1)
        left=start-a.half_width_bins
        right=stop+a.half_width_bins
        f=np.asarray(hx["frequency"][left:right],float)
        fy=np.asarray(hy["frequency"][left:right],float)
        if not np.array_equal(f,fy):
            raise ValueError("Original X/Y native frequency arrays differ")
        expected=fmin0*np.exp(np.arange(left,right)*logstep)
        if not np.allclose(f,expected,rtol=1e-8,atol=1e-7):
            raise ValueError("Actual native frequency grid differs from C metadata")
        observed={q:np.asarray(h["logPSD"][left:right],float)
                  for q,h in [("X",hx),("Y",hy)]}
        avg={q:np.asarray(h["N_avg"][left:right],float)
             for q,h in [("X",hx),("Y",hy)]}
        fitgrid={}
        for ch in ("X","Y"):
            background=LPSDJSONIO(str(backgrounds[ch]))
            table=background.get_df(background.get_label(filepath=str(spectra[ch])))
            fitgrid[ch]=fitted_background_on_native_grid(table,f,a.chi_cut)
        result=np.full((stop-start,2,7),np.nan,float)
        # seven columns: lnL0, lnLmax, muhat, nominal q0, T/B, n_avg, status
        # Status: 0 FINITE, 1 optimizer/range warning, 2 fit unavailable.
        flags=np.full((stop-start,2),2,np.int8)
        indices=np.arange(start,stop,dtype=np.int64)
        for k,center in enumerate(indices):
            sl=slice(k,k+2*a.half_width_bins+1)
            rest=float(f[k+a.half_width_bins])
            freqs=f[sl]
            modes=rest*modes_unit
            R=qwp_amplitude_per_coupling(modes,**opt)
            K=model_kernel(freqs,modes,fs=fs,fmin=fmin0,fmax=fmax0,
                           jdes=jdes,psll=a.psll)
            T=.5*a.g_ref_gev_inv**2*(K@(weights*R**2))
            if not np.isfinite(T).all() or not np.any(T>0):
                continue
            for c,ch in enumerate(("X","Y")):
                logB,pars,owners,valid=fitgrid[ch]
                if not valid[sl].all() or not np.isfinite(observed[ch][sl]).all():
                    continue
                try:
                    z=fit_signed_vendor_likelihood(observed[ch][sl],
                                                   logB[sl],pars[sl],T)
                except (ValueError,FloatingPointError,ArithmeticError):
                    continue
                if z["lnLmax"] is None:
                    continue
                result[k,c,:]=[z["lnL0"],z["lnLmax"],z["mu_hat_signed"],
                    z["q0_nominal_uncalibrated"],z["max_T_over_B"],
                    avg[ch][k+a.half_width_bins],
                    float(len(np.unique(owners[sl])))]
                flags[k,c]=0 if z["status"]=="FINITE" else 1
            if (k+1)%1000==0:
                print(f"task={a.task_id} {k+1}/{stop-start}",flush=True)
        metadata=dict(status="PROVISIONAL_FULL_NATIVE_GRID_COMPOSITE_LIKELIHOOD",
            task_id=a.task_id,total_tasks=n_tasks,first_native_index=first,
            last_native_index=last,task_start=int(start),task_stop=int(stop),
            root=str(root),settings=settings,halo=halo_report,
            n_complete_X=int(np.sum(flags[:,0]==0)),
            n_complete_Y=int(np.sum(flags[:,1]==0)),
            n_warning_X=int(np.sum(flags[:,0]==1)),
            n_warning_Y=int(np.sum(flags[:,1]==1)),
            n_invalid_X=int(np.sum(flags[:,0]==2)),
            n_invalid_Y=int(np.sum(flags[:,1]==2)),
            seconds=time.monotonic()-begin,
            vendor_likelihood="unchanged stats.log_likelihood",
            optimizer="APE bounded scipy minimize_scalar, not vendor Fisher optimizer",
            notes=["Full 15435-second measured PSD and prefit backgrounds",
                   "Every native bin in chosen range is an axion rest hypothesis",
                   "Do not combine counterfactual full-power X/Y fits",
                   "Correlated native bins invalidate naive q0 and sigma",
                   "Production Kaiser PSLL, QWP response and halo unverified",
                   "No physical coupling upper limit or detection significance"])
        # Atomic task checkpoint (tmp file in same directory, no partial .npz).
        with tempfile.NamedTemporaryFile(dir=out,suffix=".npz",
                prefix=f".native_{a.task_id:06d}_",delete=False) as tmp:
            tmppath=Path(tmp.name)
        try:
            np.savez_compressed(tmppath,native_bin_indices=indices,
                rest_frequency_hz=f[a.half_width_bins:
                                   a.half_width_bins+(stop-start)],
                channel_values=result,flags=flags,
                field_names=np.array(["lnL0","lnLmax","mu_hat_signed",
                                      "q0_nominal_uncalibrated",
                                      "max_T_over_B","n_avg_center",
                                      "background_segment_count"]),
                metadata_json=np.array(json.dumps(metadata,sort_keys=True)))
            if file.exists():
                raise FileExistsError(file)
            tmppath.rename(file)
        finally:
            tmppath.unlink(missing_ok=True)
        print(json.dumps({k:metadata[k] for k in [
            "status","task_id","total_tasks","task_start","task_stop",
            "n_complete_X","n_complete_Y","n_invalid_X","n_invalid_Y",
            "seconds"]},indent=2))


if __name__=="__main__":
    main()
