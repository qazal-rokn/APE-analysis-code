#!/usr/bin/env python3
"""Synthetic-only existing-vendor single-bin diagnostics through stages 5,6,7.

No new axion signal model, no new detector noise, no original algorithm changes.
peak_norm=measured fitted background PS is strictly a unitless RELATIVE-EXCESS
parameterization for one-bin vendor software diagnostics. It is NOT the
Göttel-Raymond/LIGO beta, APE QWP power, an injection-recovery normalization,
or a physical coupling limit. Exact original vendor functions are called.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import uuid

import numpy as np
ROOT=Path.home()/"Desktop/APE_Research/APE_40MIN_FIVE_SIGNAL_SIMULATION"
S3=ROOT/"03_Native_C_LPSD"/"native_c"
S4=ROOT/"04_Background_Fitting"
VENDOR=Path.home()/"Downloads/Scalar-Dark-Matter-LPSD-1.0.0"
APE=Path.home()/"Desktop/APE_Research/APE_Automatic_Orchestrator"
FREQS=(3200,4000,4800,5900,6800)
MLESHA="031ad21dde68b08abb33efa18ac535ff9e7ea00dbc03cf13b6ec2f479328967b"
ULSHA="8f800fed30a0327763c87f983df7f507b44a879aa6ccfa8f57558146e1d26745"
STATSHA="d41ab695a61bf51951194f4e1edc316df81582f5eae256751472d5fba93f7f42"
SOURCE_SHA="919b38ef7e14e78b82066ce87688f0315fc87840cec6622954e7c964a854a6d7"
BINARY_SHA="cb0b1df10c9863145ba8b35df60e075eac3dc08965f3ad832857a7d4fec5eec9"
CSVFIELDS=("channel","rest_frequency_hz","frequency_hz","native_ps_rad2",
           "background_ps_rad2","log_ps_over_background","mu_hat_relative_power",
           "sigma_mu_relative_power","log_likelihood_null","log_likelihood_max",
           "q0_local","sqrt_q0_local_diagnostic","expected_mu_ul_relative_power_95",
           "expected_sqrt_mu_ul_relative_amplitude_95","status","error")

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as stream:
        for part in iter(lambda:stream.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()

def save_json(p,obj):
    with p.open("x") as file:
        json.dump(obj,file,indent=2,allow_nan=False)
        file.write("\n")

def csv_write(p,rows,fields):
    with p.open("x",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def load_rows(p):
    with p.open(newline="") as file:
        return list(csv.DictReader(file))

def native_freq(p):
    rows=[]
    for raw in p.read_bytes().splitlines():
        ln=raw.strip()
        if not ln or ln.startswith(b"#"):
            continue
        cols=ln.split()
        if len(cols)!=7:
            raise RuntimeError("Bad original native PS seven-column row: "+str(p))
        nums=[float(x) for x in cols]
        if not all(math.isfinite(x) for x in nums) or nums[2]<=0:
            raise RuntimeError("Nonfinite/nonpositive PS: "+str(p))
        rows.append(nums)
    return np.asarray(rows,dtype=float)

def load_vendor():
    paths={
        "mle":VENDOR/"scalar_dark_matter_lpsd/lvk/max_lkl_executable.py",
        "ul":VENDOR/"scalar_dark_matter_lpsd/lvk/calc_upper_limit.py",
        "stats":VENDOR/"scalar_dark_matter_lpsd/lvk/stats.py"}
    for name,expected in (("mle",MLESHA),("ul",ULSHA),("stats",STATSHA)):
        if sha(paths[name])!=expected:
            raise RuntimeError("Frozen original vendor "+name+" checksum mismatch")
    sys.path.insert(0,str(VENDOR))
    sys.path.insert(0,str(APE/"src"))
    from scalar_dark_matter_lpsd.lvk import max_lkl_executable as original_mle
    from scalar_dark_matter_lpsd.lvk import calc_upper_limit as original_ul
    from scalar_dark_matter_lpsd.lvk import stats as original_stats
    from ape_c_lpsd.analysis.ape_candidates import select_candidates
    for name,module in (("mle",original_mle),("ul",original_ul),("stats",original_stats)):
        if Path(module.__file__).resolve()!=paths[name].resolve():
            raise RuntimeError("Wrong imported original vendor module: "+name)
    return original_mle,original_ul,original_stats,select_candidates,{k:sha(v) for k,v in paths.items()}

def verify_sources():
    rep=json.loads((S3/"stage03_sparse_validation.json").read_text())
    if rep["status"]!="sparse_native_ps_validated" or rep["input_sha256"]!=SOURCE_SHA or rep["c_binary_sha256"]!=BINARY_SHA:
        raise RuntimeError("Stage 3 source is not the previously verified synthetic fixture")
    batch=json.loads((S4/"remaining_nine_1ccb351c7e"/"stage04_batch_summary.json").read_text())
    if batch["status"]!="ten_vendor_fits_pass" or batch["all_vendor_quality_pass"] is not True:
        raise RuntimeError("Stage 4 batch incomplete")
    return rep,batch

def stage4_folder(ch,f):
    if ch=="rho_x" and f==3200:
        return S4/"pilot_3200Hz_X_4f653750c3"/"vendor_fit_3ea3188d76"
    return S4/"remaining_nine_1ccb351c7e"/f"{'X' if ch=='rho_x' else 'Y'}_{f}Hz"

def import_fixed_input(ch,f):
    folder=stage4_folder(ch,f)
    report=json.loads((folder/"vendor_fit_report.json").read_text())
    if report["vendor_quality_pass"] is not True or report["input_units"]!="rad^2":
        raise RuntimeError("Bad saved Stage 4 vendor PS fit")
    rows=load_rows(folder/"native_ps_background_residuals.csv")
    if len(rows)!=1600:
        raise RuntimeError("Not 1600 saved original C bins")
    fs=np.asarray([float(x["frequency_hz"]) for x in rows])
    ps=np.asarray([float(x["native_ps_rad2"]) for x in rows])
    bg=np.asarray([float(x["background_ps_rad2"]) for x in rows])
    if not (np.diff(fs)>0).all() or not ((ps>0)&(bg>0)).all() or not np.isfinite(ps).all():
        raise RuntimeError("Native PS/background invalid")
    native=S3/ch/f"{ch}_{f}.lpsd"
    p=native_freq(native)
    if len(p)==0:
        raise RuntimeError("Original sparse native PS missing")
    if not (np.diff(p[:,0])>0).all():
        raise RuntimeError("Original sparse native frequency order")
    idx=np.searchsorted(fs,p[:,0])
    if np.any(idx>=len(fs)) or np.max(np.abs(fs[idx]-p[:,0]))>2e-6:
        raise RuntimeError("Stage 3/4 frequency grid disagreement")
    if not np.allclose(ps[idx],p[:,2],rtol=2e-9,atol=0):
        raise RuntimeError("Stage 3/4 PS numerical disagreement")
    pars=report["skewnormal_parameters"]
    model=np.array([pars["loc"],pars["scale"],pars["alpha"]],float)
    if not np.isfinite(model).all() or model[1]<=0:
        raise RuntimeError("Bad original vendor background model")
    return fs[idx],ps[idx],bg[idx],model,sha(folder/"vendor_fit_report.json"),sha(native),idx

def run():
    _,_,vendor_mle,candidate, hashes = (None,None,*load_vendor()[:1],None,None) if False else (None,None,None,None,None)
    mle,ul,stats,select,hashes=load_vendor()
    verify_sources()
    run_id="synthetic_single_bin_relative_excess_"+uuid.uuid4().hex[:10]
    root=ROOT/"05_Likelihood"/run_id
    if root.exists():
        raise FileExistsError(root)
    (root/"inputs").mkdir(parents=True,exist_ok=False)
    stage6=ROOT/"06_Candidates"/run_id
    stage7=ROOT/"07_Upper_Limits"/run_id
    for path in (stage6,stage7):
        path.mkdir(parents=True,exist_ok=False)
    rows=[]
    summary={}
    for f in FREQS:
        for ch in ("rho_x","rho_y"):
            print("VENDOR_ONE_BIN",ch,f,flush=True)
            freq,ps,bkg,pars,rep_sha,native_sha,idx=import_fixed_input(ch,f)
            subset=[]
            for i,(nu,obs,background) in enumerate(zip(freq,ps,bkg)):
                row={k:"" for k in CSVFIELDS}
                row.update(channel=ch,rest_frequency_hz=f,frequency_hz=float(nu),
                           native_ps_rad2=float(obs),background_ps_rad2=float(background),
                           log_ps_over_background=float(np.log(obs/background)),
                           status="pending")
                Y=np.log(np.array([[obs]],dtype=float))
                B=np.log(np.array([[background]],dtype=float))
                normal=np.array([[background]],dtype=float)
                shape=np.ones(1)
                model=pars.reshape(1,3)
                try:
                    l0,lmax,mu,sigma=mle.calc_max_lkl(
                        Y,B,model,normal,shape,
                        min_log10mu=-6,max_log10mu=2,n_seed_lkl=100)
                    if not np.isfinite([l0,lmax,mu]).all() or lmax<l0-1e-8:
                        raise ValueError("Original vendor optimizer returned invalid likelihood or mu")
                    if not np.isfinite(sigma) or sigma<=0:
                        raise ValueError("Original vendor fitted sigma invalid")
                    q0=max(0.0,2*(lmax-l0)) if mu>0 else 0.0
                    a,da=ul.get_upper_limits_approx(sigma,alpha_CL=0.95,tol=1e-10)
                    a=float(np.asarray(a).reshape(-1)[0])
                    if not math.isfinite(a) or a<=0:
                        raise ValueError("Original vendor expected limit invalid")
                    row.update(mu_hat_relative_power=float(mu),
                               sigma_mu_relative_power=float(sigma),
                               log_likelihood_null=float(l0),
                               log_likelihood_max=float(lmax),
                               q0_local=float(q0),
                               sqrt_q0_local_diagnostic=float(np.sqrt(q0)),
                               expected_mu_ul_relative_power_95=float(a*a),
                               expected_sqrt_mu_ul_relative_amplitude_95=float(a),
                               status="vendor_single_bin_relative_excess_diagnostic_PASS")
                except Exception as exc:
                    row.update(status="original_vendor_optimizer_or_UL_failed",
                               error=f"{type(exc).__name__}: {exc}")
                subset.append(row)
                rows.append(row)
            csv_write(root/"inputs"/f"{ch}_{f}_saved_native_and_vendor.csv",[
                {"frequency_hz":float(nu),"observed_native_ps_rad2":float(obs),
                 "saved_vendor_background_ps_rad2":float(bg0),
                 "model_loc":float(pars[0]),"model_scale":float(pars[1]),
                 "model_alpha":float(pars[2]),
                 "source_stage04_index":int(index)}
                for nu,obs,bg0,index in zip(freq,ps,bkg,idx)],
                ("frequency_hz","observed_native_ps_rad2","saved_vendor_background_ps_rad2",
                 "model_loc","model_scale","model_alpha","source_stage04_index"))
            save_json(root/"inputs"/f"{ch}_{f}_source_hashes.json",{
                "stage04_vendor_report_sha256":rep_sha,"stage03_native_sha256":native_sha,
                "source_dataset_hash":SOURCE_SHA,"vendor_python_source_sha256":hashes,
                "native_ps_units":"rad^2","spectral_grid_j":847300})
            csv_write(root/f"{ch}_{f}_likelihood.csv",subset,CSVFIELDS)
            summary[f"{ch}_{f}"]={"tested":len(subset),"passed":sum(r["status"].endswith("_PASS") for r in subset),
                                  "failed":sum(not r["status"].endswith("_PASS") for r in subset)}
    csv_write(root/"all_original_vendor_single_bin_results.csv",rows,CSVFIELDS)
    numeric=[r for r in rows if r["status"].endswith("_PASS")]
    save_json(root/"STAGE05_STATUS.json",{
        "stage":5,"status":"provisional_one_bin_relative_excess_vendor_diagnostic",
        "total_bins":len(rows),"numeric_fits_passed":len(numeric),"numeric_fit_failures":len(rows)-len(numeric),
        "per_window":summary,
        "scientific_interpretation":"This tests original vendor single-bin likelihood with peak_norm=background native PS; mu is relative excess, NOT axion signal strength or g coupling.",
        "one_bin_model_only":True,"original_vendor_functions_unchanged":True,
        "vendor_source_sha256":hashes,"injection_recovery_normalization_validated":False})
    for f in FREQS:
        for ch in ("rho_x","rho_y"):
            window=[r for r in rows if r["channel"]==ch and r["rest_frequency_hz"]==f]
            if all(r["status"].endswith("_PASS") for r in window):
                freqs=np.asarray([r["frequency_hz"] for r in window],float)
                zeros=np.asarray([r["log_likelihood_null"] for r in window],float)
                maximum=np.asarray([r["log_likelihood_max"] for r in window],float)
                mu=np.asarray([r["mu_hat_relative_power"] for r in window],float)
                result=select(freqs,zeros,maximum,mu,local_sigma_threshold=5.0,
                              n_effective_trials=None)
                result["status"]="uncalibrated_local_one_bin_relative_excess_candidates"
            else:
                result={"status":"SKIPPED_INCOMPLETE_VENDOR_MLE_OR_SIGMA",
                        "frequency_hz":f,"channel":ch,
                        "failed_bins":sum(not r["status"].endswith("_PASS") for r in window)}
            result.update(channel=ch,rest_frequency_hz=f)
            save_json(stage6/f"{ch}_{f}_candidates.json",result)
    save_json(stage6/"STAGE06_STATUS.json",{
        "stage":6,"status":"provisional_local_candidate_report",
        "candidate_rule":"original existing APE single-bin local q0 grouping, threshold 5, no trials correction",
        "global_significance_validated":False,"physical_detection_claim":False,
        "source_stage05_run_id":run_id,
        "failed_windows":sum(any(not r["status"].endswith("_PASS") for r in rows
                            if r["channel"]==ch and r["rest_frequency_hz"]==f)
                            for ch in ("rho_x","rho_y") for f in FREQS)})
    limits=[{
        "channel":r["channel"],"rest_frequency_hz":r["rest_frequency_hz"],
        "frequency_hz":r["frequency_hz"],
        "sigma_mu_relative_power":r["sigma_mu_relative_power"],
        "expected_mu_ul_relative_power_95":r["expected_mu_ul_relative_power_95"],
        "expected_sqrt_mu_ul_relative_amplitude_95":r["expected_sqrt_mu_ul_relative_amplitude_95"],
        "status":r["status"]} for r in rows]
    csv_write(stage7/"all_expected_relative_excess_limits.csv",limits,
              tuple(limits[0].keys()))
    save_json(stage7/"STAGE07_STATUS.json",{
        "stage":7,"status":"provisional_vendor_expected_relative_excess_limits_only",
        "units":"dimensionless excess power / saved fitted background PS",
        "vendor_ul_method":"original vendor calc_upper_limit.get_upper_limits_approx on original MLE sigma",
        "confidence_level":0.95,"successes":len(numeric),"failures":len(rows)-len(numeric),
        "coupling_gev_inv_computed":False,"physical_ape_exclusion":False,
        "global_significance_validated":False,"coverage_validated":False,
        "source_stage05_run_id":run_id})
    render(root,stage6,stage7)
    final=ROOT/"FINAL_SYNTHETIC_REVIEW"/run_id
    final.mkdir(parents=True,exist_ok=False)
    status={
        "fixture":"original_40minute_five_injection_synthetic_20261029",
        "stage_1":"isolated_original_APE_calibration_pass",
        "stage_2":"isolated_original_APE_cleaning_pass",
        "stage_3":"ten_native_C_sparse_PS_windows_pass",
        "stage_4":"ten_original_vendor_background_fits_pass",
        "stage_5":"provisional_vendor_single_bin_relative_excess_diagnostic",
        "stage_6":"local_candidates_uncalibrated_diagnostic",
        "stage_7":"relative_expected_excess_limits_only_NO_PHYSICAL_g_LIMIT",
        "scientific_blocker":"No verified APE QWP coupling response or independent native-C signal-power template, so physical injected g recovery, global 5sigma and observed g95 exclusions are unapproved.",
        "total_tested_bins":len(rows),"numeric_bins_passed":len(numeric),
        "original_vendor_unchanged":True,"synthetic_input_sha256":SOURCE_SHA,
        "vendor_hashes":hashes,"output_root":str(ROOT),"run_id":run_id,
    }
    save_json(final/"FINAL_PIPELINE_STATUS.json",status)
    script=Path(__file__).resolve()
    if script.is_file():
        shutil.copyfile(script,final/"analysis_and_replot.py")
    files=[]
    for directory in (root,stage6,stage7,final):
        for path in sorted(directory.rglob("*")):
            if path.is_file() and path.name!="inventory.csv":
                files.append({"relative_path":str(path.relative_to(ROOT)),
                              "sha256":sha(path),"size_bytes":path.stat().st_size})
    csv_write(final/"inventory.csv",files,("relative_path","sha256","size_bytes"))
    print(json.dumps({"status":"SYNTHETIC_STAGES_05_06_07_DIAGNOSTIC_OUTPUTS_SAVED",
                      "run_id":run_id,"numeric_pass":len(numeric),"numeric_failed":len(rows)-len(numeric),
                      "final_report":str(final/"FINAL_PIPELINE_STATUS.json"),
                      "physical_coupling_result":"NOT_CALCULATED_UNVERIFIED_RESPONSE"},indent=2),flush=True)

def render(root,stage6,stage7):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for folder in (root,stage6,stage7):
        (folder/"figures").mkdir(exist_ok=False)
    rows=load_rows(root/"all_original_vendor_single_bin_results.csv")
    for f in FREQS:
        for ch in ("rho_x","rho_y"):
            selection=[r for r in rows if r["channel"]==ch and int(r["rest_frequency_hz"])==f
                       and r["status"].endswith("_PASS")]
            if not selection:
                continue
            x=np.asarray([float(r["frequency_hz"])-f for r in selection])*1e3
            q=np.asarray([float(r["sqrt_q0_local_diagnostic"]) for r in selection])
            expected=np.asarray([float(r["expected_sqrt_mu_ul_relative_amplitude_95"]) for r in selection])
            for destination,y,label,stem in (
                (stage6,q,"Local sqrt(q0) (not global significance)","local_q0"),
                (stage7,expected,"Expected 95% relative-amplitude UL (not g)","expected_relative_ul")):
                fig,ax=plt.subplots(figsize=(7.5,4))
                ax.plot(x,y,marker=".",lw=1)
                ax.set_xlabel("Observed frequency minus rest frequency (mHz)")
                ax.set_ylabel(label)
                ax.set_title(f"Synthetic {f} Hz {ch}, vendor single-bin diagnostic")
                ax.grid(alpha=0.25)
                fig.tight_layout()
                name=f"{ch}_{f}_{stem}"
                fig.savefig(destination/"figures"/(name+".png"),dpi=180)
                fig.savefig(destination/"figures"/(name+".pdf"))
                plt.close(fig)
    # Replot from saved CSV only; this script keeps the plotting logic independent of fitting.
    replot=Path(root/"plot_saved_diagnostics.py")
    replot.write_text('''#!/usr/bin/env python3
"""Replot saved Stage 5–7 numerical CSV without running any vendor calculation."""
from pathlib import Path
import csv,sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
p=Path(sys.argv[1]).resolve()
out=Path(sys.argv[2]).resolve()
out.mkdir(parents=True,exist_ok=False)
with (p/"all_original_vendor_single_bin_results.csv").open() as f:
    rows=list(csv.DictReader(f))
for target in (3200,4000,4800,5900,6800):
    for ch in ("rho_x","rho_y"):
        data=[v for v in rows if v["channel"]==ch and int(v["rest_frequency_hz"])==target and v["status"].endswith("_PASS")]
        for key,ylabel in (("sqrt_q0_local_diagnostic","Local sqrt(q0), uncalibrated"),("expected_sqrt_mu_ul_relative_amplitude_95","Expected relative amplitude UL, not g")):
            if not data:continue
            fig,ax=plt.subplots(figsize=(8,4))
            ax.plot([(float(v["frequency_hz"])-target)*1000 for v in data],[float(v[key]) for v in data],marker=".")
            ax.set_xlabel("Frequency offset (mHz)")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=.3)
            fig.tight_layout()
            fig.savefig(out/f"{ch}_{target}_{key}.png",dpi=180)
            fig.savefig(out/f"{ch}_{target}_{key}.pdf")
            plt.close(fig)
''')

if __name__=="__main__":
    run()
