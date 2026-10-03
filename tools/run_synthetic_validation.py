#!/usr/bin/env python3
"""Controlled synthetic validation only; never reads the configured cache."""
import hashlib, json, os, subprocess, sys, tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import h5py, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ape_c_lpsd.core import load_config, parse_output, sha256, ValidationError

def run(exe, work, x, tag, n, batch):
    h=work/f"{tag}.h5"; out=work/f"{tag}-{batch}.txt"; cfg=work/f"{tag}-{batch}.cfg"
    with h5py.File(h,"w") as f: f.create_dataset("/strain",data=x)
    cfg.write_text('GNUTERM "APE"\nFORMAT "fDNRb"\nEND\n')
    argv=[str(exe),"-d","-D","/strain","-f","128","-h","0","-i",str(h),"-o",str(out),"-b","0","-e",str(127/128),"-s","9","-t","36","-J","4","-n",str(n),"-N",str(batch),"-w","-2","-p","200","-r","1","-l","50","-a","100","-m","1","-x","1","-M","20"]
    r=subprocess.run(argv,text=True,capture_output=True,env={"LPSDCFN":str(cfg),**os.environ})
    (out.with_suffix(".stdout.txt")).write_text(r.stdout); (out.with_suffix(".stderr.txt")).write_text(r.stderr)
    if r.returncode: raise RuntimeError(r.stderr or r.stdout)
    return parse_output(out), argv, out

def main():
    cfg=load_config(Path(__file__).resolve().parents[1]/"configs/run.yaml"); exe=Path(cfg["lpsd_executable"])
    if not exe.is_file(): raise SystemExit(f"missing executable: {exe}")
    with tempfile.TemporaryDirectory(prefix="ape-c-lpsd-validation-") as d:
        w=Path(d); t=np.arange(128)/128; tone=np.sin(2*np.pi*18*t); rng=np.random.default_rng(42); noise=rng.normal(size=128)
        full,argv,toneout=run(exe,w,tone,"rho_x_tone",4,0); b0,_,_=run(exe,w,tone,"rho_x_tone_b",2,0); b1,_,_=run(exe,w,tone,"rho_x_tone_b",2,1); joined=np.vstack((b0,b1))
        ny,_,nyout=run(exe,w,noise,"rho_y_noise",4,0)
        if not np.allclose(full,joined,rtol=1e-9,atol=1e-12): raise AssertionError("split C batches differ from complete batch")
        if len(np.unique(joined[:,0]))!=4 or not np.all(np.diff(joined[:,0])>0): raise AssertionError("batch coverage/order failure")
        peak=full[np.argmax(full[:,1]),0]
        if not (9 <= peak <= 36): raise AssertionError("tone peak outside C frequency region")
        if not (np.all(ny[:,1]>=0) and np.all(np.isfinite(ny)) and np.all(ny[:,2]>=1)): raise AssertionError("noise PSD/averages invalid")
        try: parse_output(w/"missing.txt")
        except FileNotFoundError: pass
        bad=w/"bad.txt"; bad.write_text("1 2 3\n")
        try: parse_output(bad); raise AssertionError("malformed output accepted")
        except ValidationError: pass
        fig,ax=plt.subplots(); ax.loglog(full[:,0],full[:,1],label="rho_x PSD"); ax.loglog(ny[:,0],ny[:,1],label="rho_y PSD"); ax.set(xlabel="Frequency (Hz)",ylabel="PSD (arb$^2$/Hz)"); ax.legend(); fig.savefig(w/"synthetic_psd.png"); plt.close(fig)
        fig,ax=plt.subplots(); ax.loglog(full[:,0],np.sqrt(full[:,1]),label="rho_x ASD"); ax.loglog(ny[:,0],np.sqrt(ny[:,1]),label="rho_y ASD"); ax.set(xlabel="Frequency (Hz)",ylabel="ASD (arb/sqrt(Hz))"); ax.legend(); fig.savefig(w/"synthetic_asd.png"); plt.close(fig)
        report={"synthetic_only":True,"tone_peak_region_hz":float(peak),"tone_psd":full[:,1].tolist(),"noise_psd":ny[:,1].tolist(),"multi_batch_equivalent":True,"two_channels_separate":True,"parser_malformed_rejected":True,"asd_is_sqrt_psd":bool(np.allclose(np.sqrt(full[:,1])**2,full[:,1])),"executable_sha256":sha256(exe),"effective_command":argv,"artifacts":{"tone_output_sha256":sha256(toneout),"noise_output_sha256":sha256(nyout)},"no_real_cache_access":True}
        print(json.dumps(report,indent=2))
if __name__=="__main__": main()
