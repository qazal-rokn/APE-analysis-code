#!/usr/bin/env python3
"""Synthetic-only dense-grid/benchmark study. Never opens the APE cache."""
import json, os, resource, subprocess, sys, time, uuid
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import h5py, numpy as np
from ape_c_lpsd.core import load_config, parse_output, sha256

def job(exe,w,x,j,n,iteration,tag):
 h=w/f"{tag}.h5"; out=w/f"{tag}-{j}-{n}-{iteration}.txt"; cf=out.with_suffix('.cfg')
 with h5py.File(h,'w') as f:f.create_dataset('/strain',data=x)
 cf.write_text('GNUTERM "APE"\nFORMAT "fDNRb"\nEND\n')
 a=[str(exe),'-d','-D','/strain','-f','128','-h','0','-i',str(h),'-o',str(out),'-b','0','-e',str(127/128),'-s','9','-t','36','-J',str(j),'-n',str(n),'-N',str(iteration),'-w','-2','-p','200','-r','1','-l','50','-a','100','-m','1','-x','1','-M','20']
 t=time.perf_counter(); r=subprocess.run(a,text=True,capture_output=True,env={'LPSDCFN':str(cf),**os.environ}); sec=time.perf_counter()-t
 if r.returncode: raise RuntimeError(r.stderr or r.stdout)
 return parse_output(out), {'Jdes':j,'bins_per_job':n,'iteration':iteration,'runtime_s':sec,'output_bytes':out.stat().st_size,'command':a,'output_sha256':sha256(out)}
def main():
 root=Path(__file__).resolve().parents[1]; cfg=load_config(root/'configs/run.yaml'); exe=Path(cfg['lpsd_executable']); outdir=Path(cfg['output_directory'])/'synthetic-benchmarks'/f"run-{uuid.uuid4().hex}"
 outdir.mkdir(parents=True,exist_ok=False); rng=np.random.default_rng(20260924); t=np.arange(128)/128; tones=[9,18,34]
 results=[]
 for f in tones:
  a,m=job(exe,outdir,np.sin(2*np.pi*f*t),16,16,0,f'tone-{f}'); peak=a[np.argmax(a[:,1])]; results.append({'tone_hz':f,'peak_hz':float(peak[0]),'peak_psd':float(peak[1]),'peak_asd':float(np.sqrt(peak[1])),'frequency_offset_hz':float(peak[0]-f),**m})
 noise=[]
 for seed in range(4):
  a,m=job(exe,outdir,np.random.default_rng(seed).normal(size=128),16,16,0,f'noise-{seed}'); noise.append({'seed':seed,'band_mean_psd':float(a[:,1].mean()),'band_std_psd':float(a[:,1].std()),**m})
 bench=[]
 for j in (4,8,16,32):
  _,m=job(exe,outdir,rng.normal(size=128),j,j,0,f'bench-{j}'); bench.append(m)
 batch_signal=rng.normal(size=128); full,_=job(exe,outdir,batch_signal,16,16,0,'batch-full'); b0,_=job(exe,outdir,batch_signal,16,8,0,'batch-split'); b1,_=job(exe,outdir,batch_signal,16,8,1,'batch-split')
 report={'synthetic_only':True,'tone_results':results,'noise_results':noise,'noise_mean_psd_mean':float(np.mean([x['band_mean_psd'] for x in noise])),'noise_mean_psd_std_across_realizations':float(np.std([x['band_mean_psd'] for x in noise],ddof=1)),'resource_benchmark':bench,'batch_equivalent':bool(np.allclose(full,np.vstack((b0,b1)))),'executable_sha256':sha256(exe),'run_directory':str(outdir),'config_snapshot':cfg,'no_real_cache_access':True}
 (outdir/'run_record.json').write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))
if __name__=='__main__':main()
