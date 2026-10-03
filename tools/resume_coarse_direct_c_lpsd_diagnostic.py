#!/usr/bin/env python3
"""Safely resume one interrupted COARSE DIAGNOSTIC ONLY run."""
from __future__ import annotations
import argparse, hashlib, json, os, resource, subprocess, sys, time, uuid
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, yaml
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from ape_c_lpsd.c_grid import c_frequency
from ape_c_lpsd.c_output_validation import validate_frequency_rows
LABEL='COARSE DIAGNOSTIC ONLY — NOT FOR CANDIDATES OR LIMITS'
def digest(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def expected_command(exe,meta,out,fs,n,s,g,b):
 return [str(exe),'-d','-D',meta['dataset'],'-f',repr(fs),'-h','0','-i',meta['path'],'-o',str(out),'-b','0','-e',repr((n-1)/fs),'-s',repr(g['fmin_hz']),'-t',repr(g['fmax_hz']),'-J',str(g['jdes']),'-n',str(g['bins_per_batch']),'-N',str(b),'-w',str(s['window']),'-p',str(s['psll_db']),'-r',str(s['linear_regression']),'-l',str(s['overlap_percent']),'-a',str(s['desired_averages']),'-m',str(s['minimum_averages']),'-x',str(s['input_scale']),'-M',str(s['max_memory_power'])]
def expected_frequency(g):
 return np.asarray([c_frequency(i,g['fmin_hz'],g['fmax_hz'],g['jdes']) for i in range(g['jdes'])])
def validate_saved(run,name,b,g,freq):
 text=run/f'{name}-batch-{b:03d}.txt';stdout=text.with_suffix('.stdout.txt');stderr=text.with_suffix('.stderr.txt');cf=text.with_suffix('.cfg')
 if not(text.is_file() and stdout.is_file() and stderr.is_file() and cf.is_file()):raise RuntimeError(f'missing completed artifact for {name} batch {b:03d}')
 rows,detail=validate_frequency_rows(text.read_text().splitlines(),freq[b*g['bins_per_batch']:(b+1)*g['bins_per_batch']])
 if f' -N {b} ' not in text.read_text():raise RuntimeError(f'provenance command does not identify batch {b:03d}')
 if 'FORMAT "fDNRb"' not in cf.read_text():raise RuntimeError(f'wrong output format config for {name} batch {b:03d}')
 return rows,{'artifact':text.name,'sha256':digest(text),'frequency_validation':detail,'stdout':stdout.name,'stderr':stderr.name,'config':cf.name,'revalidated_saved':True}
def run_new(run,name,b,exe,meta,fs,n,s,g,freq):
 attempt=uuid.uuid4().hex; out=run/f'{name}-batch-{b:03d}-retry-{attempt}.txt';cf=out.with_suffix('.cfg');cf.write_text('GNUTERM "APE coarse diagnostic resume"\nFORMAT "fDNRb"\nEND\n')
 argv=expected_command(exe,meta,out,fs,n,s,g,b);t=time.monotonic();r=subprocess.run(argv,text=True,capture_output=True,env={**os.environ,'LPSDCFN':str(cf)});elapsed=time.monotonic()-t
 out.with_suffix('.stdout.txt').write_text(r.stdout);out.with_suffix('.stderr.txt').write_text(r.stderr)
 if r.returncode:raise RuntimeError(f'C failed {name} batch {b}: {r.stderr or r.stdout}')
 rows,detail=validate_frequency_rows(out.read_text().splitlines(),freq[b*g['bins_per_batch']:(b+1)*g['bins_per_batch']])
 return rows,{'artifact':out.name,'sha256':digest(out),'frequency_validation':detail,'stdout':out.with_suffix('.stdout.txt').name,'stderr':out.with_suffix('.stderr.txt').name,'config':cf.name,'command':argv,'elapsed_seconds':elapsed,'max_child_rss_kb':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'revalidated_saved':False}
def main():
 ap=argparse.ArgumentParser(description=LABEL);ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args();run=args.run_directory.resolve()
 if not run.is_dir() or (run/'coarse_diagnostic_report.json').exists():raise RuntimeError('run must be an incomplete coarse diagnostic directory')
 cfg=yaml.safe_load((run/'coarse_config.yaml').read_text());g=cfg['coarse_grid'];s=cfg['lpsd'];n=308_717_000;fs=20_000.0;exe=Path(cfg['lpsd_executable'])
 if g['jdes']!=1000 or g['bins_per_batch']!=100 or not exe.is_file():raise RuntimeError('unexpected coarse configuration')
 if digest(Path(cfg['frozen_dense_grid']))!='7a82c9b78a22b0672189d811e25c2bb5605d9c45ac4ef550f9fad29097cccc6e':raise RuntimeError('frozen dense manifest changed')
 freq=expected_frequency(g);z=np.load(run/'coarse_c_native_grid.npz',allow_pickle=False)
 if not np.array_equal(freq,z['frequency_hz']) or np.any(z['required_dft_samples']>n):raise RuntimeError('coarse grid manifest identity/preflight failure')
 inputs={k:{'path':str(Path(v['hdf5'])),'dataset':v['dataset']} for k,v in cfg['inputs'].items()}
 for k,v in cfg['inputs'].items():
  if digest(Path(v['hdf5']))!=v['expected_file_sha256']:raise RuntimeError(f'{k} calibrated input hash mismatch')
 records=[];values={}
 # Existing complete X batches are immutable and revalidated. x-batch-009.txt is intentionally never read or overwritten.
 xrows=[]
 for b in range(9):
  rows,record=validate_saved(run,'x',b,g,freq);xrows.append(rows);record.update(channel='x',batch=b);records.append(record)
 rows,record=run_new(run,'x',9,exe,inputs['x'],fs,n,s,g,freq);xrows.append(rows);record.update(channel='x',batch=9,reason='retry after preserved incomplete x-batch-009.txt');records.append(record);values['x']=np.vstack(xrows)
 yrows=[]
 for b in range(10):
  rows,record=run_new(run,'y',b,exe,inputs['y'],fs,n,s,g,freq);yrows.append(rows);record.update(channel='y',batch=b);records.append(record)
 values['y']=np.vstack(yrows)
 for name,a in values.items():
  if a.shape!=(1000,5) or not np.array_equal(a[:,0],freq) or not np.all(np.diff(a[:,0])>0) or not np.all(np.isfinite(a)) or not np.all(a[:,1]>=0) or not np.all(a[:,2]>0):raise RuntimeError(f'complete grid coverage/validity failure for {name}')
  np.save(run/f'{name}_coarse_complete_fDNRb.npy',a)
  for q,v,lab in [('psd',a[:,1],'PSD (ellipticity²/Hz)'),('asd',np.sqrt(a[:,1]),'ASD (ellipticity/sqrt(Hz))')]:
   fig,ax=plt.subplots(figsize=(9,5));ax.loglog(a[:,0],v,lw=.8);ax.set(xlabel='Frequency (Hz)',ylabel=lab,title=f'{LABEL}: {name.upper()} {q.upper()}');ax.grid(True,which='both',alpha=.3);fig.tight_layout();fig.savefig(run/f'{name}_{q}_unsmoothed_complete.png',dpi=150);plt.close(fig)
 report={'label':LABEL,'resume_completed':True,'preserved_incomplete_artifact':'x-batch-009.txt','frozen_dense_manifest_sha256':'7a82c9b78a22b0672189d811e25c2bb5605d9c45ac4ef550f9fad29097cccc6e','grid_coverage':{'x':[0,999],'y':[0,999]},'records':records,'complete_outputs':{k:{'numeric':f'{k}_coarse_complete_fDNRb.npy','psd_plot':f'{k}_psd_unsmoothed_complete.png','asd_plot':f'{k}_asd_unsmoothed_complete.png'} for k in values},'not_performed':['X/Y combination','FFT','candidates','limits']}
 (run/'coarse_diagnostic_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'run_directory':str(run),'records':len(records)},indent=2))
if __name__=='__main__':main()
