#!/usr/bin/env python3
"""COARSE DIAGNOSTIC ONLY — NOT FOR CANDIDATES OR LIMITS."""
from __future__ import annotations
import argparse, hashlib, json, os, resource, subprocess, sys, time, uuid
from pathlib import Path
import h5py, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, yaml
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from ape_c_lpsd.c_grid import c_frequency,c_required_dft_samples
from ape_c_lpsd.c_output_validation import validate_frequency_rows

LABEL='COARSE DIAGNOSTIC ONLY — NOT FOR CANDIDATES OR LIMITS'
def digest(path:Path)->str:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def snap(path:Path)->dict:
 return {'path':str(path),'sha256':digest(path)}
def verify_input(spec:dict,n:int,fs:float)->dict:
 p=Path(spec['hdf5']);want=spec['expected_file_sha256'];got=digest(p)
 if got!=want:raise RuntimeError(f'HDF5 SHA-256 mismatch for {p}: {got} != {want}')
 with h5py.File(p,'r') as f:
  if spec['dataset'] not in f:raise RuntimeError(f'missing dataset {spec["dataset"]} in {p}')
  ds=f[spec['dataset']]; attrs={str(k):str(v) for k,v in ds.attrs.items()}
  if ds.shape!=(n,) or ds.dtype!=np.dtype('float64'):raise RuntimeError(f'HDF5 shape/dtype mismatch for {p}')
  required={'unit':'1 (dimensionless ellipticity)','sample_rate_hz':fs,'sample_count':n,'recorded_channel':spec['recorded_channel'],'calibration_method':'reference_scalar_denominator_separate_quadrature','raw_to_reduced_mapping':'raw_index // 1000'}
  for k,v in required.items():
   if ds.attrs.get(k)!=v:raise RuntimeError(f'HDF5 attribute mismatch {p}:{k}={ds.attrs.get(k)!r}')
 return {'path':str(p),'file_sha256':got,'dataset':spec['dataset'],'shape':[n],'dtype':'float64','attributes':attrs}
def command(exe:Path,input_meta:dict,out:Path,cfg:Path,fs:float,n:int,settings:dict,grid:dict,batch:int)->list[str]:
 return [str(exe),'-d','-D',input_meta['dataset'],'-f',repr(fs),'-h','0','-i',input_meta['path'],'-o',str(out),'-b','0','-e',repr((n-1)/fs),'-s',repr(grid['fmin_hz']),'-t',repr(grid['fmax_hz']),'-J',str(grid['jdes']),'-n',str(grid['bins_per_batch']),'-N',str(batch),'-w',str(settings['window']),'-p',str(settings['psll_db']),'-r',str(settings['linear_regression']),'-l',str(settings['overlap_percent']),'-a',str(settings['desired_averages']),'-m',str(settings['minimum_averages']),'-x',str(settings['input_scale']),'-M',str(settings['max_memory_power'])]
def main():
 ap=argparse.ArgumentParser(description=LABEL);ap.add_argument('--config',type=Path,default=ROOT/'configs/coarse_diagnostic.yaml');args=ap.parse_args()
 cfg=yaml.safe_load(args.config.read_text());grid=dict(cfg['coarse_grid']);s=cfg['lpsd'];n_expected=308_717_000;fs_expected=20_000.0
 if grid['jdes']!=1000 or grid['bins_per_batch']<=0 or grid['jdes']%grid['bins_per_batch']:raise RuntimeError('coarse batching must exactly partition Jdes=1000')
 if {k:s[k] for k in ('method','window','psll_db','linear_regression','overlap_percent','desired_averages','minimum_averages','input_scale')} != {'method':0,'window':-2,'psll_db':200,'linear_regression':1,'overlap_percent':50,'desired_averages':100,'minimum_averages':1,'input_scale':1}:raise RuntimeError('coarse diagnostic settings differ from frozen direct profile')
 exe=Path(cfg['lpsd_executable']);
 if not exe.is_file() or not os.access(exe,os.X_OK):raise RuntimeError(f'missing executable: {exe}')
 segment=json.loads(Path(cfg['segment_manifest']).read_text())
 if segment.get('sample_count')!=n_expected or float(segment.get('sample_rate_hz'))!=fs_expected or len(segment.get('parts',[]))!=3:raise RuntimeError('selected segment identity mismatch')
 calibration=json.loads(Path(cfg['calibration_profile']).read_text())
 if calibration.get('method')!='reference_scalar_denominator_applied_separately_to_pem1_x_and_pem1_Y':raise RuntimeError('calibration provenance mismatch')
 dense=Path(cfg['frozen_dense_grid']); dense_hash=digest(dense)
 with np.load(dense,allow_pickle=False) as z:
  if len(z['frequency_hz'])!=4_947_016:raise RuntimeError('frozen dense manifest identity mismatch')
 run=Path(cfg['output_root'])/f'coarse-direct-c-lpsd-diagnostic-{uuid.uuid4().hex}';run.mkdir(parents=True,exist_ok=False)
 (run/'coarse_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=True));(run/'selected_segment.json').write_text(json.dumps(segment,indent=2,sort_keys=True)+'\n');(run/'calibration_provenance.json').write_text(json.dumps({'profile':calibration,'source':snap(Path(cfg['calibration_profile']))},indent=2,sort_keys=True)+'\n')
 g=np.log(grid['fmax_hz']/grid['fmin_hz']);step=np.expm1(g/(grid['jdes']-1));freq=np.asarray([c_frequency(i,grid['fmin_hz'],grid['fmax_hz'],grid['jdes']) for i in range(grid['jdes'])]);length=np.asarray([c_required_dft_samples(float(f),float(step),fs_expected) for f in freq],dtype=np.int64)
 bad=np.flatnonzero(length>n_expected)
 if bad.size:raise RuntimeError(f'coarse preflight failed node={int(bad[0])} f={freq[bad[0]]} required={length[bad[0]]} available={n_expected}')
 np.savez_compressed(run/'coarse_c_native_grid.npz',frequency_hz=freq,required_dft_samples=length)
 grid_record={'label':LABEL,'formula':'lpsd.c calc_params','jdes':1000,'fmin_hz':65.0,'fmax_hz':9000.0,'relative_step':float(step),'first_frequency_hz':float(freq[0]),'last_frequency_hz':float(freq[-1]),'maximum_required_samples':int(length.max()),'maximum_required_index':int(length.argmax()),'nodes_exceeding_segment':0,'npz_sha256':digest(run/'coarse_c_native_grid.npz'),'frozen_dense_grid_sha256':dense_hash}
 (run/'coarse_grid_manifest.json').write_text(json.dumps(grid_record,indent=2,sort_keys=True)+'\n')
 inputs={k:verify_input(v,n_expected,fs_expected) for k,v in cfg['inputs'].items()}
 records=[]
 for name,meta in inputs.items():
  rows=[]
  for batch in range(grid['jdes']//grid['bins_per_batch']):
   out=run/f'{name}-batch-{batch:03d}.txt';cf=out.with_suffix('.cfg');cf.write_text('GNUTERM "APE coarse diagnostic"\nFORMAT "fDNRb"\nEND\n')
   argv=command(exe,meta,out,cf,fs_expected,n_expected,s,grid,batch);started=time.monotonic();r=subprocess.run(argv,text=True,capture_output=True,env={**os.environ,'LPSDCFN':str(cf)});elapsed=time.monotonic()-started
   out.with_suffix('.stdout.txt').write_text(r.stdout);out.with_suffix('.stderr.txt').write_text(r.stderr)
   if r.returncode:raise RuntimeError(f'lpsd-exec failed {name} batch {batch}: {r.stderr or r.stdout}')
   expected=freq[batch*grid['bins_per_batch']:(batch+1)*grid['bins_per_batch']];a,validation=validate_frequency_rows(out.read_text().splitlines(),expected);rows.append(a)
   records.append({'channel':name,'batch':batch,'global_indices':[batch*grid['bins_per_batch'],(batch+1)*grid['bins_per_batch']-1],'command':argv,'elapsed_seconds':elapsed,'max_child_rss_kb':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'output':out.name,'output_sha256':digest(out),'frequency_validation':validation})
  allrows=np.vstack(rows)
  if allrows.shape!=(1000,5) or not np.all(np.diff(allrows[:,0])>0):raise RuntimeError(f'coarse merge validation failed for {name}')
  np.save(run/f'{name}_coarse_fDNRb.npy',allrows)
  for quantity,values,ylabel in [('psd',allrows[:,1],'PSD (ellipticity²/Hz)'),('asd',np.sqrt(allrows[:,1]),'ASD (ellipticity/sqrt(Hz))')]:
   fig,ax=plt.subplots(figsize=(9,5));positive=values>0;ax.loglog(allrows[positive,0],values[positive],lw=0.8);ax.set(xlabel='Frequency (Hz)',ylabel=ylabel,title=f'{LABEL}: {name} {quantity.upper()}');ax.grid(True,which='both',alpha=.3);fig.tight_layout();fig.savefig(run/f'{name}_{quantity}_unsmoothed.png',dpi=150);plt.close(fig)
 report={'label':LABEL,'configuration':snap(args.config),'executable':snap(exe),'selected_segment':snap(Path(cfg['segment_manifest'])),'calibration_profile':snap(Path(cfg['calibration_profile'])),'input_verification':inputs,'grid':grid_record,'batches':records,'outputs':{name:{'numeric':f'{name}_coarse_fDNRb.npy','psd_plot':f'{name}_psd_unsmoothed.png','asd_plot':f'{name}_asd_unsmoothed.png'} for name in inputs},'not_performed':['dense production calculation','X/Y combination','FFT','candidate search','limits']}
 (run/'coarse_diagnostic_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'run_directory':str(run),'report':report},indent=2))
if __name__=='__main__':main()
