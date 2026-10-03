#!/usr/bin/env python3
"""Approved non-overwriting calibrated real-data C-LPSD pilot only."""
from __future__ import annotations
import argparse, hashlib, json, os, resource, subprocess, sys, time, uuid
from pathlib import Path
import h5py, numpy as np, yaml
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--segment',type=Path,required=True);ap.add_argument('--grid',type=Path,required=True);ap.add_argument('--config',type=Path,default=ROOT/'configs/run.yaml');ap.add_argument('--output-root',type=Path,default=ROOT/'results');args=ap.parse_args()
 cfg=yaml.safe_load(args.config.read_text()); seg=json.loads(args.segment.read_text()); gridj=json.loads(args.grid.with_suffix('.json').read_text()); gridhash=digest(args.grid)
 if gridhash!=gridj['manifest_npz_sha256']:raise RuntimeError(f'grid manifest hash mismatch: {gridhash} != {gridj["manifest_npz_sha256"]}')
 z=np.load(args.grid,allow_pickle=False); freq=np.asarray(z['frequency_hz'],float)
 if len(freq)!=4947016 or not np.all(np.diff(freq)>0) or freq[0]!=65. or not np.isclose(freq[-1],9000.):raise RuntimeError('frozen grid validation failed')
 if seg['sample_count']!=308717000 or seg['sample_rate_hz']!=20000. or len(seg['parts'])!=3:raise RuntimeError('frozen selected segment validation failed')
 ref=Path('/Users/user/Desktop/APE Research/APE-data-analysis');sys.path.insert(0,str(ref/'src'))
 from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths
 from ape_c_lpsd.verified_segment import Part,_quality,_read,_check_metadata,_mask
 from ape_c_lpsd.calibrated_review import calibrated_quadratures,calibration_profile
 from ape_analysis.calibration import PolarizationCalibrationParameters
 exe=Path(cfg['lpsd_executable']);
 if not exe.is_file() or not os.access(exe,os.X_OK):raise RuntimeError(f'missing executable: {exe}')
 c=gridj['c_command_options']; required={'-J':4947016,'-h':0,'-w':-2,'-p':200,'-r':1,'-l':50,'-a':100,'-m':1,'-x':1,'-s':65.0,'-t':9000.0}
 for k,v in required.items():
  if c.get(k)!=v:raise RuntimeError(f'frozen grid setting mismatch {k}: {c.get(k)!r}')
 dp=resolve_reference_dewesoft_paths(ref);params=PolarizationCalibrationParameters();run=args.output_root/f'calibrated-c-lpsd-pilot-{uuid.uuid4().hex}';run.mkdir(parents=True)
 (run/'config_snapshot.yaml').write_text(yaml.safe_dump(cfg,sort_keys=True));(run/'segment_manifest.json').write_text(json.dumps(seg,indent=2,sort_keys=True)+'\n');(run/'grid_manifest.json').write_text(json.dumps(gridj,indent=2,sort_keys=True)+'\n')
 n=seg['sample_count'];fs=seg['sample_rate_hz'];chunk=1_000_000; exports={}
 for label,chan,dset in [('x_ellipticity','pem1_x','/x_ellipticity'),('y_ellipticity','pem1_Y','/y_ellipticity')]:
  hp=run/f'{label}.h5'; h=hashlib.sha256();count=0;started=time.monotonic();valid=0
  with h5py.File(hp,'w') as f:
   ds=f.create_dataset(dset,shape=(n,),dtype='f8',chunks=(chunk,));ds.attrs.update(unit='1 (dimensionless ellipticity)',sample_rate_hz=fs,sample_count=n,recorded_channel=chan,calibration_method='reference_scalar_denominator_separate_quadrature',raw_to_reduced_mapping='raw_index // 1000')
   for p in seg['parts']:
    path=Path(p['path']);start0=int(p['raw_start']);stop=int(p['raw_stop']);q=_quality(path,(stop+999)//1000,dewesoft_paths=dp)
    for start in range(start0,stop,chunk):
     ct=min(chunk,stop-start);x=_read(path,start,ct,'pem1_x',dewesoft_paths=dp);y=_read(path,start,ct,'pem1_Y',dewesoft_paths=dp);_check_metadata(x,y,Part(path,start0,stop));xv=np.asarray(x.samples,float);yv=np.asarray(y.samples,float);m=_mask(xv,yv,start,q)
     if not np.all(m):raise RuntimeError(f'invalid retained sample {path.name}:{start+int(np.flatnonzero(~m)[0])}')
     den,xe,ye=calibrated_quadratures(xv,yv,q,start,params,m);out=xe if label.startswith('x_') else ye
     if not np.all(np.isfinite(out)) or not np.all(np.isfinite(den)):raise RuntimeError(f'nonfinite calibrated value {path.name}:{start}')
     ds[count:count+ct]=out;h.update(np.ascontiguousarray(out,dtype='<f8').tobytes());count+=ct;valid+=int(m.sum())
  stream_hash=h.hexdigest()
  with h5py.File(hp,'r') as f:
   ds=f[dset];verify=hashlib.sha256()
   if ds.shape!=(n,) or ds.dtype!=np.dtype('float64') or ds.attrs['unit']!='1 (dimensionless ellipticity)':raise RuntimeError('HDF5 metadata/type verification failed')
   for a in range(0,n,chunk):
    v=np.asarray(ds[a:min(a+chunk,n)]); 
    if not np.all(np.isfinite(v)):raise RuntimeError('HDF5 finite verification failed')
    verify.update(np.ascontiguousarray(v,dtype='<f8').tobytes())
  if verify.hexdigest()!=stream_hash:raise RuntimeError('HDF5 order/value hash verification failed')
  exports[label]={'path':str(hp),'dataset':dset,'count':count,'dtype':'float64','unit':'1 (dimensionless ellipticity)','stream_value_sha256':stream_hash,'file_sha256':digest(hp),'valid_samples':valid,'elapsed_seconds':time.monotonic()-started,'size_bytes':hp.stat().st_size}
 batches=[0,4947016//16,4947016//8-1];batch_n=8;logs=[]
 for label in exports:
  for b in batches:
   out=run/f'{label}-batch-{b}.txt'; cf=out.with_suffix('.cfg');cf.write_text('GNUTERM "APE C-LPSD"\nFORMAT "fDNRb"\nEND\n');argv=[str(exe),'-d','-D',exports[label]['dataset'],'-f','20000','-h','0','-i',exports[label]['path'],'-o',str(out),'-b','0','-e',repr((n-1)/fs),'-s','65','-t','9000','-J','4947016','-n',str(batch_n),'-N',str(b),'-w','-2','-p','200','-r','1','-l','50','-a','100','-m','1','-x','1','-M',str(c['-M'])]
   env={**os.environ,'LPSDCFN':str(cf)};t=time.monotonic();r=subprocess.run(argv,text=True,capture_output=True,env=env);(out.with_suffix('.stdout.txt')).write_text(r.stdout);(out.with_suffix('.stderr.txt')).write_text(r.stderr)
   if r.returncode:raise RuntimeError(f'lpsd-exec failed {label} batch {b}: {r.stderr}')
   rows=[]
   for line in out.read_text().splitlines():
    if line.strip() and not line.startswith('#'):
     v=line.split()
     if len(v)!=5:raise RuntimeError('malformed C output')
     rows.append([float(x) for x in v])
   a=np.asarray(rows,float);expected=freq[b*batch_n:(b+1)*batch_n]
   if a.shape!=(batch_n,5) or not np.all(np.isfinite(a)) or not np.all(a[:,1]>=0) or not np.all(a[:,2]>0) or not np.all(np.diff(a[:,0])>0) or not np.allclose(a[:,0],expected,rtol=2e-12,atol=2e-10):raise RuntimeError(f'strict C output validation failed {label} batch {b}')
   np.save(run/f'{label}-batch-{b}.npy',a);logs.append({'channel':label,'batch':b,'node_indices':[b*batch_n,(b+1)*batch_n-1],'command':argv,'elapsed_seconds':time.monotonic()-t,'max_rss_kb':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'output':str(out),'output_sha256':digest(out),'returned_frequency_hz':[float(a[0,0]),float(a[-1,0])],'psd_unit':'ellipticity^2/Hz','asd_unit':'ellipticity/sqrt(Hz)'})
 estimate=sum(x['elapsed_seconds'] for x in logs)/len(logs)*4947016/batch_n*2
 report={'kind':'approved_calibrated_real_data_c_lpsd_pilot','segment_samples':n,'sample_rate_hz':fs,'grid_npz_sha256':gridhash,'executable_sha256':digest(exe),'exports':exports,'pilot_batches':logs,'calibration_profile':calibration_profile(params),'full_production_runtime_estimate_seconds':estimate,'estimate_note':'linear extrapolation from six eight-node pilot invocations; substantial low-frequency/runtime uncertainty remains','prohibited_not_performed':['full production grid','X/Y combination','FFT','candidates','limits']}
 (run/'pilot_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'run_directory':str(run),'report':report},indent=2))
if __name__=='__main__':main()
