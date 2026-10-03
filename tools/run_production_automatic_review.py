#!/usr/bin/env python3
"""Review-only chronology-gated automatic acquisition inspection.
Never writes HDF5 or calls an estimator. Chronology is accepted only from the
researcher-provided manifest written in this run directory.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time, uuid
from pathlib import Path
import numpy as np
import yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
TAG='REVIEW-ONLY WHOLE-ACQUISITION AUTOMATIC INGESTION'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def env(v, maximum):
 v=np.asarray(v,float); w=max(1,(len(v)+maximum-1)//maximum); q=np.pad(v,(0,(-len(v))%w),constant_values=np.nan).reshape(-1,w)
 with np.errstate(all='ignore'): return np.arange(len(q))*w,np.nanmin(q,1),np.nanmax(q,1),w

def masked(v,m,maximum):
 i,lo,hi,w=env(v,maximum); a=np.asarray(v,float); m=np.asarray(m,bool)
 def b(q):
  q=np.pad(q,(0,(-len(q))%w),constant_values=np.nan).reshape(-1,w); good=np.any(np.isfinite(q),1); x=np.full(len(q),np.nan); y=np.full(len(q),np.nan)
  if good.any(): x[good]=np.nanmin(q[good],1); y[good]=np.nanmax(q[good],1)
  return x,y
 vl,vh=b(np.where(m,a,np.nan)); return i,lo,hi,vl,vh,w

def plot_raw(i,lo,hi,fs,p,title,unit,w):
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(figsize=(11,4)); ax.fill_between(i/fs,lo,hi,color='.75',alpha=.7,label='raw data: light-gray original-time envelope');ax.set(title=title,xlabel='Original acquisition time [s]',ylabel=unit);ax.legend(fontsize=7);ax.grid(alpha=.2);ax.text(.01,.01,f'{TAG}; {w:,} original samples/bin; no resampling',transform=ax.transAxes,fontsize=7);fig.tight_layout();fig.savefig(p,dpi=160);plt.close(fig)
 return {'path':str(p),'sha256':sha(p),'display':'raw original-time envelope','bin_samples':int(w)}
def plot_clean(i,lo,hi,fs,p,title,unit,w,valid,total):
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(figsize=(11,4));ax.fill_between(i/fs,lo,hi,color='tab:blue',alpha=.8,label='retained samples only: original-time envelope');ax.set(title=title,xlabel='Original acquisition time [s]',ylabel=unit);ax.legend(fontsize=7);ax.grid(alpha=.2);ax.text(.01,.01,f'valid fraction={valid/total:.9f}; rejected={total-valid:,}; gaps are not joined; {w:,} original samples/bin',transform=ax.transAxes,fontsize=7);fig.tight_layout();fig.savefig(p,dpi=160);plt.close(fig)
 return {'path':str(p),'sha256':sha(p),'display':'retained-only original-time envelope','bin_samples':int(w),'valid_fraction':valid/total,'rejected':int(total-valid)}
def concat(parts,k): return np.concatenate([p[k] for p in parts]) if parts else np.empty(0)
def true_runs(m,offset):
 e=np.flatnonzero(np.r_[False,m,False][1:]!=np.r_[False,m,False][:-1]);return [(offset+int(a),offset+int(b)) for a,b in zip(e[::2],e[1::2])]
def main():
 ap=argparse.ArgumentParser(description=TAG); ap.add_argument('--config',type=Path,default=ROOT/'configs/production_automatic_review.yaml'); args=ap.parse_args(); c=yaml.safe_load(args.config.read_text())
 root=Path(c['input_root']); prefix=c['filename_prefix']; first,last=int(c['first_file_index']),int(c['last_file_index']); chunk=int(c['chunk_samples'])
 if chunk<=0 or chunk%1000: raise RuntimeError('chunk_samples must be positive and divisible by 1000')
 paths=[root/f'{prefix}{n:04d}.dxd' for n in range(first,last+1)]
 missing=[str(p) for p in paths if not p.is_file()]
 if missing: raise RuntimeError('missing researcher-manifest source(s): '+', '.join(missing))
 ref=Path(c['reference_project']);sys.path.insert(0,str(ref/'src'))
 from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths
 from ape_c_lpsd.verified_segment import Part,_quality,_read,_check_metadata,_mask
 from ape_c_lpsd.calibrated_review import calibrated_quadratures, calibration_profile
 from ape_c_lpsd.ingestion.models import ValidRun
 from ape_c_lpsd.ingestion.segments import build_global_segments,choose_longest_global_segment
 from ape_analysis.calibration import PolarizationCalibrationParameters
 dp=resolve_reference_dewesoft_paths(ref); params=PolarizationCalibrationParameters(); run=Path(c['output_root'])/f'production-automatic-review-{uuid.uuid4().hex}';run.mkdir(parents=True)
 manifest={'kind':'ordered_source_manifest','label':TAG,'chronology_basis':'researcher_provided_or_explicit_manifest','chronology_statement':'Researcher provided numeric order 0000.dxd through 0013.dxd and continuity across adjacent parts. Dewesoft timestamps and filesystem times are not used to establish chronology.','sources':[{'ordinal':o,'relative_path':p.name,'chronology_evidence':'researcher-provided numerical acquisition order and continuity','boundary_to_next_verified':o<len(paths)-1,'file_size_bytes':p.stat().st_size,'sha256':None,'input_hash_note':'not computed to avoid a second full raw-data pass'} for o,p in enumerate(paths)]}
 (run/'researcher_provided_source_manifest.yaml').write_text(yaml.safe_dump(manifest,sort_keys=False));(run/'config_snapshot.json').write_text(json.dumps(c,indent=2,sort_keys=True)+'\n')
 allenv={n:[] for n in ('cavity_trans','pem2_X','pem2_Y','pem1_x','pem1_Y')}; qualities=[];infos=[];runs=[];counts={};global_offset=0;started=time.monotonic()
 for ordinal,path in enumerate(paths):
  # Probe metadata exactly; the one sample is part of the bounded raw inspection.
  probe=_read(path,0,1,'pem1_x',dewesoft_paths=dp); total=int(probe.metadata['total_raw_sample_count']); fs=float(probe.sample_rate)
  if total%1000: raise RuntimeError(f'{path.name}: raw count {total} violates frozen 1000:1 mapping')
  if fs!=20000.: raise RuntimeError(f'{path.name}: sample-rate {fs} differs from frozen reference 20000 Hz')
  quality=_quality(path,total//1000,dewesoft_paths=dp); qualities.append(quality); counts[ordinal]=total
  rm=np.asarray(quality.combined_good,bool); red={'cavity_trans':np.asarray([r.cavity_transmission_average_v for r in quality.records]),'pem2_X':np.asarray([r.pem2_x_average_v for r in quality.records]),'pem2_Y':np.asarray([r.pem2_y_average_v for r in quality.records])}
  for name,v in red.items():
   i,lo,hi,vl,vh,w=masked(v,rm,int(c['display_max_points_per_reduced_file']));allenv[name].append({'i':global_offset+i*1000,'lo':lo,'hi':hi,'vl':vl,'vh':vh,'w':w*1000})
  file_valid=0; file_runs=[]; active_run=None
  for start in range(0,total,chunk):
   count=min(chunk,total-start);x=_read(path,start,count,'pem1_x',dewesoft_paths=dp);y=_read(path,start,count,'pem1_Y',dewesoft_paths=dp);_check_metadata(x,y,Part(path,0,total)); xv=np.asarray(x.samples,float);yv=np.asarray(y.samples,float);m=_mask(xv,yv,start,quality);file_valid+=int(m.sum())
   for name,v in (('pem1_x',xv),('pem1_Y',yv)):
    i,lo,hi,vl,vh,w=masked(v,m,int(c['display_envelopes_per_raw_chunk']));allenv[name].append({'i':global_offset+start+i,'lo':lo,'hi':hi,'vl':vl,'vh':vh,'w':w})
   for a,b in true_runs(m,start):
    if active_run is not None and a==active_run[1]: active_run=(active_run[0],b)
    else:
     if active_run is not None: file_runs.append(active_run)
     active_run=(a,b)
   if not m[-1] and active_run is not None:
    file_runs.append(active_run); active_run=None
   print(f'{TAG}: {path.name} chunk [{start},{start+count}) valid={int(m.sum())}',flush=True)
  if active_run is not None: file_runs.append(active_run)
  for a,b in file_runs:runs.append(ValidRun(ordinal,a,b,fs,ordinal>0,path))
  infos.append({'ordinal':ordinal,'path':str(path),'raw_sample_count':total,'sample_rate_hz':fs,'reduced_record_count':len(quality.records),'quality_good_raw_equivalent':int(rm.sum())*1000,'raw_calibration_valid_and_good':file_valid,'raw_valid_fraction':file_valid/total,'file_runs':len(file_runs)})
  global_offset+=total
 segments=build_global_segments(tuple(runs),source_sample_counts=counts); selected=choose_longest_global_segment(segments)
 # Whole-acquisition raw/cleaned plots. Plot source-boundary time positions without joining across sources.
 plots={}
 for name,items in allenv.items():
  i=concat(items,'i');lo=concat(items,'lo');hi=concat(items,'hi');vl=concat(items,'vl');vh=concat(items,'vh');w=max(x['w'] for x in items);valid=sum(x['raw_calibration_valid_and_good'] for x in infos) if name.startswith('pem1') else sum(x['quality_good_raw_equivalent'] for x in infos);total=sum(counts.values())
  d=run/'raw'/name;d.mkdir(parents=True);plots[f'raw/{name}']=plot_raw(i,lo,hi,20000,d/'raw_time_series.png',f'{TAG}: raw {name}','V',w)
  d=run/'cleaned'/name;d.mkdir(parents=True);plots[f'cleaned/{name}']=plot_clean(i,vl,vh,20000,d/'cleaned_time_series.png',f'{TAG}: cleaned {name}','V',w,valid,total)
 # Re-read only selected valid pieces to create selected and calibrated review envelopes.
 selenv={n:[] for n in ('cavity_trans','pem2_X','pem2_Y','pem1_x','pem1_Y','X_ellipticity','Y_ellipticity','calibration_denominator_D')};sel_global=0
 selected_parts=[]
 for r in selected.parts:
  q=qualities[r.source_ordinal]; rr={'cavity_trans':np.asarray([z.cavity_transmission_average_v for z in q.records]),'pem2_X':np.asarray([z.pem2_x_average_v for z in q.records]),'pem2_Y':np.asarray([z.pem2_y_average_v for z in q.records])}; lo_block=r.raw_start//1000; hi_block=(r.raw_stop+999)//1000
  for dn,dv in rr.items():
   di,dl,dh,dw=env(dv[lo_block:hi_block],int(c['display_max_points_per_reduced_file']));selenv[dn].append({'i':sel_global+di*1000,'lo':dl,'hi':dh,'w':dw*1000})
  selected_parts.append({'ordinal':r.source_ordinal,'path':str(r.source_path),'raw_start':r.raw_start,'raw_stop':r.raw_stop,'sample_count':r.sample_count})
  for start in range(r.raw_start,r.raw_stop,chunk):
   count=min(chunk,r.raw_stop-start);x=_read(r.source_path,start,count,'pem1_x',dewesoft_paths=dp);y=_read(r.source_path,start,count,'pem1_Y',dewesoft_paths=dp);xv=np.asarray(x.samples,float);yv=np.asarray(y.samples,float);m=_mask(xv,yv,start,q)
   if not np.all(m): raise RuntimeError(f'selected run became invalid at {r.source_path.name}:{start+int(np.flatnonzero(~m)[0])}')
   den,xe,ye=calibrated_quadratures(xv,yv,q,start,params,m)
   for name,v in (('pem1_x',xv),('pem1_Y',yv),('X_ellipticity',xe),('Y_ellipticity',ye),('calibration_denominator_D',den)):
    i,lo,hi,w=env(v,int(c['display_envelopes_per_raw_chunk']));selenv[name].append({'i':sel_global+i,'lo':lo,'hi':hi,'w':w})
   sel_global+=count
 for name,items in selenv.items():
  i=concat(items,'i');lo=concat(items,'lo');hi=concat(items,'hi');w=max(x['w'] for x in items);unit='V' if name in ('cavity_trans','pem2_X','pem2_Y','pem1_x','pem1_Y','calibration_denominator_D') else '1 (dimensionless ellipticity)'; d=run/('calibrated_review' if name.startswith(('X_','Y_','calibration')) else 'selected_longest_segment')/name;d.mkdir(parents=True);plots[f'selected/{name}']=plot_clean(i,lo,hi,20000,d/'time_series.png',f'{TAG}: selected {name}',unit,w,sel_global,sel_global)
 frozen={'kind':'frozen_validated_segment','selection_mode':'automatic_researcher_chronology_review','sample_rate_hz':20000.,'unit':'V','sample_count':selected.sample_count,'parts':selected_parts,'continuity':'Joined only across endpoint-adjacent valid runs at researcher-provided verified boundaries; no Dewesoft or filesystem time ordering claim.'}
 (run/'frozen_selected_segment.json').write_text(json.dumps(frozen,indent=2,sort_keys=True)+'\n');(run/'calibration_profile.json').write_text(json.dumps(calibration_profile(params),indent=2,sort_keys=True)+'\n')
 manual=[{'path':str(root/f'{prefix}{n:04d}.dxd'),'raw_start':154903000 if n==11 else 0,'raw_stop':216000000 if n in (11,12) else 31620000} for n in (11,12,13)]
 comparison={'previous_verified_manual_segment':{'parts':manual,'sample_count':308717000},'automatic_selected_segment':frozen,'same_parts_and_boundaries':[{k:v for k,v in x.items() if k!='ordinal'} for x in selected_parts]==manual}
 (run/'manual_segment_comparison.json').write_text(json.dumps(comparison,indent=2,sort_keys=True)+'\n')
 report={'label':TAG,'elapsed_seconds':time.monotonic()-started,'source_count':len(paths),'total_raw_samples':sum(counts.values()),'per_file':infos,'valid_run_count':len(runs),'global_segment_count':len(segments),'selected_segment':frozen,'comparison_to_verified_manual_segment':comparison,'prohibited_actions_not_performed':['HDF5 export','FFT','C-LPSD','X/Y combination','candidate search','limits'],'calibration_profile':calibration_profile(params),'plot_count':len(plots)}
 (run/'plot_metadata.json').write_text(json.dumps({'label':TAG,'plots':plots,'selected_segment':frozen,'calibration_profile':calibration_profile(params)},indent=2,sort_keys=True)+'\n');(run/'review_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'run_directory':str(run),'selected':frozen,'elapsed_seconds':report['elapsed_seconds']},indent=2))
if __name__=='__main__':main()
