#!/usr/bin/env python3
"""TEMPORARY SINGLE-FILE SANITY CHECK: review-only; no HDF5/FFT/C-LPSD."""
from __future__ import annotations
import argparse, hashlib, json, os, sys, time, uuid
from pathlib import Path
import numpy as np, yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
TAG='TEMPORARY SINGLE-FILE SANITY CHECK'
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(1<<20),b''):h.update(c)
 return h.hexdigest()
def envelope(v, maximum=10000):
 v=np.asarray(v,float); w=max(1,int(np.ceil(len(v)/maximum))); q=np.pad(v,(0,(-len(v))%w),constant_values=np.nan).reshape(-1,w); return np.arange(len(q))*w,np.nanmin(q,axis=1),np.nanmax(q,axis=1),w
def plot_env(v,fs,path,title,unit):
 import matplotlib.pyplot as plt
 i,lo,hi,w=envelope(v); fig,ax=plt.subplots(figsize=(10,4));ax.fill_between(i/fs,lo,hi,alpha=.45,label=f'{TAG}; display min/max envelope; {w} original samples/bin; no resampling')
 ax.set(title=title,xlabel='Time [s]',ylabel=unit);ax.legend(fontsize=7);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)
 return {'path':str(path),'sha256':digest(path),'display_bin_width_samples':w}
def masked_envelope(v, mask, maximum=10000):
    """Envelope preserving original bin positions; NaNs leave visible cleaning gaps."""
    v=np.asarray(v,float); mask=np.asarray(mask,bool)
    i,lo,hi,w=envelope(v,maximum)
    valid=np.where(mask,v,np.nan); rejected=np.where(~mask,v,np.nan)
    def bounds(a):
        q=np.pad(a,(0,(-len(a))%w),constant_values=np.nan).reshape(-1,w)
        with np.errstate(all='ignore'): return np.nanmin(q,axis=1),np.nanmax(q,axis=1)
    vlo,vhi=bounds(valid); rlo,rhi=bounds(rejected)
    return i,lo,hi,vlo,vhi,rlo,rhi,w

def plot_raw_cleaned_envelope(i,lo,hi,vlo,vhi,rlo,rhi,fs,path,title,unit,valid_count,total,display_width):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(11,4.25))
    ax.fill_between(i/fs,lo,hi,color='0.75',alpha=.65,label='raw data (light-gray display envelope)')
    ax.fill_between(i/fs,vlo,vhi,color='tab:blue',alpha=.75,label='retained valid samples (original-time envelope)')
    if total-valid_count:
        ax.fill_between(i/fs,rlo,rhi,color='red',alpha=.85,label='rejected samples (original-time envelope)')
    else:
        ax.text(.5,.94,'100% valid: cleaned trace equals raw trace.',transform=ax.transAxes,ha='center',va='top',color='tab:blue',weight='bold')
    ax.set(title=title,xlabel='Original time [s]',ylabel=unit)
    ax.legend(fontsize=7,loc='best'); ax.grid(alpha=.2)
    ax.text(.01,.01,f'valid fraction = {valid_count/total:.9f}; rejected = {total-valid_count:,}; display envelope = {display_width:,} original samples/bin; no resampling',transform=ax.transAxes,fontsize=7,va='bottom')
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)
    return {'path':str(path),'sha256':digest(path),'display_method':'original-time min/max envelopes','display_bin_width_samples':int(display_width),'valid_fraction':valid_count/total,'rejected_sample_count':total-valid_count}


def plot_raw_envelope(i,lo,hi,fs,path,title,unit,display_width):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(11,4.25))
    ax.fill_between(i/fs,lo,hi,color='0.75',alpha=.70,label='raw data (light-gray display envelope)')
    ax.set(title=title,xlabel='Original time [s]',ylabel=unit); ax.legend(fontsize=7); ax.grid(alpha=.2)
    ax.text(.01,.01,f'display envelope = {display_width:,} original samples/bin; no resampling',transform=ax.transAxes,fontsize=7,va='bottom')
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)
    return {'path':str(path),'sha256':digest(path),'display_method':'original-time min/max envelopes','display_bin_width_samples':int(display_width)}

def plot_cleaned_envelope(i,vlo,vhi,fs,path,title,unit,valid_count,total,display_width):
    """Plot only retained values; NaN bins retain original-time gaps."""
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(11,4.25))
    ax.fill_between(i/fs,vlo,vhi,color='tab:blue',alpha=.80,label='retained valid samples (original-time envelope)')
    ax.set(title=title,xlabel='Original time [s]',ylabel=unit); ax.legend(fontsize=7); ax.grid(alpha=.2)
    ax.text(.01,.01,f'valid fraction = {valid_count/total:.9f}; rejected = {total-valid_count:,}; display envelope = {display_width:,} original samples/bin; no resampling; gaps are not joined',transform=ax.transAxes,fontsize=7,va='bottom')
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)
    return {'path':str(path),'sha256':digest(path),'display_method':'retained-only original-time min/max envelopes','display_bin_width_samples':int(display_width),'valid_fraction':valid_count/total,'rejected_sample_count':total-valid_count}

def plot_no_invalid_boundary(path,title,unit):
    """Required boundary-zoom artifact when no invalid/valid boundary exists."""
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,3.5));ax.axis('off')
    ax.text(.5,.58,'100% valid: cleaned trace equals raw trace.',ha='center',va='center',fontsize=15,weight='bold',color='tab:blue')
    ax.text(.5,.35,'No invalid/valid boundary exists in this single-file review; therefore no boundary zoom can be drawn.',ha='center',va='center',fontsize=10)
    ax.set_title(title);fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)
    return {'path':str(path),'sha256':digest(path),'kind':'no-boundary notice','unit':unit}
def main():
 ap=argparse.ArgumentParser(description=TAG);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output-root',type=Path,default=ROOT/'results');ap.add_argument('--config',type=Path,default=ROOT/'configs/ingestion.yaml');ap.add_argument('--chunk-samples',type=int,default=1_000_000);args=ap.parse_args()
 if not args.input.is_file():raise RuntimeError('missing single-file source')
 if args.chunk_samples<=0 or args.chunk_samples%1000:raise RuntimeError('chunk samples must be positive and 1000-aligned')
 cfg=yaml.safe_load(args.config.read_text())
 reference_root=Path(cfg.get('reference_project',''))
 if not reference_root: raise RuntimeError(f'missing reference_project in {args.config}')
 sys.path.insert(0,str(reference_root/'src'))
 from ape_c_lpsd.ingestion.reference_paths import resolve_reference_dewesoft_paths
 from ape_c_lpsd.verified_segment import Part,_quality,_read,_check_metadata,_mask
 from ape_c_lpsd.calibrated_review import calibrated_quadratures, calibration_profile
 dewesoft_paths=resolve_reference_dewesoft_paths(reference_root)
 run=args.output_root/f'temporary-single-file-sanity-{uuid.uuid4().hex}';run.mkdir(parents=True)
 (run/'TEMPORARY_SINGLE_FILE_source_manifest.yaml').write_text(yaml.safe_dump({'chronology_basis':'researcher_provided_or_explicit_manifest','label':TAG,'sources':[{'relative_path':args.input.name,'raw_sample_count':216000000,'chronology_evidence':'single-file sanity check only; no cross-file continuity asserted','boundary_to_next_verified':False,'sha256':digest(args.input)}]},sort_keys=False))
 (run/'config_snapshot.json').write_text(json.dumps({'label':TAG,'ingestion_config':cfg,'reference_project_root':str(dewesoft_paths.reference_project_root),'dewesoft_wrapper_path':str(dewesoft_paths.wrapper_path),'dewesoft_library_path':str(dewesoft_paths.library_path),'single_source':str(args.input),'source_sha256':digest(args.input),'no_cross_file_continuity':True,'no_hdf5':True,'no_fft':True,'no_c_lpsd':True},indent=2,sort_keys=True)+'\n')
 started=time.monotonic(); n=216_000_000; part=Part(args.input,0,n);quality=_quality(args.input,n//1000,dewesoft_paths=dewesoft_paths)
 from ape_analysis.calibration import PolarizationCalibrationParameters
 calibration_parameters=PolarizationCalibrationParameters()
 red={'cavity_trans':np.array([r.cavity_transmission_average_v for r in quality.records]),'pem2_X':np.array([r.pem2_x_average_v for r in quality.records]),'pem2_Y':np.array([r.pem2_y_average_v for r in quality.records])}; reduced_mask=np.array(quality.combined_good,dtype=bool)
 # streaming validity and display summaries for PEM1 quadratures
 raw_env={'pem1_x':[],'pem1_Y':[]}; calibrated_env={'X_ellipticity':[],'Y_ellipticity':[],'calibration_denominator_D':[]}; valid_count=0;best=(0,0);active=None;chunks=0
 for start in range(0,n,args.chunk_samples):
  count=min(args.chunk_samples,n-start);x,y=_read(args.input,start,count,'pem1_x',dewesoft_paths=dewesoft_paths),_read(args.input,start,count,'pem1_Y',dewesoft_paths=dewesoft_paths);_check_metadata(x,y,part); xv=np.asarray(x.samples,float);yv=np.asarray(y.samples,float);mask=_mask(xv,yv,start,quality)
  denominator,xe,ye=calibrated_quadratures(xv,yv,quality,start,calibration_parameters,mask)
  valid_count+=int(mask.sum()); chunks+=1
  for name,v in (('pem1_x',xv),('pem1_Y',yv)):
   i,lo,hi,vlo,vhi,rlo,rhi,w=masked_envelope(v,mask,100);raw_env[name].append((start+i,lo,hi,vlo,vhi,rlo,rhi,w))
  for calibrated_name,calibrated_values in (('X_ellipticity',xe),('Y_ellipticity',ye),('calibration_denominator_D',np.where(mask,denominator,np.nan))):
   ci,clo,chi,cvlo,cvhi,crlo,crhi,cw=masked_envelope(calibrated_values,mask,100)
   calibrated_env[calibrated_name].append((start+ci,clo,chi,cvlo,cvhi,crlo,crhi,cw))
  padded=np.r_[False,mask,False];edges=np.flatnonzero(padded[1:]!=padded[:-1])
  for a,b in zip(edges[::2],edges[1::2]):
   a+=start;b+=start
   if active is not None and a==active[1]:active=(active[0],b)
   else:
    if active and active[1]-active[0]>best[1]-best[0]:best=active
    active=(a,b)
  if not mask[-1] and active:
   if active[1]-active[0]>best[1]-best[0]:best=active
   active=None
  print(f'{TAG}: chunk {chunks}, raw [{start},{start+count}), valid={int(mask.sum())}',flush=True)
 if active and active[1]-active[0]>best[1]-best[0]:best=active
 if best[1]<=best[0]:raise RuntimeError('cleaning refusal: no valid raw interval')
 # compact outputs and plots; raw X/Y use envelope points already accumulated
 meta={'label':TAG,'source':str(args.input),'source_sha256':digest(args.input),'raw_sample_count':n,'sample_rate_hz':20000,'valid_sample_count':valid_count,'valid_fraction':valid_count/n,'selected_start':int(best[0]),'selected_stop':int(best[1]),'selected_sample_count':int(best[1]-best[0]),'selected_duration_seconds':float((best[1]-best[0])/20000),'chunks':chunks,'elapsed_seconds':time.monotonic()-started,'no_cross_file_continuity':True,'calibration_profile':calibration_profile(calibration_parameters)}
 (run/'frozen_temporary_selected_segment.json').write_text(json.dumps({'kind':'frozen_validated_segment','label':TAG,'selection_mode':'automatic_single_file_sanity','sample_rate_hz':20000,'unit':'V','parts':[{'path':str(args.input),'raw_start':int(best[0]),'raw_stop':int(best[1])}],'sample_count':int(best[1]-best[0]),'continuity':'one source only; no cross-file assertion'},indent=2)+'\n')
 plots={}
 for name,v in red.items():
  valid=int(reduced_mask.sum()); total=len(reduced_mask)
  i,lo,hi,vlo,vhi,rlo,rhi,w=masked_envelope(v,reduced_mask,10000)
  d=run/'raw'/name;d.mkdir(parents=True);plots[f'raw/{name}']=plot_raw_envelope(i,lo,hi,20,d/'raw_time_series.png',f'{TAG}: raw {name}','V',w)
  d=run/'cleaned'/name;d.mkdir(parents=True);plots[f'cleaned/{name}']=plot_cleaned_envelope(i,vlo,vhi,20,d/'cleaned_time_series.png',f'{TAG}: cleaned {name}','V',valid,total,w)
  plots[f'cleaned-boundary/{name}']=plot_no_invalid_boundary(d/'boundary_zoom.png',f'{TAG}: full-resolution boundary zoom {name}','V') if valid==total else {'status':'boundary zoom required for invalid data'}
  d=run/'selected_longest_segment'/name;d.mkdir(parents=True);a=best[0]//1000;b=(best[1]+999)//1000;plots[f'selected/{name}']=plot_env(v[a:b],20,d/'overview.png',f'{TAG}: selected {name}','V')
 for name in ('pem1_x','pem1_Y'):
  q=raw_env[name]; i=np.concatenate([x[0] for x in q]);lo=np.concatenate([x[1] for x in q]);hi=np.concatenate([x[2] for x in q]);vlo=np.concatenate([x[3] for x in q]);vhi=np.concatenate([x[4] for x in q]);rlo=np.concatenate([x[5] for x in q]);rhi=np.concatenate([x[6] for x in q]);w=q[0][7]
  d=run/'raw'/name;d.mkdir(parents=True);plots[f'raw/{name}']={'path':str(d/'overview.png'),'sha256':None,'display_method':'100 original-sample envelopes'}
  plot_raw_envelope(i,lo,hi,20000,d/'raw_time_series.png',f'{TAG}: raw {name}','V',w);plots[f'raw/{name}']['path']=str(d/'raw_time_series.png');plots[f'raw/{name}']['sha256']=digest(d/'raw_time_series.png')
  d=run/'cleaned'/name;d.mkdir(parents=True);plots[f'cleaned/{name}']=plot_cleaned_envelope(i,vlo,vhi,20000,d/'cleaned_time_series.png',f'{TAG}: cleaned {name}','V',valid_count,n,w)
  plots[f'cleaned-boundary/{name}']=plot_no_invalid_boundary(d/'boundary_zoom.png',f'{TAG}: full-resolution boundary zoom {name}','V') if valid_count==n else {'status':'boundary zoom required for invalid data'}
  d=run/'selected_longest_segment'/name;d.mkdir(parents=True);keep=(i>=best[0])&(i<best[1]);plots[f'selected/{name}']=plot_raw_cleaned_envelope(i[keep],lo[keep],hi[keep],vlo[keep],vhi[keep],rlo[keep],rhi[keep],20000,d/'overview.png',f'{TAG}: selected {name}','V',valid_count,n,w)
 # Separate calibrated-review plots; gaps are retained as NaN envelope gaps on the original axis.
 for name in ('X_ellipticity','Y_ellipticity','calibration_denominator_D'):
  q=calibrated_env[name]; i=np.concatenate([x[0] for x in q]);lo=np.concatenate([x[1] for x in q]);hi=np.concatenate([x[2] for x in q]);vlo=np.concatenate([x[3] for x in q]);vhi=np.concatenate([x[4] for x in q]);rlo=np.concatenate([x[5] for x in q]);rhi=np.concatenate([x[6] for x in q]);w=q[0][7]
  d=run/'calibrated_review'/name;d.mkdir(parents=True)
  unit='V' if name=='calibration_denominator_D' else '1 (dimensionless ellipticity)'
  plots[f'calibrated_review/{name}']=plot_raw_cleaned_envelope(i,lo,hi,vlo,vhi,rlo,rhi,20000,d/'overview.png',f'{TAG}: {name} versus time',unit,valid_count,n,w)
  plots[f'calibrated_review-boundary/{name}']=plot_no_invalid_boundary(d/'boundary_zoom.png',f'{TAG}: full-resolution boundary zoom {name}',unit) if valid_count==n else {'status':'boundary zoom required for invalid data'}
 calibration_provenance=calibration_profile(calibration_parameters)
 (run/'calibration_profile.json').write_text(json.dumps(calibration_provenance,indent=2,sort_keys=True)+'\n')
 (run/'plot_metadata.json').write_text(json.dumps({'label':TAG,'plots':plots,'source_identifier':str(args.input),'selected_segment':meta,'calibration_profile':calibration_provenance},indent=2,sort_keys=True)+'\n');(run/'report.json').write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n');print(json.dumps({'run_directory':str(run),'report':meta},indent=2))
if __name__=='__main__':main()
