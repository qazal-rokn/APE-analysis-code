"""Non-overwriting, display-only ingestion-review plots for arbitrary channels/stages."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import hashlib, json, uuid
from pathlib import Path
from typing import Mapping
import numpy as np

class PlottingError(RuntimeError): pass

@dataclass(frozen=True)
class PlotOptions:
    display_max_points: int
    zoom_max_samples: int
    dpi: int
    raw_enabled: bool
    cleaned_enabled: bool
    selected_enabled: bool

@dataclass(frozen=True)
class ChannelPlotInput:
    values: np.ndarray
    sample_rate_hz: float
    source_identifier: str
    unit: str
    valid_mask: np.ndarray | None = None
    selected_start: int | None = None
    selected_stop: int | None = None

def _hash(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for c in iter(lambda:f.read(1<<20),b''): h.update(c)
    return h.hexdigest()

def _envelope(values: np.ndarray, maximum: int) -> tuple[np.ndarray,np.ndarray,np.ndarray,int]:
    x=np.asarray(values,dtype=float)
    if x.ndim!=1 or not len(x): raise PlottingError('plot values must be a non-empty vector')
    width=max(1,int(np.ceil(len(x)/maximum))); padded=np.pad(x,(0,(-len(x))%width),constant_values=np.nan); blocks=padded.reshape(-1,width)
    indices=np.arange(blocks.shape[0])*width
    return indices,np.nanmin(blocks,axis=1),np.nanmax(blocks,axis=1),width

def _save_envelope(*, values:np.ndarray, sample_rate_hz:float, title:str, ylabel:str, destination:Path, options:PlotOptions, mask:np.ndarray|None=None, focus:tuple[int,int]|None=None) -> dict:
    import matplotlib.pyplot as plt
    x=np.asarray(values); offset=0
    if focus is not None:
        start,stop=focus
        if start<0 or stop<=start or stop>len(x): raise PlottingError('invalid focus interval')
        if stop-start>options.zoom_max_samples: raise PlottingError('focused plot exceeds zoom_max_samples; caller must select a bounded full-resolution interval')
        x=x[start:stop]; offset=start; width=1; indices=np.arange(len(x)); lo=hi=x
        label='full-resolution selected/zoom samples; no decimation'
    else:
        indices,lo,hi,width=_envelope(x,options.display_max_points)
        label=f'display min/max envelope, {width} original samples per bin; not analysis resampling'
    fig,ax=plt.subplots(figsize=(10,4.8)); time=(indices+offset)/sample_rate_hz
    ax.fill_between(time,lo,hi,alpha=.45,label=label)
    ax.set(title=title,xlabel='Time [s]',ylabel=ylabel); ax.grid(alpha=.25); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(destination,dpi=options.dpi); plt.close(fig)
    return {'path':str(destination),'sha256':_hash(destination),'display_method':label,'sample_offset':offset,'display_bin_width_samples':width}

def _save_mask(*, mask:np.ndarray, sample_rate_hz:float, title:str, destination:Path, options:PlotOptions) -> dict:
    import matplotlib.pyplot as plt
    m=np.asarray(mask,dtype=bool)
    index,lo,hi,width=_envelope(m.astype(float),options.display_max_points)
    # Mean valid fraction is distinct from min/max envelope and remains display-only.
    padded=np.pad(m.astype(float),(0,(-len(m))%width),constant_values=np.nan); fraction=np.nanmean(padded.reshape(-1,width),axis=1)
    fig,ax=plt.subplots(figsize=(10,3.8)); ax.step(index/sample_rate_hz,fraction,where='post',label=f'valid fraction per {width}-sample display block')
    ax.set(title=title,xlabel='Time [s]',ylabel='Valid fraction',ylim=(-.02,1.02)); ax.grid(alpha=.25); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(destination,dpi=options.dpi); plt.close(fig)
    return {'path':str(destination),'sha256':_hash(destination),'display_method':'valid-fraction display blocks; mask unchanged','display_bin_width_samples':width}

def generate_ingestion_review_plots(*, output_root:Path, config_snapshot:dict, diagnostic_channels:tuple[str,...], raw:Mapping[str,ChannelPlotInput], cleaned:Mapping[str,ChannelPlotInput], selected:Mapping[str,ChannelPlotInput], options:PlotOptions) -> Path:
    """Write one immutable review directory; input arrays are never modified."""
    expected=set(diagnostic_channels)
    if set(raw)!=expected or set(cleaned)!=expected or set(selected)!=expected: raise PlottingError('each configured diagnostic channel requires raw, cleaned, and selected inputs')
    run=output_root/f'ingestion-plots-{uuid.uuid4().hex}'
    if run.exists(): raise FileExistsError(f'refusing overwrite: {run}')
    run.mkdir(parents=True); (run/'config_snapshot.json').write_text(json.dumps(config_snapshot,indent=2,sort_keys=True)+'\n')
    metadata={'kind':'ingestion_review_plots','channels':{},'config_snapshot':'config_snapshot.json','stages':['raw','cleaned','selected_longest_segment']}
    for channel in diagnostic_channels:
        metadata['channels'][channel]={}
        for stage,item,enabled in (('raw',raw[channel],options.raw_enabled),('cleaned',cleaned[channel],options.cleaned_enabled),('selected_longest_segment',selected[channel],options.selected_enabled)):
            if not enabled: continue
            d=run/stage/channel; d.mkdir(parents=True)
            overview=_save_envelope(values=item.values,sample_rate_hz=item.sample_rate_hz,title=f'{stage}: {channel} ({item.source_identifier})',ylabel=item.unit,destination=d/'overview.png',options=options)
            record={'source_identifier':item.source_identifier,'sample_rate_hz':item.sample_rate_hz,'unit':item.unit,'overview':overview}
            if stage=='cleaned' and item.valid_mask is not None: record['quality_mask']=_save_mask(mask=item.valid_mask,sample_rate_hz=item.sample_rate_hz,title=f'cleaning mask: {channel}',destination=d/'mask_fraction.png',options=options)
            if stage=='selected_longest_segment' and item.selected_start is not None and item.selected_stop is not None: record['focus']=_save_envelope(values=item.values,sample_rate_hz=item.sample_rate_hz,title=f'selected segment focus: {channel}',ylabel=item.unit,destination=d/'focus.png',options=options,focus=(item.selected_start,item.selected_stop))
            metadata['channels'][channel][stage]=record
    (run/'plot_metadata.json').write_text(json.dumps(metadata,indent=2,sort_keys=True)+'\n'); return run
