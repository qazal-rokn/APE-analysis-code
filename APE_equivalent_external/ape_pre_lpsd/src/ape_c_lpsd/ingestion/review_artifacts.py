"""Future non-overwriting review artifacts; plotting uses labelled display envelopes only."""
from __future__ import annotations
import json,uuid
from pathlib import Path
import numpy as np

def envelope_decimate(values: np.ndarray, max_points: int) -> tuple[np.ndarray,np.ndarray]:
    """Return min/max display envelopes; never a data product or input to cleaning."""
    x=np.asarray(values); width=max(1,int(np.ceil(x.size/max_points))); padded=np.pad(x,(0,(-x.size)%width),constant_values=np.nan); blocks=padded.reshape(-1,width)
    return np.nanmin(blocks,axis=1),np.nanmax(blocks,axis=1)

def write_review_summary(*, output_directory:Path, payload:dict) -> Path:
    output_directory.mkdir(parents=True,exist_ok=True); path=output_directory/f"ingestion-review-{uuid.uuid4().hex}.json"
    if path.exists(): raise FileExistsError("refusing to overwrite review artifact")
    path.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n"); return path

def save_overview_envelope_plot(*, values:np.ndarray, sample_rate_hz:float, max_points:int, title:str, output_path:Path) -> None:
    """Future display-only raw/cleaned overview. Min/max envelope is labelled."""
    import matplotlib.pyplot as plt
    lo,hi=envelope_decimate(values,max_points); width=max(1,int(np.ceil(len(values)/max_points))); time=(np.arange(lo.size)*width)/sample_rate_hz
    fig,ax=plt.subplots(); ax.fill_between(time,lo,hi,alpha=.45,label="display min/max envelope; not analysis resampling")
    ax.set(title=title,xlabel="Time [s]",ylabel="Recorded value"); ax.legend(); fig.savefig(output_path,dpi=160); plt.close(fig)
