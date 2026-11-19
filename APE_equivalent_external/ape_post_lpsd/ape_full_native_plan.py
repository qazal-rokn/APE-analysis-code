#!/usr/bin/env python3
"""Print exact native four-hour APE likelihood job-array size; no computation."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import h5py

def plan(h5, fmin, fmax, half_width, task_bins):
    with h5py.File(h5,"r") as h:
        attrs=h["logPSD"].attrs
        n=len(h["frequency"])
        f0=float(attrs["fmin"]);fn=float(attrs["fmax"])
        jdes=int(attrs["Jdes"])
    if n!=jdes or not (0<f0<fmin<fmax<fn):
        raise ValueError("Native spectrum configuration invalid")
    step=math.log(fn/f0)/(jdes-1)
    first=max(half_width,int(math.ceil(math.log(fmin/f0)/step-1e-9)))
    last=min(n-half_width-1,int(math.floor(math.log(fmax/f0)/step+1e-9)))
    count=last-first+1
    if count<=0 or task_bins<=0:
        raise ValueError("Empty mass grid / invalid task-bins")
    tasks=(count+task_bins-1)//task_bins
    return dict(first_index=first,last_index=last,
        masses=count,task_bins=task_bins,tasks=tasks,
        slurm_array=f"0-{tasks-1}",native_log_spacing=step)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spectrum",required=True,type=Path)
    p.add_argument("--fmin",type=float,default=65.1)
    p.add_argument("--fmax",type=float,default=8000.)
    p.add_argument("--half-width-bins",type=int,default=40)
    p.add_argument("--task-bins",type=int,default=10000)
    a=p.parse_args()
    print(json.dumps(plan(a.spectrum,a.fmin,a.fmax,a.half_width_bins,
                          a.task_bins),indent=2))
if __name__=="__main__":
    main()
