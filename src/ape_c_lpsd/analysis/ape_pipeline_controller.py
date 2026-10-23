"""Synthetic-only APE PS execution plan and guarded postprocessing.

This controller cannot manufacture or deploy the external PS C binary.
An exact operator supplied executable SHA256 is required. No real-data mode,
no hidden C execution, no vendor modification, no overwriting any run output.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ape_c_lpsd.spectral.constantq import plan_commands
from ape_c_lpsd.spectral.ps_optin import plan_ps_optin, split_ps_batch, PSError
from ape_c_lpsd.backends.local import run_local


def make_synthetic_ps_plan(*, original_executable, ps_executable,
                           ps_binary_sha256, synthetic_hdf5, dataset,
                           fft_hdf5, output_dir):
    """Frozen production profile; for small fixtures use existing fixture plan."""
    original=plan_commands(original_executable,synthetic_hdf5,dataset,fft_hdf5,output_dir)
    return plan_ps_optin(original,ps_executable=ps_executable,
                         expected_sha256=ps_binary_sha256)


def finalize_ps_plan(plan):
    """Finalize native PS raw7 batch files, retaining originals read-only.

    Requires C execution to have succeeded already. Does not join across
    partitions: the original six-column batch combiner stays authoritative.
    """
    if not plan.get("ps_enabled") or plan.get("format")!="fDSNbri":
        raise PSError("Require explicit PS-enabled C plan")
    report=[]
    for index,argv in plan["batches"]:
        out=Path(argv[argv.index("-o")+1])
        start=int(argv[argv.index("-N")+1])
        count=int(argv[argv.index("-n")+1])
        if not out.is_file():
            raise FileNotFoundError(f"Missing completed C batch: {out}")
        report.append((index,start,count,out))
    results=[]
    for index,start,count,out in report:
        base=Path(str(out)[:-5]) if str(out).endswith(".raw7") else None
        if base is None:
            raise PSError("PS output must use .raw7 suffix")
        rows=split_ps_batch(out,base,Path(str(base)+".ps"),expected_rows=count)
        results.append(dict(batch=index,start_index=start,rows=rows,ps_file=str(base)+".ps"))
    return results


def read_candidate_ps_window(ps_files_by_quadrature, *, start_index, bin_count):
    """Read identical candidate windows of already finalized native C PS."""
    if start_index<0 or bin_count<=0:
        raise ValueError("Invalid candidate window index/count")
    freqs, values=[],[]
    for label,path in ps_files_by_quadrature.items():
        data=np.loadtxt(path,ndmin=2)
        if data.ndim!=2 or data.shape[1]!=2:
            raise ValueError(f"{label}: expected frequency,PS columns")
        end=start_index+bin_count
        if end>len(data):
            raise ValueError("Candidate window extends beyond PS data")
        part=data[start_index:end]
        if (not np.isfinite(part).all() or np.any(part[:,0]<=0)
                or np.any(part[:,1]<=0) or np.any(np.diff(part[:,0])<=0)):
            raise ValueError(f"{label}: invalid candidate PS")
        freqs.append(part[:,0])
        values.append(part[:,1])
    if not freqs or any(not np.array_equal(f,freqs[0]) for f in freqs[1:]):
        raise ValueError("Quadrature grids do not match")
    return freqs[0],np.stack(values,axis=0)
