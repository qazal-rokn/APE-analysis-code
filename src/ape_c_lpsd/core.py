from __future__ import annotations

import hashlib, json, math, os, re, shutil, subprocess, sys, tempfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import yaml


class ValidationError(RuntimeError): pass

@dataclass(frozen=True)
class Segment:
    start: int; stop: int; blocks: tuple[int, ...]
    @property
    def count(self): return self.stop - self.start

def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path).expanduser().resolve()
    if not path.is_file(): raise ValidationError(f"run configuration does not exist: {path}")
    with path.open() as f: value = yaml.safe_load(f)
    if not isinstance(value, dict): raise ValidationError("run configuration must be a YAML mapping")
    for key in ("cache_directory", "output_directory", "lpsd_executable", "channels", "lpsd"):
        if key not in value: raise ValidationError(f"configuration is missing required key: {key}")
    if value["lpsd"].get("minimum_averages") != 1: raise ValidationError("approved profile requires lpsd.minimum_averages: 1")
    return value

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""): h.update(part)
    return h.hexdigest()

def cache_metadata(cache: Path) -> dict[str, Any]:
    p = cache / "cache_metadata.json"
    if not p.is_file(): raise ValidationError(f"cache metadata missing: {p}")
    with p.open() as f: d = json.load(f)
    if not isinstance(d.get("blocks"), list) or not d["blocks"]: raise ValidationError("cache metadata has no blocks")
    if float(d.get("sample_rate_hz", 0)) <= 0: raise ValidationError("cache metadata has no positive sample_rate_hz")
    return d

def inspect_cache(cache_dir: str | Path, channels: list[str]) -> dict[str, Any]:
    cache = Path(cache_dir).expanduser().resolve(); d = cache_metadata(cache)
    previous = -1
    for b in d["blocks"]:
        if int(b["start_idx"]) != previous + 1 and previous >= 0:
            # Index gaps are allowed cache facts; they must terminate a segment.
            pass
        previous = int(b["stop_idx"]) - 1
        for ch in channels:
            key = f"{ch}_file"
            if key not in b: raise ValidationError(f"block {b['block_id']} has no {key}")
            p = cache / b[key]
            if not p.is_file(): raise ValidationError(f"missing cache channel file: {p}")
        vp = cache / b["valid_file"]
        if not vp.is_file(): raise ValidationError(f"missing cache mask file: {vp}")
    return {"cache": str(cache), "sample_rate_hz": d["sample_rate_hz"], "blocks": len(d["blocks"]), "channels": channels, "metadata_sha256": sha256(cache / "cache_metadata.json")}

def longest_segment(cache_dir: str | Path) -> tuple[dict[str, Any], Segment]:
    cache = Path(cache_dir).expanduser().resolve(); d = cache_metadata(cache)
    best: Segment | None = None; active_start: int | None = None; active_blocks: list[int] = []; expected: int | None = None
    for b in d["blocks"]:
        start, stop, bid = int(b["start_idx"]), int(b["stop_idx"]), int(b["block_id"])
        valid = np.load(cache / b["valid_file"], mmap_mode="r")
        if valid.ndim != 1 or len(valid) != stop-start: raise ValidationError(f"mask shape mismatch in block {bid}")
        # A metadata index discontinuity is a hard boundary.
        if expected is not None and start != expected:
            if active_start is not None:
                candidate=Segment(active_start, expected, tuple(active_blocks)); best=candidate if best is None or candidate.count>best.count else best; active_start=None; active_blocks=[]
        for i, ok in enumerate(valid):
            g = start+i
            if bool(ok):
                if active_start is None: active_start=g; active_blocks=[bid]
                elif not active_blocks or active_blocks[-1] != bid: active_blocks.append(bid)
            elif active_start is not None:
                candidate=Segment(active_start,g,tuple(active_blocks)); best=candidate if best is None or candidate.count>best.count else best; active_start=None; active_blocks=[]
        expected=stop
    if active_start is not None:
        candidate=Segment(active_start,expected,tuple(active_blocks)); best=candidate if best is None or candidate.count>best.count else best
    if best is None: raise ValidationError("cache contains no valid samples")
    return d,best

def stream_segment(cache_dir: str | Path, metadata: dict[str, Any], segment: Segment, channel: str) -> np.ndarray:
    cache=Path(cache_dir); chunks=[]
    for b in metadata["blocks"]:
        a,z=int(b["start_idx"]),int(b["stop_idx"]); lo=max(a,segment.start); hi=min(z,segment.stop)
        if lo>=hi: continue
        arr=np.load(cache/b[f"{channel}_file"], mmap_mode="r"); mask=np.load(cache/b["valid_file"], mmap_mode="r")
        part=np.asarray(arr[lo-a:hi-a], dtype=np.float64)
        if not np.all(mask[lo-a:hi-a]) or not np.all(np.isfinite(part)): raise ValidationError("selected segment changed or is not fully valid/finite")
        chunks.append(part)
    out=np.concatenate(chunks) if chunks else np.empty(0,dtype=np.float64)
    if len(out)!=segment.count: raise ValidationError("selected data order/count mismatch")
    return out

def derived_grid(n: int, fs: float, q: float, sbin: float = 1.0, overlap: float = 0.0, minimum_averages: int = 1) -> tuple[float,float,int]:
    if n < 2: raise ValidationError("segment must contain at least two samples")
    one_cycle=fs/n
    c_minimum=sbin*one_cycle*(1+(1-overlap/100)*(minimum_averages-1))
    return c_minimum, np.nextafter(fs/2, 0.0), math.ceil(math.log(n/2)/math.log(1+1/q))

def preflight_nodes(n: int, fs: float, fmin: float, fmax: float, jdes: int, settings: dict[str, Any]) -> None:
    """Mirror C calc_params: reject, never alter, infeasible requested nodes."""
    if jdes < 2: raise ValidationError("C Jdes must be at least 2")
    g=math.log(fmax/fmin); step=math.exp(g/(jdes-1))-1.0
    for index in range(jdes):
        frequency=fmin*math.exp(index*g/(jdes-1)); required=round(fs/(frequency*step))
        if required > n:
            raise ValidationError("C-LPSD feasibility preflight refused before HDF5/C execution: " + json.dumps({"first_failing_index":index,"first_failing_frequency_hz":frequency,"requested_Jdes":jdes,"required_dft_window_samples":required,"available_segment_samples":n,"settings":settings},sort_keys=True))

def export_hdf5(samples: np.ndarray, path: Path, dataset: str, fs: float, channel: str, unit: str) -> str:
    if path.exists(): raise ValidationError(f"refusing to overwrite HDF5 input: {path}")
    with h5py.File(path,"w") as f:
        ds=f.create_dataset(dataset, data=samples, dtype="f8"); ds.attrs.update(sample_rate_hz=fs, channel=channel, unit=unit, sample_count=len(samples))
    with h5py.File(path,"r") as f:
        ds=f[dataset]
        if ds.dtype != np.dtype("float64") or ds.shape != (len(samples),) or not np.array_equal(ds[:],samples): raise ValidationError("HDF5 export verification failed")
    return sha256(path)

def c_config(path: Path, fmt: str) -> None:
    path.write_text(f'GNUTERM "APE C-LPSD"\nFORMAT "{fmt}"\nEND\n')

def probe_sbin(exe: Path, cfg: dict[str,Any]) -> float:
    """Ask the actual C binary for the effective window SBIN using synthetic data."""
    l=cfg["lpsd"]
    with tempfile.TemporaryDirectory(prefix="ape-c-lpsd-sbin-") as d:
        d=Path(d); h=d/"probe.h5"; export_hdf5(np.zeros(128),h,"/strain",128.,"synthetic","arb")
        out=d/"probe.txt"; c_config(d/"probe.cfg",str(l["output_format"]))
        env={"LPSDCFN":str(d/"probe.cfg"),**os.environ}
        argv=[str(exe),"-d","-D","/strain","-f","128","-h","0","-i",str(h),"-o",str(out),"-b","0","-e",str(127/128),"-s","32","-t","48","-J","2","-n","2","-N","0","-w",str(l["window"]),"-p",str(l["psll_db"]),"-r",str(l["linear_regression"]),"-l",str(l["overlap_percent"]),"-a",str(l["desired_averages"]),"-m",str(l["minimum_averages"]),"-x",str(l["input_scale"]),"-M",str(l["max_memory_power"])]
        r=subprocess.run(argv,text=True,capture_output=True,env=env)
        if r.returncode: raise ValidationError(f"C SBIN probe failed: {r.stderr or r.stdout}")
        found=re.search(r"SBIN \(bins\):\s*([0-9.eE+-]+)",r.stdout)
        if not found: raise ValidationError("C SBIN probe did not report SBIN")
        return float(found.group(1))

def invoke(exe: Path, h5: Path, out: Path, cfg: dict[str,Any], fs: float, n:int, dataset:str, batch:int=0, count:int|None=None) -> list[str]:
    l=cfg["lpsd"]; sbin=probe_sbin(exe,cfg); fmin,fmax,jdes=derived_grid(n,fs,float(l["q_axion"]),sbin,float(l["overlap_percent"]),int(l["minimum_averages"])); count=count or int(l["bins_per_batch"])
    preflight_nodes(n,fs,fmin,fmax,jdes,{"window":l["window"],"psll_db":l["psll_db"],"overlap_percent":l["overlap_percent"],"minimum_averages":l["minimum_averages"],"sbin":sbin})
    argv=[str(exe),"-d","-D",dataset,"-f",repr(fs),"-h","0","-i",str(h5),"-o",str(out),"-b","0","-e",repr((n-1)/fs),"-s",repr(fmin),"-t",repr(fmax),"-J",str(jdes),"-n",str(count),"-N",str(batch),"-w",str(l["window"]),"-p",str(l["psll_db"]),"-r",str(l["linear_regression"]),"-l",str(l["overlap_percent"]),"-a",str(l["desired_averages"]),"-m",str(l["minimum_averages"]),"-x",str(l["input_scale"]),"-M",str(l["max_memory_power"])]
    env={"LPSDCFN":str(out.with_suffix(".cfg")), "MPLCONFIGDIR":tempfile.gettempdir(), **os.environ}; c_config(Path(env["LPSDCFN"]),str(l["output_format"]))
    r=subprocess.run(argv,text=True,capture_output=True,env=env)
    (out.with_suffix(".stdout.txt")).write_text(r.stdout); (out.with_suffix(".stderr.txt")).write_text(r.stderr)
    if r.returncode: raise ValidationError(f"lpsd-exec failed ({r.returncode}); see {out.with_suffix('.stderr.txt')}")
    return argv

def parse_output(path: Path) -> np.ndarray:
    rows=[]
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip(): continue
        vals=line.split()
        if len(vals)!=5: raise ValidationError(f"unexpected C output columns in {path}")
        rows.append([float(v) for v in vals])
    a=np.asarray(rows,float)
    if not len(a) or not np.all(np.isfinite(a)) or np.any(np.diff(a[:,0])<=0) or np.any(a[:,1]<0): raise ValidationError("invalid C spectrum output")
    return a
