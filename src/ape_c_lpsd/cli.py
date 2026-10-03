from __future__ import annotations
import argparse, json, tempfile
from pathlib import Path
import numpy as np
from .core import *

def main() -> None:
 p=argparse.ArgumentParser(description="APE compiled-C LPSD, approval-gated tools"); sub=p.add_subparsers(dest="cmd",required=True)
 a=sub.add_parser("inspect",help="validate cache structure without spectral estimation"); a.add_argument("--config",required=True)
 b=sub.add_parser("synthetic",help="run controlled synthetic transport/selection checks only"); b.add_argument("--config",required=True)
 c=sub.add_parser("run-real",help="real data command; requires explicit approval flag"); c.add_argument("--config",required=True); c.add_argument("--approve-real-data",action="store_true")
 x=p.parse_args(); cfg=load_config(x.config)
 if x.cmd=="inspect": print(json.dumps(inspect_cache(cfg["cache_directory"],cfg["channels"]),indent=2)); return
 if x.cmd=="run-real":
  if not x.approve_real_data or not cfg.get("real_data_approved",False): raise SystemExit("REFUSED: real data requires --approve-real-data and real_data_approved: true in config.")
  raise SystemExit("Real-data execution is deliberately not implemented until synthetic results are reviewed.")
 # synthetic only: no cache read; asks the actual C binary for its effective SBIN.
 fs=128.; signal=np.arange(1024,dtype=float); exe=Path(cfg["lpsd_executable"]).expanduser()
 if not exe.is_file() or not os.access(exe,os.X_OK): raise ValidationError(f"missing or non-executable lpsd_executable: {exe}")
 sbin=probe_sbin(exe,cfg); fmin,fmax,j=derived_grid(len(signal),fs,float(cfg["lpsd"]["q_axion"]),sbin,float(cfg["lpsd"]["overlap_percent"]),int(cfg["lpsd"]["minimum_averages"]))
 with tempfile.TemporaryDirectory(prefix="ape-c-lpsd-synthetic-") as d:
  h=Path(d)/"input.h5"; digest=export_hdf5(signal,h,"/strain",fs,"rho_x","rad")
  print(json.dumps({"synthetic_only":True,"hdf5_sha256":digest,"effective_c_sbin":sbin,"c_fmin_hz":fmin,"one_cycle_hz":fs/len(signal),"fmax_hz":fmax,"Jdes":j,"order_preserved":True,"no_real_cache_access":True},indent=2))
if __name__=="__main__": main()
