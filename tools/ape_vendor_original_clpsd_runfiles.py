"""Use the ORIGINAL vendor create_lpsd_runfile.create_runfile for X and Y.

This only supplies two APE input channels and separate output paths to the
unchanged vendor runfile generator. Original C executable, templates, method,
batching and execution semantics are retained. It generates scripts but does
not submit or execute them without explicit user action.

The only new spectral operation after execution is linear PS_X + PS_Y.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

def generate(config):
    from scalar_dark_matter_lpsd.create_lpsd_runfile import create_runfile
    cfg=json.loads(Path(config).resolve(strict=True).read_text())
    base=dict(cfg["vendor_lpsd_args"])
    required=("input_file","fsample","fmin","resolution","method",
              "path_to_lpsd_exec","data_length")
    missing=[k for k in required if k not in base]
    if missing:
        raise ValueError(f"Missing original vendor C-LPSD parameters: {missing}")
    out=Path(cfg["runfile_directory"]).resolve()
    if not out.is_dir():
        raise NotADirectoryError(out)
    scripts=[]
    for channel in ("X","Y"):
        args=base.copy()
        args["channel"]=cfg["quadrature_channels"][channel]
        stem=cfg["runfile_prefix"]+"_"+channel
        args["run_file"]=str(out/(stem+".sh"))
        args["output_file"]=str(out/(stem+".lpsd"))
        args["fft_file"]=str(out/(stem+"_fft.h5"))
        args["submit_file"]=str(out/(stem+".submit"))
        args["outdir"]=str(out)
        if any(Path(args[k]).exists() for k in ("run_file","output_file","fft_file","submit_file")):
            raise FileExistsError("Refusing to overwrite an existing vendor run artifact")
        args.setdefault("use_condor",False)
        args.setdefault("fmax",None)
        args.setdefault("Jdes",None)
        args.setdefault("n_batches",1)
        args.setdefault("block_epsilon",10)
        args.setdefault("constQ_epsilon",1e-15)
        args.setdefault("max_memory_power",26)
        args.setdefault("request_cpus",1)
        args.setdefault("request_disk_GB",16)
        args.setdefault("request_memory_GB",16)
        args.setdefault("accounting_group","ligo.prod.o4.cw.darkmatter.lpsd")
        args.setdefault("accounting_group_user","")
        args.setdefault("python_executable","python")
        if args["method"] == 5 and not args["use_condor"]:
            raise ValueError("Original vendor method=5 requires Condor; use original method=2/3 locally or configure Condor")
        create_runfile(args)
        scripts.append(args["run_file"])
    return {"status":"GENERATED_NOT_EXECUTED","original_vendor_function":"create_runfile",
            "scripts":scripts,"note":"Run both vendor-generated scripts, then sum X/Y linear PS."}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--config",required=True)
    print(json.dumps(generate(p.parse_args().config),indent=2))
