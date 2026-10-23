"""Non-destructive preflight for native APE synthetic pipeline execution.

Runs on Mac or VIPER without executing compiled C and without writing files.
Ensures binary, vendor and science configuration are present before work.
"""
import argparse
import json
import os
import platform
from pathlib import Path
import h5py
from tools.run_ape_pipeline import sha256, validate_input


def check(config):
    failures=[]
    passed={}
    def inspect(label, fn):
        try:
            passed[label]=fn()
        except Exception as exc:
            failures.append({"check":label,"error":f"{type(exc).__name__}: {exc}"})
    inspect("synthetic_input",lambda: str(validate_input(config)))
    def binary():
        path=Path(config["ps_executable"]).expanduser().resolve(strict=True)
        if not path.is_file() or not os.access(path,os.X_OK):
            raise ValueError("Missing/nonexecutable C binary")
        digest=sha256(path)
        if digest.lower()!=config["ps_sha256"].lower():
            raise ValueError("Executable SHA256 does not match config")
        return {"path":str(path),"sha256":digest,"host_system":platform.system()}
    inspect("native_c_executable",binary)
    def vendor():
        root=Path(config["external_vendor_root"]).expanduser().resolve(strict=True)
        parts=("scalar_dark_matter_lpsd/lvk/stats.py",
               "scalar_dark_matter_lpsd/lvk/fit_background.py")
        return {name:sha256(root/name) for name in parts}
    inspect("original_vendor_files",vendor)
    def shape():
        if config.get("peak_shape_file") is None:
            return "original vendor default [1]"
        import numpy as np
        path=Path(config["peak_shape_file"]).expanduser().resolve(strict=True)
        a=np.load(path,allow_pickle=False)
        if a.ndim!=1 or not a.size or not np.isfinite(a).all() or np.any(a<0) or not np.isclose(a.sum(),1.,rtol=1e-12,atol=1e-12):
            raise ValueError("Not a normalized finite vendor peak-shape vector")
        return {"path":str(path),"bins":len(a),"sha256":sha256(path)}
    inspect("vendor_peak_shape",shape)
    def output():
        p=Path(config["run_dir"]).expanduser().resolve()
        if p.exists() or not p.name.startswith("ape-pipeline-"):
            raise ValueError("Run directory must be new and named ape-pipeline-*")
        if p.parent.is_dir() and not os.access(p.parent,os.W_OK):
            raise ValueError("No write permission in run parent")
        return str(p)
    inspect("fresh_output",output)
    return {"ready_for_attempt":not failures,"passed":passed,"failures":failures,
            "note":"Preflight is not a test of C numerical output or scientific inference"}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--config",required=True)
    a=p.parse_args()
    outcome=check(json.loads(Path(a.config).read_text()))
    print(json.dumps(outcome,indent=2))
    if not outcome["ready_for_attempt"]:
        raise SystemExit(2)


if __name__=="__main__":
    main()
