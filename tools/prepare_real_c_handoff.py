"""One-command, guarded handoff of the completed real APE calibration to Docker.

Pure preflight by default; --transfer copies once, then verifies the bytes.
Does not run native C or vendor inference, and never overwrites a target.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

from tools.ape_one_run import verify_completed_real_calibration

FROZEN_PS_SHA256 = "cb0b1df10c9863145ba8b35df60e075eac3dc08965f3ad832857a7d4fec5eec9"
PS_BINARY = "/work/ape-ps-method2-20261021/build/lpsd-exec"
DOCKER_ROOT = "/work/ape-real-calibrated-input"


def docker(*args, check=True):
    return subprocess.run(["docker", *args], check=check, text=True, capture_output=True)


def inspection(container, path, source):
    script = """import hashlib,json,os,shutil,sys
p=sys.argv[1];dst=sys.argv[2]
h=hashlib.sha256()
with open(p,'rb') as file:
 for block in iter(lambda:file.read(1024*1024),b''):h.update(block)
result=dict(binary=p,binary_size=os.stat(p).st_size,binary_sha256=h.hexdigest(),
            binary_executable=os.access(p,os.X_OK),
            destination=dst,destination_exists=os.path.lexists(dst),
            destination_parent_exists=os.path.isdir(os.path.dirname(dst)),
            free_bytes=shutil.disk_usage(os.path.dirname(dst) if os.path.isdir(os.path.dirname(dst)) else '/work').free)
if os.path.isfile(dst):
 result['destination_size']=os.stat(dst).st_size
print(json.dumps(result))
"""
    res = docker("exec", container, "python3", "-c", script, PS_BINARY, path)
    return json.loads(res.stdout)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage-result",type=Path,required=True)
    p.add_argument("--container",default="lpsd-linux")
    p.add_argument("--transfer",action="store_true",help="Explicit 9+ GiB Docker copy; no C execution")
    args=p.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",args.container):
        raise ValueError("Invalid Docker container name")
    record = verify_completed_real_calibration(args.stage_result)
    source=Path(record["calibrated_hdf5"])
    if source.name!="calibrated_ape.h5":
        raise ValueError("Unexpected calibrated HDF5 filename")
    run_name=source.parent.name
    if not re.fullmatch(r"ape-run-[A-Za-z0-9_.-]+",run_name):
        raise ValueError("Unexpected run identifier")
    destination=f"{DOCKER_ROOT}/{run_name}/calibrated_ape.h5"
    stat=inspection(args.container,destination,source)
    if stat["binary_sha256"]!=FROZEN_PS_SHA256 or stat["binary_size"]<=1024 or not stat["binary_executable"]:
        raise RuntimeError("Frozen PS executable SHA256/permissions mismatch; STOP")
    required=os.stat(source).st_size + 512*1024*1024
    preflight=dict(status="ready_for_explicit_transfer" if not stat["destination_exists"] else "already_present_verify_before_reuse",
                   mac_source=str(source),mac_size=os.stat(source).st_size,
                   container=args.container,docker_path=destination,
                   verified_frozen_C_binary=PS_BINARY,binary_sha256=FROZEN_PS_SHA256,
                   docker_free_bytes=stat["free_bytes"],
                   docker_destination_exists=stat["destination_exists"],
                   transfer_requested=bool(args.transfer),native_C_executed=False,
                   physical_limits_authorized=False)
    if stat["destination_exists"]:
        raise FileExistsError(f"Refusing to overwrite any existing target: {destination}")
    if stat["free_bytes"] < required:
        raise RuntimeError(f"Docker storage too small: need {required} bytes")
    print(json.dumps(preflight,indent=2),flush=True)
    if not args.transfer:
        print("APE_REAL_C_HANDOFF_PREFLIGHT_PASS_NO_COPY_NO_C")
        return
    # docker cp supports only copying into a path inside a running container.
    # Unique target folder; if copy fails preserve it for diagnosis (never auto delete).
    parent=str(Path(destination).parent)
    docker("exec",args.container,"mkdir","-p",parent)
    if inspection(args.container,destination,source)["destination_exists"]:
        raise FileExistsError("Target appeared during preflight; refusing overwrite")
    print("COPYING_APPROVED_HDF5_ONCE",flush=True)
    subprocess.run(["docker","cp",str(source),f"{args.container}:{destination}"],check=True)
    # Verify full digest against immutable stage record, from the Docker copy.
    script = """import hashlib,json,os,sys
p=sys.argv[1]; h=hashlib.sha256()
with open(p,'rb') as f:
 for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
print(json.dumps(dict(path=p,sha256=h.hexdigest(),size=os.stat(p).st_size)))
"""
    copied=json.loads(docker("exec",args.container,"python3","-c",script,destination).stdout)
    if copied["sha256"]!=record["hdf5_sha256_recorded"] or copied["size"]!=os.stat(source).st_size:
        raise RuntimeError("Docker copied file digest/length mismatch: QUARANTINE; do not execute C")
    print(json.dumps(dict(status="real_calibrated_HDF5_transferred_verified",
        docker_hdf5=destination,verified_sha256=copied["sha256"],
        frozen_C_binary=PS_BINARY,frozen_C_binary_sha256=FROZEN_PS_SHA256,
        native_C_executed=False,physical_limits_authorized=False),indent=2))
    print("APE_REAL_C_INPUT_TRANSFER_VERIFIED_NO_C")


if __name__=="__main__":
    main()
