"""One-shot host-to-Docker handoff for the already calibrated real APE segment.

Read-only preflight by default. --copy-hdf5 explicitly copies into the
existing lpsd-linux container, verifies complete SHA-256 there, and leaves
original files untouched. NEVER runs native C or vendor likelihood. This is
a transport stage, not a real-data inference approval.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

from tools.ape_one_run import verify_completed_real_calibration

EXPECTED_BINARY_SHA = "cb0b1df10c9863145ba8b35df60e075eac3dc08965f3ad832857a7d4fec5eec9"
BINARY = "/work/ape-ps-method2-20261021/build/lpsd-exec"
TARGET_DIR = "/work/ape-real-handoff-c0284-v1"
TARGET = TARGET_DIR + "/calibrated_ape.h5"


def docker(container, *args):
    return subprocess.run(["docker", *args], capture_output=True, text=True, check=True)


def exec_python(container, script, *args):
    res = docker(container, "exec", container, "python3", "-c", script, *args)
    return json.loads(res.stdout)


def check_container(container):
    return exec_python(container, """import os,hashlib,json,sys
p=sys.argv[1]
with open(p,'rb') as f: sha=hashlib.sha256(f.read()).hexdigest()
print(json.dumps(dict(path=p,sha256=sha,executable=os.access(p,os.X_OK))))
""", BINARY)


def check_target(container):
    return exec_python(container, """import os,json,sys
p=sys.argv[1]
print(json.dumps(dict(exists=os.path.lexists(p),is_file=os.path.isfile(p),
size_bytes=os.stat(p).st_size if os.path.isfile(p) else None)))
""", TARGET)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage-result",type=Path,required=True)
    ap.add_argument("--container",default="lpsd-linux")
    ap.add_argument("--copy-hdf5",action="store_true")
    a=ap.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",a.container):
        raise ValueError("Unsafe container name")
    src=verify_completed_real_calibration(a.stage_result)
    h5=Path(src["calibrated_hdf5"])
    binary=check_container(a.container)
    if binary["sha256"]!=EXPECTED_BINARY_SHA or not binary["executable"]:
        raise RuntimeError("Frozen validated C executable missing or changed; refusing copy")
    target=check_target(a.container)
    available=exec_python(a.container, """import os,json
s=os.statvfs("/work")
print(json.dumps(dict(available_bytes=s.f_bavail*s.f_frsize)))
""")
    if h5.stat().st_size + 512*1024*1024 > available["available_bytes"]:
        raise RuntimeError("Insufficient free Docker storage for a safe copy")
    report=dict(status="handoff_plan_only",source_hdf5=str(h5),
                source_sha256_recorded=src["hdf5_sha256_recorded"],
                docker_container=a.container,docker_target=TARGET,
                validated_native_binary_sha256=binary["sha256"],
                destination_preexists=target["exists"],
                docker_free_bytes=available["available_bytes"],native_c_executed=False,
                vendor_modified=False,physical_limits_authorized=False)
    if not a.copy_hdf5:
        print(json.dumps(report,indent=2))
        print("HANDOFF_PREFLIGHT_PASS_NO_DATA_COPIED")
        return
    if target["exists"]:
        raise FileExistsError("Refusing existing Docker destination: "+TARGET)
    # Copy once, not per frequency or quadrature. The completed source is
    # immutable. A new dedicated directory is required; never remove/reuse.
    docker(a.container,"exec",a.container,"mkdir",TARGET_DIR)
    # Docker cp creates its destination file; no shell or pipe injection.
    subprocess.run(["docker","cp",str(h5),a.container+":"+TARGET],check=True)
    verified=exec_python(a.container,"""import hashlib,json,os,sys
p=sys.argv[1]
h=hashlib.sha256()
with open(p,'rb') as f:
 for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
print(json.dumps(dict(path=p,sha256=h.hexdigest(),bytes=os.stat(p).st_size)))
""",TARGET)
    if verified["sha256"] != src["hdf5_sha256_recorded"]:
        raise RuntimeError("Docker copy did not match recorded source SHA-256")
    report.update(status="verified_hdf5_handoff_complete_native_C_not_executed",
                  docker_copy=verified)
    print(json.dumps(report,indent=2))
    print("HANDOFF_COPY_VERIFIED_C_NOT_EXECUTED")


if __name__=="__main__":
    main()
