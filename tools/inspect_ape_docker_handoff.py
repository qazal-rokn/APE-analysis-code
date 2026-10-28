"""Read-only preflight for real APE HDF5 -> Docker C-LPSD handoff."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from tools.ape_one_run import verify_completed_real_calibration

# Hash found in frozen C-LPSD validation execution reports.
FROZEN_PS_C_SHA256 = "cb0b1df10c9863145ba8b35df60e075eac3dc08965f3ad832857a7d4fec5eec9"


def inspect(stage_result, container):
    checked = verify_completed_real_calibration(stage_result)
    source = Path(checked["calibrated_hdf5"]).resolve(strict=True)
    if not container or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for c in container):
        raise ValueError("Invalid container name")
    mounts = json.loads(subprocess.run(
        ["docker", "inspect", "--format", "{{json .Mounts}}", container],
        check=True, capture_output=True, text=True).stdout)
    visible = []
    for mount in mounts:
        if mount.get("Type") != "bind":
            continue
        host = Path(mount["Source"]).resolve()
        if source == host or host in source.parents:
            visible.append(str(Path(mount["Destination"]) / source.relative_to(host)))
    script = """import os,hashlib,json
result=[]
for directory in ('/work/ape-ps-method2-20261021/build','/work/lpsd-project/LPSD','/work/lpsd-project/LPSD/build','/work/lpsd-project/LPSD/install/bin','/source/LPSD','/work/ape-pipeline'):
 if not os.path.isdir(directory): continue
 for name in ('lpsd-exec','lpsd-exec-ps','lpsd','lpsd_ps','lpsd-ps'):
  path=os.path.join(directory,name)
  if os.path.isfile(path):
   with open(path,'rb') as f: data=f.read()
   result.append(dict(path=path,size_bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),is_executable=os.access(path,os.X_OK)))
print(json.dumps(result))
"""
    result = subprocess.run(["docker", "exec", container, "python3", "-c", script],
                            check=True, capture_output=True, text=True)
    candidates=json.loads(result.stdout)
    matched=[b for b in candidates if b["sha256"] == FROZEN_PS_C_SHA256 and b["size_bytes"] > 1024 and b["is_executable"]]
    return dict(status="read_only_handoff_inspection",container=container,
                calibrated_hdf5=str(source),mounted_paths=visible,
                hdf5_mounted=bool(visible),
                candidate_binaries_not_yet_approved=candidates,
                frozen_validation_binary_sha256=FROZEN_PS_C_SHA256,
                exact_frozen_binary_matches=matched,
                matched_hash_not_a_substitute_for_real_data_run_approval=True,
                ps_binary_approved=False,native_c_executed=False,
                source_data_modified=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-result",type=Path,required=True)
    parser.add_argument("--container",default="lpsd-linux")
    args=parser.parse_args()
    print(json.dumps(inspect(args.stage_result,args.container),indent=2))


if __name__=="__main__":
    main()
