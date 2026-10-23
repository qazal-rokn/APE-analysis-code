"""Opt-in transport for previously validated experimental C-LPSD direct-PS binary.

No PS is estimated in Python. Original C algorithms and planner remain unchanged.
"""
from __future__ import annotations
import copy
import hashlib
import math
import os
from pathlib import Path


class PSError(ValueError):
    pass


def plan_ps_optin(legacy_plan, *, ps_executable, expected_sha256):
    """Copy an existing PSD C plan, use verified PS binary and 7-column format."""
    binary = Path(ps_executable).expanduser().resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise PSError("PS executable missing or not executable")
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    if (not isinstance(expected_sha256, str) or len(expected_sha256) != 64
            or digest.lower() != expected_sha256.lower()):
        raise PSError("PS C executable SHA-256 mismatch")
    if (legacy_plan.get("format") != "fDNbri"
            or 'FORMAT "fDNbri"' not in legacy_plan.get("c_config", "")):
        raise PSError("Expected the unchanged historical six-column plan")
    plan = copy.deepcopy(legacy_plan)
    plan["fft"][0] = str(binary)
    for _, argv in plan["batches"]:
        argv[0] = str(binary)
        out_idx = argv.index("-o")+1
        argv[out_idx] = argv[out_idx]+".raw7"
    plan["c_config"] = plan["c_config"].replace(
        'FORMAT "fDNbri"', 'FORMAT "fDSNbri"')
    plan["format"] = "fDSNbri"
    plan["ps_enabled"] = True
    plan["ps_binary_sha256"] = digest
    return plan


def split_ps_batch(raw7_path, legacy6_path, ps_path, *, expected_rows):
    """Stream seven C columns into unchanged legacy six plus (frequency,PS).

    Preserve original numeric tokens. Exclusive partials; no target overwrites.
    Output six fields: f,PSD,Navg,Nkernel,Re,Im; two fields: f,PS.
    """
    raw, legacy, ps = (Path(v) for v in (raw7_path, legacy6_path, ps_path))
    if len({raw.resolve(), legacy.resolve(), ps.resolve()}) != 3:
        raise PSError("Input and output paths must be distinct")
    parts = [Path(str(legacy)+".partial"), Path(str(ps)+".partial")]
    if any(p.exists() for p in (legacy, ps, *parts)):
        raise PSError("Output already exists")
    if not isinstance(expected_rows, int) or expected_rows <= 0:
        raise PSError("Positive expected_rows required")
    count, previous = 0, None
    try:
        with raw.open(encoding="utf-8") as inp, parts[0].open("x") as a, parts[1].open("x") as b:
            for lineno, line in enumerate(inp, 1):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.split()
                if len(fields) != 7:
                    raise PSError(f"Invalid seven-column row at {lineno}")
                try:
                    nums = [float(x) for x in fields]
                except ValueError as exc:
                    raise PSError(f"Nonnumeric row at {lineno}") from exc
                if (not all(math.isfinite(x) for x in nums)
                        or nums[0] <= 0 or nums[1] < 0 or nums[2] < 0
                        or nums[3] < 0 or nums[4] < 0
                        or (previous is not None and nums[0] <= previous)):
                    raise PSError(f"Invalid PS/PSD or frequency at {lineno}")
                count += 1
                if count > expected_rows:
                    raise PSError("Too many output rows")
                previous = nums[0]
                a.write(" ".join((fields[0], fields[1], *fields[3:]))+"\n")
                b.write(fields[0]+" "+fields[2]+"\n")
        if count != expected_rows:
            raise PSError(f"Expected {expected_rows} rows; found {count}")
        os.link(parts[0], legacy)
        try:
            os.link(parts[1], ps)
        except BaseException:
            legacy.unlink(missing_ok=True)
            raise
        return count
    finally:
        for p in parts:
            p.unlink(missing_ok=True)
