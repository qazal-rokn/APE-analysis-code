"""One-command, fail-closed APE orchestration of existing frozen stage programs.

Scientific methods are supplied as immutable external commands. This module
does not replace or modify cleaning, calibration, C-LPSD or vendor likelihood.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

STAGES = ("ingest_clean_select", "calibrate", "native_c_ps", "background_fit",
          "likelihood_mle", "candidate_search", "expected_upper_limits", "final_report")

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()

def _check(outputs, *, required):
    result = {}
    for item in outputs:
        p = Path(item["path"]).expanduser().resolve()
        if not p.is_file():
            if required:
                raise FileNotFoundError(p)
            result[str(p)] = None
            continue
        got = digest(p)
        expected = item["sha256"]
        if expected is not None and got.lower() != expected.lower():
            raise RuntimeError(f"Immutable source hash mismatch: {p}")
        result[str(p)] = got
    return result

def execute(config_path, *, preflight=False):
    path = Path(config_path).expanduser().resolve(strict=True)
    cfg = json.loads(path.read_text())
    if cfg.get("schema_version") != 1:
        raise ValueError("Expected schema_version=1")
    stages = cfg.get("stages", {})
    if set(stages) != set(STAGES):
        raise ValueError("All eight stages must be configured explicitly")
    directory = Path(cfg["run_dir"]).expanduser().resolve()
    record = directory / "orchestrator_record.json"
    if record.exists():
        raise FileExistsError(f"Already completed or attempted: {record}")
    plan = []
    for name in STAGES:
        s = stages[name]
        if s.get("mode") not in ("reuse", "run"):
            raise ValueError(f"{name}: require reuse or run")
        outputs = s.get("outputs", [])
        if not isinstance(outputs, list) or not outputs:
            raise ValueError(f"{name}: explicit expected outputs required")
        for item in outputs:
            if not isinstance(item.get("path"), str):
                raise ValueError(f"{name}: invalid artifact path")
            sha = item.get("sha256")
            if sha is not None and (not isinstance(sha, str) or len(sha) != 64
                                    or any(x not in "0123456789abcdef" for x in sha.lower())):
                raise ValueError(f"{name}: invalid SHA256")
        before = _check(outputs, required=s["mode"] == "reuse")
        if s["mode"] == "run":
            argv = s.get("argv")
            if not isinstance(argv, list) or not argv or not all(isinstance(x,str) and x for x in argv):
                raise ValueError(f"{name}: exact original executable argv required")
            if any(value is not None for value in before.values()):
                raise FileExistsError(f"{name}: refuse to overwrite existing stage output")
        elif any(item["sha256"] is None for item in outputs):
            raise ValueError(f"{name}: reused artifacts need SHA256")
        plan.append({"stage": name, "mode": s["mode"]})
    if preflight:
        return {"status": "PREFLIGHT_OK_NO_EXECUTION", "plan": plan}
    directory.mkdir(parents=True, exist_ok=True)
    finished = []
    try:
        for item in plan:
            name = item["stage"]
            stage = stages[name]
            if item["mode"] == "run":
                env = dict(os.environ)
                env["PYTHONDONTWRITEBYTECODE"] = "1"
                env.update(stage.get("env", {}))
                subprocess.run(stage["argv"], cwd=stage.get("cwd"), env=env, check=True)
            outputs = _check(stage["outputs"], required=True)
            finished.append({"stage": name, "mode": item["mode"], "sha256": outputs})
    except BaseException as exc:
        report = {"status": "FAILED_NOT_SCIENTIFIC_RESULT", "error": repr(exc),
                  "completed": finished, "physical_ape_limits_validated": False}
        with record.open("x") as f:
            json.dump(report, f, indent=2)
        raise
    report = {"status": "STAGES_COMPLETED_PROVISIONAL", "config_sha256": digest(path),
              "completed": finished, "scientific_algorithms_modified": False,
              "physical_ape_limits_validated": False,
              "global_significance_validated": False,
              "observed_limits_validated": False}
    with record.open("x") as f:
        json.dump(report, f, indent=2)
    return report

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--preflight", action="store_true")
    a = parser.parse_args()
    print(json.dumps(execute(a.config, preflight=a.preflight), indent=2))

if __name__ == "__main__":
    main()
