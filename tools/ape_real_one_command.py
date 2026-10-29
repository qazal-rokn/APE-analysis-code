"""One-command APE provisional real-data chain, reusing frozen scientific code.

Accepts existing frozen calibration or calls the unchanged calibration stage,
then uses verified native C PS (or runs the same C via Docker) and calls the
original vendor background, likelihood, candidate, and expected-limit stages.
No synthetic provenance bypass, no PSD-to-PS guess, no physical-exclusion claim.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def docker(container, *args, capture=False):
    return subprocess.run(["docker", "exec", container, *map(str, args)],
                          check=True, capture_output=capture, text=True)


def choose_calibration(source, *, dry_run):
    from tools.ape_one_run import verify_completed_real_calibration
    if ("stage_result" in source) == ("run_config" in source):
        raise ValueError("Provide exactly one stage_result or frozen run_config")
    if "stage_result" in source:
        return verify_completed_real_calibration(Path(source["stage_result"]))
    from tools.run_ape_automatic import prepare, execute
    config = Path(source["run_config"]).expanduser().resolve(strict=True)
    _, _, folder, _, _, _ = prepare(config, require_mac=not dry_run)
    report = folder / "stage_result.json"
    if report.is_file():
        return verify_completed_real_calibration(report)
    if dry_run:
        return {"status": "calibration_planned", "stage_result": str(report)}
    execute(config, execute_calibration=True)
    return verify_completed_real_calibration(report)


def reuse_native(native, calibrated):
    if native.get("calibrated_hdf5_sha256") != calibrated["hdf5_sha256_recorded"]:
        raise ValueError("Native C PS cannot be linked to this calibrated HDF5")
    spectra = native.get("spectra", {})
    if set(spectra) != {"rho_x", "rho_y"}:
        raise ValueError("Need both independent calibrated X/Y native PS products")
    paths = {}
    for label, artifact in spectra.items():
        if artifact.get("unit") != "rad^2" or artifact.get("schema") != "frequency_PS_2_columns":
            raise ValueError("Require native C PS, not PSD or ASD")
        path = Path(artifact["path"]).expanduser().resolve(strict=True)
        if sha256(path) != artifact["sha256"]:
            raise ValueError(f"Native PS SHA256 mismatch: {label}")
        paths[label] = str(path)
    return paths


def docker_sha(container, path):
    return docker(container, "sha256sum", path, capture=True).stdout.split()[0]


def docker_exists(container, path):
    return subprocess.run(["docker", "exec", container, "test", "-e", str(path)],
                          capture_output=True).returncode == 0


def ensure_docker_hdf5(container, calibrated, target, *, dry_run=False):
    source = Path(calibrated["calibrated_hdf5"]).resolve(strict=True)
    expected = calibrated["hdf5_sha256_recorded"]
    if docker_exists(container, target):
        if docker_sha(container, target) != expected:
            raise RuntimeError("Existing Docker HDF5 differs: do not overwrite")
        return
    if dry_run:
        return
    docker(container, "mkdir", "-p", str(Path(target).parent))
    if docker_exists(container, target):
        raise FileExistsError("Docker destination appeared during preflight")
    subprocess.run(["docker", "cp", str(source), f"{container}:{target}"], check=True)
    if docker_sha(container, target) != expected:
        raise RuntimeError("Docker transfer checksum mismatch")


def run_native_c(config, calibrated, output_dir, *, dry_run):
    container = config["container"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", container):
        raise ValueError("Invalid Docker container name")
    binary = config["binary"]
    if docker_sha(container, binary) != config["binary_sha256"]:
        raise ValueError("Original C binary hash mismatch")
    settings = config["c_settings"]
    required = {"fs_hz", "duration_s", "f_min_hz", "f_max_hz", "grid_bins",
                "overlap_percent", "window", "psll_db", "regression", "input_scale",
                "epsilon", "memory_power", "kernel_threshold", "a_option", "batch_count"}
    if set(settings) != required:
        raise ValueError("All original C numerical settings must be specified")
    from ape_c_lpsd.spectral.constantq import _c_argv
    if settings["fs_hz"] != calibrated["sample_rate_hz"]:
        raise ValueError("C sampling frequency does not match calibrated input")
    if settings["duration_s"] > calibrated["samples_per_channel"] / calibrated["sample_rate_hz"]:
        raise ValueError("C requests longer observation than calibrated HDF5")
    batches, count = int(settings["batch_count"]), int(settings["grid_bins"])
    if not 1 <= batches <= count:
        raise ValueError("Invalid C frequency partition")
    root = str(Path(config["docker_work_dir"]) / output_dir.parent.name / output_dir.name)
    source_name = Path(calibrated["calibrated_hdf5"]).parent.name
    target = str(Path(config["docker_hdf5_root"]) / source_name / "calibrated_ape.h5")
    if not dry_run:
        ensure_docker_hdf5(container, calibrated, target)
        if docker_exists(container, root):
            raise FileExistsError("C destination exists: refusing overwrite " + root)
        docker(container, "mkdir", "-p", root)
        output_dir.mkdir(parents=True, exist_ok=True)
        config_file = output_dir / "lpsd.cfg"
        config_file.write_text('GNUTERM "APE direct PS"\nFORMAT "fDSNbri"\nEND\n')
        subprocess.run(["docker", "cp", str(config_file), f"{container}:{root}/lpsd.cfg"],
                       check=True)
    spectra = {}
    for channel in ("rho_x", "rho_y"):
        folder = root + "/" + channel
        if not dry_run:
            docker(container, "mkdir", "-p", folder)
        kw = dict(exe=binary, inp=target, dataset="/" + channel,
                  fft=folder + "/fft.h5", fs=settings["fs_hz"],
                  t_stop=settings["duration_s"], f_min=settings["f_min_hz"],
                  f_max=settings["f_max_hz"], grid_bins=count,
                  overlap=settings["overlap_percent"], regression=settings["regression"],
                  window=settings["window"], psll=settings["psll_db"],
                  scale=settings["input_scale"], epsilon=settings["epsilon"],
                  memory_power=settings["memory_power"],
                  kernel_threshold=settings["kernel_threshold"],
                  a_option=settings["a_option"])
        fft = _c_argv(output=folder + "/fft_stage", method=3, start=0, count=count, **kw)
        if not dry_run:
            docker(container, "env", "LPSDCFN=" + root + "/lpsd.cfg", *fft)
        path = output_dir / (channel + ".native.ps")
        if not dry_run and path.exists():
            raise FileExistsError(path)
        previous, nrows = 0.0, 0
        handle = None if dry_run else path.open("x")
        try:
            for i in range(batches):
                a, b = i * count // batches, (i + 1) * count // batches
                raw = folder + f"/{channel}.{i}.raw7"
                command = _c_argv(output=raw, method=2, start=a, count=b-a, **kw)
                if dry_run:
                    continue
                docker(container, "env", "LPSDCFN=" + root + "/lpsd.cfg", *command)
                extractor = (
                    'import sys\n'
                    'with open(sys.argv[1]) as stream:\n'
                    ' for line in stream:\n'
                    '  if not line.strip() or line.lstrip().startswith("#"): continue\n'
                    '  columns=line.split()\n'
                    '  if len(columns)!=7: raise ValueError("Not native fDSNbri")\n'
                    '  print(columns[0],columns[2])\n')
                p = subprocess.Popen(["docker", "exec", container, "python3", "-c",
                                      extractor, raw], stdout=subprocess.PIPE, text=True)
                part = 0
                for line in p.stdout:
                    pair = line.split()
                    if len(pair) != 2:
                        raise ValueError("Bad native C PS transport")
                    import math
                    frequency, power = map(float, pair)
                    if not (math.isfinite(frequency) and math.isfinite(power)
                            and frequency > previous and power >= 0):
                        raise ValueError("Bad C PS frequency, power, or ordering")
                    previous = frequency
                    handle.write(line)
                    part += 1
                if p.wait() != 0 or part != b-a:
                    raise RuntimeError(f"Invalid native C batch: {channel}/{i}")
                nrows += part
            if not dry_run and nrows != count:
                raise RuntimeError("Incorrect total native C rows")
        finally:
            if handle is not None:
                handle.close()
        spectra[channel] = str(path)
    return spectra


def run(config_path, *, dry_run=False, resume=False):
    path = Path(config_path).expanduser().resolve(strict=True)
    cfg = json.loads(path.read_text())
    if cfg.get("scientific_scope") != "provisional_real_data_only":
        raise ValueError("This is not an authorized physical exclusion measurement")
    if cfg.get("signal_model") != "existing_vendor_single_bin":
        raise ValueError("No unapproved multibin halo substitution")
    physics_keys = {"reference_coupling_gev_inv", "quadrature_power_fraction",
                    "finesse_s", "tau_s", "tau1_s", "rho_dm_gev4", "unit_conversion"}
    if set(cfg.get("physics", {})) != physics_keys:
        raise ValueError("Explicit already-existing APE signal model parameters required")
    out = Path(cfg["output_dir"]).expanduser().resolve()
    if out.exists():
        if not resume or dry_run:
            raise FileExistsError("Existing output; use --resume for checked downstream restart")
        snapshot = out / "orchestration_inputs.json"
        if not snapshot.is_file():
            raise RuntimeError("No native-PS checkpoint; refusing unsafe C restart")
        saved = json.loads(snapshot.read_text())
        if saved["config_sha256"] != sha256(path):
            raise RuntimeError("Resume configuration hash mismatch")
        verified = choose_calibration(cfg["source"], dry_run=False)
        if saved["calibrated_hdf5_sha256"] != verified["hdf5_sha256_recorded"]:
            raise RuntimeError("Resume calibrated input hash mismatch")
        for channel, file in saved["downstream"]["native_ps"].items():
            if sha256(file) != saved["native_ps_sha256"][channel]:
                raise RuntimeError("Resume native C PS hash mismatch")
        downstream = out / "downstream"
        if (downstream / "summary.json").is_file():
            return {"status": "ALREADY_COMPLETED_REPORT_REUSED", "output": str(downstream)}
        from tools.ape_real_full_pipeline import run as downstream_run
        result = downstream_run(saved["downstream"], snapshot)
        return {"status": "PROVISIONAL_DOWNSTREAM_RESUMED", "downstream": result}
    calibrated = choose_calibration(cfg["source"], dry_run=dry_run)
    if dry_run and calibrated.get("status") == "calibration_planned":
        return {"status": "DRY_RUN_CALIBRATION_REQUIRED", "calibration": calibrated}
    native = cfg["native_ps"]
    if native["mode"] == "reuse":
        paths = reuse_native(native, calibrated)
    elif native["mode"] == "run_docker":
        if dry_run:
            run_native_c(native, calibrated, out / "native_c", dry_run=True)
            return {"status": "DRY_RUN_C_PLAN_OK", "output": str(out)}
        out.mkdir(parents=True, exist_ok=False)
        paths = run_native_c(native, calibrated, out / "native_c", dry_run=False)
    else:
        raise ValueError("native_ps.mode must be reuse or run_docker")
    if dry_run:
        return {"status": "DRY_RUN_INPUTS_VERIFIED", "native_ps": paths}
    downstream = dict(
        scientific_scope="real_data_provisional_diagnostic",
        input_unit="native_PS_rad2", peak_shape="vendor_single_bin",
        external_vendor_root=cfg["external_vendor_root"], native_ps=paths,
        output_dir=str(out / "downstream"), background=cfg["background"],
        physics=cfg["physics"], optimizer=cfg["optimizer"],
        local_sigma_threshold=cfg["local_sigma_threshold"],
        scan_stride=cfg.get("scan_stride", 1))
    if not out.exists():
        out.mkdir(parents=True, exist_ok=False)
    snapshot = out / "orchestration_inputs.json"
    with snapshot.open("x") as h:
        json.dump({"config_sha256": sha256(path),
                   "calibrated_hdf5_sha256": calibrated["hdf5_sha256_recorded"],
                   "physical_limits_validated": False, "downstream": downstream,
                   "native_ps_sha256": {ch: sha256(f) for ch, f in paths.items()}},
                  h, indent=2)
    from tools.ape_real_full_pipeline import run as downstream_run
    downstream_result = downstream_run(downstream, snapshot)
    return {"status": "PROVISIONAL_REAL_DATA_RUN_COMPLETE",
            "downstream": downstream_result,
            "physical_exclusion_validated": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    a = parser.parse_args()
    print(json.dumps(run(a.config, dry_run=a.dry_run, resume=a.resume), indent=2))


if __name__ == "__main__":
    main()
