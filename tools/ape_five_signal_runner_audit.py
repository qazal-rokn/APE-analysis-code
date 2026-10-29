#!/usr/bin/env python3
"""Read-only Stage 0–7 status and plot-ready archive for the existing five-signal fixture.

No scientific computation, original C execution, vendor refitting or data modifications.
Only small synthetic outputs are copied to an ephemeral GitHub Actions artifact.
The original HDF5, FFT caches and archival stage files remain on the Mac.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil

HOME = Path.home()
ARCHIVE = HOME / "Desktop/APE_Research/APE_40MIN_FIVE_SIGNAL_SIMULATION"
STAGE3 = ARCHIVE / "03_Native_C_LPSD"
STAGE4 = ARCHIVE / "04_Background_Fitting"
FREQS = {3200: 31, 4000: 26, 4800: 23, 5900: 19, 6800: 18}
INPUT_SHA = "919b38ef7e14e78b82066ce87688f0315fc87840cec6622954e7c964a854a6d7"
BINARY_SHA = "cb0b1df10c9863145ba8b35df60e075eac3dc08965f3ad832857a7d4fec5eec9"
FIG_NAMES = (
    "01_ps_and_background_full",
    "02_ps_and_background_50mHz",
    "03_log_residual_vs_frequency",
    "04_residual_histogram_skewnormal",
)
EXPORT_SUFFIXES = {".csv", ".json", ".npz", ".png", ".pdf", ".py", ".md", ".log", ".cfg", ".lpsd"}
ALLOWED_NATIVE = {"native_1600bins.lpsd", "stage03_reference.native.lpsd"}
ALLOWED_NATIVE.update(f"{ch}_{f}.lpsd" for ch in ("rho_x", "rho_y") for f in FREQS)
MAX_ARTIFACT_FILE_BYTES = 32 * 1024 * 1024


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def require(cond, message):
    if not cond:
        raise RuntimeError(message)


def check_csv(path: Path, *, count: int, required: tuple[str, ...]) -> None:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        require(all(n in (reader.fieldnames or []) for n in required),
                f"Missing CSV columns: {path}")
        prev = -math.inf
        total = 0
        for row in reader:
            f_hz = float(row["frequency_hz"])
            require(math.isfinite(f_hz) and f_hz > prev,
                    f"Bad frequency order: {path}")
            prev = f_hz
            if "native_ps_rad2" in row:
                power = float(row["native_ps_rad2"])
                bkg = float(row["background_ps_rad2"])
                resid = float(row["log_ps_over_background"])
                require(all(map(math.isfinite, (power, bkg, resid)))
                        and power > 0 and bkg > 0
                        and math.isclose(resid, math.log(power / bkg),
                                         rel_tol=1e-8, abs_tol=1e-8),
                        f"Invalid saved fit values: {path}")
            total += 1
    require(total == count, f"Expected {count} rows in {path}, got {total}")


def export_selected(source: Path, target: Path, manifest: list) -> None:
    require(source.is_dir(), f"Missing export source {source}")
    for p in sorted(source.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(source)
        if p.suffix not in EXPORT_SUFFIXES:
            continue
        if p.suffix == ".lpsd" and p.name not in ALLOWED_NATIVE:
            continue
        if p.stat().st_size > MAX_ARTIFACT_FILE_BYTES:
            print(f"SKIP_LARGE_FILE_RETAINED_ON_MAC: {p} ({p.stat().st_size} bytes)")
            continue
        if p.name == "file_inventory.csv":
            # Preserve the original inventory as one of the independently saved reports.
            pass
        dst = target / rel
        require(not dst.exists(), f"Refusing to overwrite {dst}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, dst)
        require(sha(dst) == sha(p), f"Copy mismatch: {p}")
        manifest.append({
            "stage_relative_source": str(p.relative_to(ARCHIVE)),
            "artifact_relative_path": str(dst.relative_to(next(p for p in target.parents if p.name == "ape_five_signal_report"))),
            "sha256": sha(dst),
            "bytes": dst.stat().st_size,
        })


def check_stage3() -> dict:
    source = STAGE3 / "native_c"
    report = load_json(source / "stage03_sparse_validation.json")
    require(report.get("status") == "sparse_native_ps_validated", "Stage 3 status mismatch")
    require(report.get("input_sha256") == INPUT_SHA, "Stage 3 source hash mismatch")
    require(report.get("c_binary_sha256") == BINARY_SHA, "Stage 3 C hash mismatch")
    require(report.get("frequency_grid_J") == 847300, "Stage 3 grid mismatch")
    expected = {(c, f) for c in ("rho_x", "rho_y") for f in FREQS}
    actual = {(x["channel"], int(x["target_hz"])) for x in report["windows"]}
    require(len(report["windows"]) == 10 and actual == expected, "Stage 3 window set mismatch")
    for entry in report["windows"]:
        ch, freq = entry["channel"], int(entry["target_hz"])
        native = source / ch / f"{ch}_{freq}.lpsd"
        require(sha(native) == entry["sha256"], f"Stage 3 SHA mismatch: {native}")
        total, last, positive = 0, -math.inf, 0
        for line in native.read_bytes().splitlines():
            if not line.strip() or line.lstrip().startswith(b"#"):
                continue
            columns = line.split()
            require(len(columns) == 7, f"Stage 3 native columns: {native}")
            vals = [float(x) for x in columns]
            require(all(map(math.isfinite, vals)) and vals[0] > last
                    and vals[1] >= 0 and vals[2] > 0,
                    f"Stage 3 native PS invalid: {native}")
            last = vals[0]
            total += 1
            positive += vals[2] > 0
        require(total == FREQS[freq] and int(entry["rows"]) == total and positive == total,
                f"Stage 3 incomplete: {native}")
    plot_data = STAGE3 / "plot_data"
    for freq, count in FREQS.items():
        p = plot_data / f"XY_{freq}Hz.csv"
        check_csv(p, count=count, required=("frequency_hz", "PS_X_rad2",
                                          "PS_Y_rad2", "PS_X_plus_Y_rad2"))
        with p.open(newline="") as f:
            for row in csv.DictReader(f):
                x, y, total = (float(row[v]) for v in
                               ("PS_X_rad2", "PS_Y_rad2", "PS_X_plus_Y_rad2"))
                require(x > 0 and y > 0 and math.isclose(x + y, total, rel_tol=1e-9),
                        f"Stage 3 plotted X/Y sum mismatch {p}")
    require((STAGE3 / "plot_stage03.py").is_file(), "Stage 3 replot script missing")
    return {"status": "pass", "native_ps_windows": 10,
            "frequencies_hz": list(FREQS), "units": "rad^2",
            "note": "Native spectra computed, not five statistical detections"}


def check_fit(folder: Path, channel: str, freq: int) -> dict:
    report = load_json(folder / "vendor_fit_report.json")
    status = load_json(folder / "STATUS.json")
    require(report.get("channel") == channel and
            int(report.get("rest_frequency_hz", -1)) == freq and
            report.get("number_of_bins") == 1600 and
            report.get("vendor_quality_pass") is True, f"Fit report failed: {folder}")
    require(status.get("status") == "vendor_background_fit_pass",
            f"Fit status not pass: {folder}")
    chi = float(report["vendor_reduced_chi2"])
    require(math.isfinite(chi) and 0.5 <= chi <= 10.0, f"Invalid chi2 {folder}")
    require((folder / "vendor_fit_arrays.npz").is_file(), f"NPZ missing {folder}")
    require((folder / "plot_vendor_background.py").is_file(), f"Replot missing {folder}")
    check_csv(folder / "native_ps_background_residuals.csv", count=1600,
              required=("frequency_hz", "native_ps_rad2",
                        "background_ps_rad2", "log_ps_over_background"))
    for name in FIG_NAMES:
        for ext in (".png", ".pdf"):
            require((folder / "figures" / (name + ext)).is_file(),
                    f"Stage 4 figure missing: {folder}/{name}{ext}")
    return {"channel": channel, "frequency_hz": freq, "chi2": chi,
            "fit_directory": str(folder.relative_to(ARCHIVE)),
            "status": "pass"}


def check_stage4() -> tuple[dict, list[Path]]:
    pilot = STAGE4 / "pilot_3200Hz_X_4f653750c3" / "vendor_fit_3ea3188d76"
    batch = STAGE4 / "remaining_nine_1ccb351c7e"
    summary = load_json(batch / "stage04_batch_summary.json")
    require(summary.get("status") == "ten_vendor_fits_pass"
            and summary.get("all_vendor_quality_pass") is True
            and summary.get("new_fits_completed") == 9
            and summary.get("total_fits_including_pilot") == 10,
            "Stage 4 ten-fit batch summary is not validated")
    require(summary.get("input_sha256") == INPUT_SHA
            and summary.get("executable_sha256") == BINARY_SHA,
            "Stage 4 source hash differs")
    require(load_json(batch / "STAGE_STATUS.json").get("status") == "pass",
            "Stage 4 batch STATUS not pass")
    results = [check_fit(pilot, "rho_x", 3200)]
    for ch in ("rho_x", "rho_y"):
        for f in FREQS:
            if (ch, f) == ("rho_x", 3200):
                continue
            results.append(check_fit(batch / f"{'X' if ch == 'rho_x' else 'Y'}_{f}Hz", ch, f))
    require({(r["channel"], r["frequency_hz"]) for r in results}
            == {(ch, f) for ch in ("rho_x", "rho_y") for f in FREQS},
            "Missing Stage 4 channel/frequency")
    return {"status": "pass", "vendor_background_fits": 10,
            "fits": results, "note": "Background fits only; likelihood pending"}, [pilot, batch]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--destination", type=Path, required=True)
    opts = ap.parse_args()
    dest = opts.destination.expanduser().resolve()
    require(not dest.exists(), f"Refusing overwrite of {dest}")
    require(ARCHIVE.is_dir(), f"Missing Mac archive: {ARCHIVE}")
    stage3 = check_stage3()
    stage4, fit_sources = check_stage4()
    dest.mkdir(parents=True, exist_ok=False)
    manifest = []
    for name in ("00_Simulation_Input", "01_Calibration", "02_Cleaning"):
        src = ARCHIVE / name
        if src.is_dir():
            export_selected(src, dest / name, manifest)
    export_selected(STAGE3 / "native_c", dest / "03_Native_C_LPSD" / "native_c", manifest)
    export_selected(STAGE3 / "plot_data", dest / "03_Native_C_LPSD" / "plot_data", manifest)
    script = STAGE3 / "plot_stage03.py"
    target = dest / "03_Native_C_LPSD" / "plot_stage03.py"
    shutil.copyfile(script, target)
    manifest.append({"stage_relative_source": str(script.relative_to(ARCHIVE)),
                     "artifact_relative_path": str(target.relative_to(dest)),
                     "sha256": sha(target), "bytes": target.stat().st_size})
    export_selected(STAGE4 / "pilot_3200Hz_X_4f653750c3", dest / "04_Background_Fitting" / "pilot_3200_X", manifest)
    export_selected(fit_sources[1], dest / "04_Background_Fitting" / "remaining_nine", manifest)
    report = {
        "experiment": "synthetic_40min_five_signal_20261029",
        "full_pipeline_complete": False,
        "source_mac_archive": str(ARCHIVE),
        "source_hdf5_sha256": INPUT_SHA,
        "c_binary_sha256": BINARY_SHA,
        "stage_status": {
            "00_Simulation_Input": {"status": "fixture_exists_not_rerun"},
            "01_Calibration": {"status": "original_APE_execution_unverified"},
            "02_Cleaning": {"status": "original_APE_execution_unverified"},
            "03_Native_C_LPSD": stage3,
            "04_Background_Fitting": stage4,
            "05_Likelihood": {"status": "not_executed", "blocker": "validated signal power and one-bin versus multi-bin scientific contract"},
            "06_Candidates": {"status": "not_executed", "blocker": "requires approved local likelihood"},
            "07_Upper_Limits": {"status": "not_executed", "blocker": "physical normalization, coverage and boundary convention unverified"},
        },
        "interpretation": "Synthetic software diagnostics only; no APE discovery or observed physical exclusion",
    }
    (dest / "simulation_stage_status.json").write_text(json.dumps(report, indent=2) + "\n")
    (dest / "file_inventory.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (dest / "file_inventory.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=("stage_relative_source", "artifact_relative_path", "sha256", "bytes"))
        writer.writeheader()
        writer.writerows(manifest)
    lines = ["# APE 40-minute five-injection synthetic fixture: archived stage report", "",
             "## Stage status", ""]
    for stage, value in report["stage_status"].items():
        lines.append(f"- **{stage}:** {value['status']}")
    lines += ["", "Stages 3–4 are checked against saved numerical artifacts and fits.",
              "Stages 1–2 are not yet independently confirmed as passing the original APE adapters.",
              "Stages 5–7 are not executed; global significance and physical limits are not validated.",
              "The Mac holds the original HDF5 and FFT caches. GitHub artifact contains plot-ready synthetic data.",
              "Plotting scripts are saved beside the scientific CSV and NPZ products.", ""]
    (dest / "README.md").write_text("\n".join(lines))
    print(json.dumps({"status": "APE_FIVE_SIGNAL_READONLY_STAGE_AUDIT_PASS",
                      "destination": str(dest), "exported_files": len(manifest),
                      "stage03": "pass", "stage04": "pass",
                      "stage05_through_07": "not_executed"}, indent=2))


if __name__ == "__main__":
    main()
