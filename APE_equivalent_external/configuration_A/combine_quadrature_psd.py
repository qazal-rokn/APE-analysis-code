"""Configuration A: sum calibrated X/Y *linear PSDs* before vendor fitting.

Input: two accepted APE handoff HDF5 files, exact same native C-LPSD grid.
Output: one APE_<t0>_<t1>_Q.h5 read directly by original LPSDOutput and
fit_background.py. Does not alter vendor source or refit anything by itself.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np


class CombinedPSDInputError(ValueError):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def combine_psd_hdf5(x_path, y_path, output_dir, *, chunk_size=16384):
    x_path, y_path, output_dir = (Path(p).expanduser().resolve()
                                  for p in (x_path, y_path, output_dir))
    if x_path == y_path or not x_path.is_file() or not y_path.is_file():
        raise CombinedPSDInputError("Require two distinct existing X/Y handoff files")
    if not isinstance(chunk_size, int) or chunk_size < 1:
        raise CombinedPSDInputError("chunk_size must be positive")
    with h5py.File(x_path, "r") as x, h5py.File(y_path, "r") as y:
        for data, channel in ((x, "X"), (y, "Y")):
            if data.attrs.get("channel") != channel:
                raise CombinedPSDInputError("Channel labels do not match X/Y files")
            if data.attrs.get("spectrum_unit") != "rad^2/Hz":
                raise CombinedPSDInputError("Both channels must be calibrated to rad^2/Hz")
            if not all(k in data for k in ("frequency", "PSD", "logPSD", "N_avg", "N_kernel")):
                raise CombinedPSDInputError("Missing original native-grid handoff datasets")
        common = ("segment_t0", "segment_t1", "sample_rate_hz")
        if any(x.attrs.get(a) != y.attrs.get(a) for a in common):
            raise CombinedPSDInputError("X/Y segment or sampling metadata mismatch")
        t0, t1 = int(x.attrs["segment_t0"]), int(x.attrs["segment_t1"])
        n = len(x["frequency"])
        if n < 2 or any(len(data[k]) != n for data in (x,y)
                         for k in ("frequency", "PSD", "logPSD", "N_avg", "N_kernel")):
            raise CombinedPSDInputError("Incompatible X/Y native frequency sizes")
        name = f"APE_{t0}_{t1}_Q.h5"
        dest = output_dir / name
        report = output_dir / f"APE_{t0}_{t1}_Q_sum_report.json"
        stage = output_dir / f"{name}.partial"
        if any(p.exists() for p in (dest, report, stage)):
            raise CombinedPSDInputError("Combined output already exists; refusing overwrite")
        output_dir.mkdir(parents=True, exist_ok=True)
        previous = None
        try:
            with h5py.File(stage, "x") as out:
                out.attrs["channel"] = "Q"
                out.attrs["spectrum_unit"] = "rad^2/Hz"
                out.attrs["construction"] = "PSD_Q = PSD_X + PSD_Y, BEFORE log and vendor background fit"
                for a in common:
                    out.attrs[a] = x.attrs[a]
                for name_ds in ("frequency", "PSD", "logPSD", "N_avg", "N_kernel"):
                    out.create_dataset(name_ds, shape=(n,), dtype="f8",
                                       chunks=(min(chunk_size,n),))
                for start in range(0, n, chunk_size):
                    stop = min(start+chunk_size, n)
                    fx, fy = np.asarray(x["frequency"][start:stop]), np.asarray(y["frequency"][start:stop])
                    if not np.array_equal(fx, fy) or not np.isfinite(fx).all() or np.any(fx <= 0) or np.any(np.diff(fx) <= 0):
                        raise CombinedPSDInputError("X/Y native grids differ, are not finite, or are not increasing")
                    if previous is not None and fx[0] <= previous:
                        raise CombinedPSDInputError("Frequency ordering across chunks failed")
                    previous = fx[-1]
                    ax, ay = np.asarray(x["N_avg"][start:stop]), np.asarray(y["N_avg"][start:stop])
                    kx, ky = np.asarray(x["N_kernel"][start:stop]), np.asarray(y["N_kernel"][start:stop])
                    if (not np.array_equal(ax, ay) or not np.array_equal(kx, ky)
                            or not np.isfinite(ax).all() or np.any(ax <= 0)
                            or not np.isfinite(kx).all() or np.any(kx < 0)):
                        raise CombinedPSDInputError("Different estimator averaging/window metadata: cannot sum blindly")
                    px, py = np.asarray(x["PSD"][start:stop]), np.asarray(y["PSD"][start:stop])
                    if (not np.isfinite(px).all() or not np.isfinite(py).all()
                            or np.any(px <= 0) or np.any(py <= 0)):
                        raise CombinedPSDInputError("Invalid linear PSD inputs")
                    combined = px + py
                    if not np.isfinite(combined).all() or np.any(combined <= 0):
                        raise CombinedPSDInputError("Invalid combined PSD")
                    if (not np.allclose(np.asarray(x["logPSD"][start:stop]), np.log(px), atol=1e-12, rtol=1e-12)
                            or not np.allclose(np.asarray(y["logPSD"][start:stop]), np.log(py), atol=1e-12, rtol=1e-12)):
                        raise CombinedPSDInputError("Input logPSD and PSD inconsistent")
                    for label, values in (("frequency",fx), ("PSD",combined),
                                          ("logPSD",np.log(combined)),
                                          ("N_avg",ax), ("N_kernel",kx)):
                        out[label][start:stop] = values
                for key,val in x["logPSD"].attrs.items():
                    if key in y["logPSD"].attrs and y["logPSD"].attrs[key] != val:
                        raise CombinedPSDInputError("X/Y vendor LPSD metadata differ")
                    out["logPSD"].attrs[key] = val
            stage.rename(dest)
            details = {
                "status": "SUMMED_PSD_REQUIRES_NEW_VENDOR_BACKGROUND_FIT",
                "new_channel": "Q",
                "construction": "PSD_Q(f)=PSD_X(f)+PSD_Y(f); logPSD_Q=log(PSD_Q)",
                "n_bins": n,
                "native_grid_preserved": True,
                "inputs": {"X":str(x_path),"Y":str(y_path)},
                "input_sha256": {"X":sha256(x_path),"Y":sha256(y_path)},
                "output_file": str(dest),
                "output_sha256": sha256(dest),
                "fitting_requirement": "New original-vendor background fit on Q. Do not reuse background_X or background_Y",
                "statistical_model": "Original background fitting and skew-normal likelihood require fit diagnostics and injection validation",
            }
            with report.open("x", encoding="utf-8") as stream:
                json.dump(details, stream, indent=2)
                stream.write("\n")
            return details
        except BaseException:
            stage.unlink(missing_ok=True)
            dest.unlink(missing_ok=True)
            raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--x-h5", required=True)
    p.add_argument("--y-h5", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--chunk-size", type=int, default=16384)
    a = p.parse_args()
    print(json.dumps(combine_psd_hdf5(a.x_h5,a.y_h5,a.output_dir,
                                      chunk_size=a.chunk_size),indent=2))


if __name__ == "__main__":
    main()
