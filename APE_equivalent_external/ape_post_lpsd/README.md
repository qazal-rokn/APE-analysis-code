# Step 1 after C-LPSD: source-preserving APE spectrum handoff

This folder is an **APE-only adapter** added to the external-equivalent working
copy. Neither the original external C-LPSD nor any external Python statistical
source file is edited.

## What is preserved

* **Two independent spectra:** `rho_x` and `rho_y` are not summed.
* **Original native C frequencies:** no regenerated logarithmic grid, interpolation
  or rebinning. X and Y must have numerically identical frequency arrays.
* **Original one-sided C PSD:** no new calibration, normalization or averaging.
  The original C `PSD`, `N_avg`, and `N_kernel` are saved as float64 datasets.
* **Unchanged external reader:** the vendor `LPSDOutput` reads datasets
  `frequency` and `logPSD`. The latter is exactly `np.log(PSD)`, consistent
  with its own text-input transformation.
* **Provenance:** input SHA-256 values, sample rate, channel labels, and source
  segment identifiers are recorded in the HDF5 roots and a JSON report.
* **No overwrites:** both output channels are staged before final filenames.

## Input contract

Each source is a **combined four-column C method-2 spectrum**, where rows have
`frequency_Hz  PSD_rad2_per_Hz  N_avg  N_kernel`. This is the format
produced by the existing APE Constant-Q combiner. Zero or nonfinite PSD cannot
be passed through the unchanged vendor log-PSD background model and is rejected,
not silently repaired. The frequency grid must be strictly increasing and not
exceed Nyquist.

Do not feed a six-column raw C batch file directly into this adapter.

## Execute on an approved local spectrum

From the repository's root:

```bash
python3 APE_equivalent_external/ape_post_lpsd/prepare_vendor_lpsd.py \
  --x-psd /path/to/combined_rho_x.lpsd \
  --y-psd /path/to/combined_rho_y.lpsd \
  --output-dir /path/to/new/ape_vendor_handoff \
  --fs 20000 \
  --t0 ACTUAL_SEGMENT_START_INTEGER \
  --t1 ACTUAL_SEGMENT_STOP_INTEGER \
  --expected-bins 5035956
```

Replace the paths and `t0`/`t1` with verified metadata; do not invent
timestamps. The adapter does **not** run or modify the original C estimator.
It writes files `APE_<t0>_<t1>_X.h5` and
`APE_<t0>_<t1>_Y.h5` for use by the vendor IO and background fitter.

For the original unchanged background fitting program, pass **one resulting
HDF5 file at a time** to `--data-path` and send JSON fit output to a
separate path. Do not pass the handoff directory indiscriminately because it
also contains the provenance JSON.

## Validate

From `APE_equivalent_external/` in an environment with NumPy, h5py, SciPy,
pandas and tqdm:

```bash
python3 -m unittest discover -s ape_post_lpsd/tests -v
```

The tests cover synthetic X/Y preservation, metadata, malformed inputs,
frequency mismatch, length mismatch, overwrite refusal, and compatibility with
the original external `LPSDOutput` and `LPSDDataGroup` interfaces.

**Scope:** verified synthetic handoff and original-reader compatibility only.
No real 5,035,956-bin APE handoff or background fit was executed for this
change. The original likelihood still expects LIGO H1/L1 response metadata;
that is a **later step**, not changed here.


## Optional APE path using the *original* external converter

`prepare_vendor_lpsd_original.py` is a separate experiment that routes the
APE four-column spectra through the **unchanged**
`scalar_dark_matter_lpsd.LPSDIO.LPSDOutput` and its original
`save_to_hdf5()` method. The adapter supplies the vendor text header and
tab delimiter expected by that reader; it does not change spectral mathematics.
The original reader reconstructs the logarithmic grid with `np.logspace`.

Unlike `prepare_vendor_lpsd.py`, this path deliberately uses the vendor's
reconstructed frequencies. It produces a JSON comparison of those frequencies
against the original C output, with maximum absolute/relative differences.
**Evaluate this comparison on the full real APE spectra before replacing the
current production handoff.** The test with a deliberately non-logarithmic
three/five-bin grid illustrates possible differences, not a demonstrated
difference in the production grid.

```bash
python3 APE_equivalent_external/ape_post_lpsd/prepare_vendor_lpsd_original.py \
  --x-psd /path/to/combined_rho_x.lpsd \
  --y-psd /path/to/combined_rho_y.lpsd \
  --output-dir /path/to/new/original_vendor_handoff \
  --fs 20000 --t0 VERIFIED_T0 --t1 VERIFIED_T1 \
  --expected-bins 5035956
```

The original converter's HDF5 files contain only `frequency` and `logPSD`,
with original vendor metadata. Ancillary `N_avg` and `N_kernel` remain
in source spectra; they are not copied into these HDF5 files.
The adapter omits only optional `None` metadata attributes unsupported by
HDF5, preserving the original frequency and PSD calculations.
Keep X and Y separate. No vendor code was edited and the existing handoff
was not removed. Synthetic tests: `tests/test_prepare_vendor_lpsd_original.py`.
