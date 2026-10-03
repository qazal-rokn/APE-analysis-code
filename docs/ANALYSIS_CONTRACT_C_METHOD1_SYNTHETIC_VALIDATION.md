# Analysis Contract: Synthetic C Method-1 Validation

## Approved scope

Compare the unchanged `lpsd-exec` direct LPSD path (`-h 0`) with the
program's documented FFT-approximation path (`-h 1`) using generated data
only.  This contract authorizes neither real-data use of `-h 1` nor a change
to the frozen direct-production profile.

## Matched calculation

Each method receives the same generated float64 HDF5 input, sample rate,
time interval, frequency limits, `Jdes`, output format, Kaiser window
(`-w -2`, `-p 200`), linear-regression option (`-r 1`), overlap (`-l 50`),
averaging options (`-a 100`, `-m 1`), scale (`-x 1`), and memory option
(`-M 20`).  Method 1 additionally receives an explicit `-E 1` percent
approximation tolerance.  The latter has no direct-method counterpart.

Three frequency regions are tested independently: low, middle, and high.
In each region the generated inputs are: a tone on an actual C-generated
node, seeded white noise, and their sum.  The study records node identity,
PSD and ASD values, recovered tone peak, noise summaries, wall time, child
maximum RSS, stdout/stderr, commands, hashes, and parsed outputs.

## Acceptance is not pre-approved

No numerical tolerance is pre-approved for scientific use.  The resulting
report will propose evidence-based tolerances separately for: exact generated
node identity, tone-peak location, relative PSD difference, and white-noise
band statistics.  A researcher decision remains required before `-h 1` can
be selected for any real data.

## Exclusions

No real APE signal/cache/raw data is read.  No existing HDF5 export, grid
manifest, C source/executable, calibration, cleaning definition, or frozen
real-data result is changed.
