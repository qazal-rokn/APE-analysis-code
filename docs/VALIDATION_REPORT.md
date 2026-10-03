# Controlled Synthetic Validation Report

Date: 2026-09-24. Scope: temporary synthetic HDF5 only; no APE cache samples.

- A float64 HDF5 round trip preserved sample order and values.
- Synthetic masks spanning two blocks selected the globally longest valid run and
  did not stitch across an invalid sample.
- The actual C executable reported Kaiser `SBIN=8.15` under the approved window
  settings; the wrapper uses it to derive the C feasibility lower bound.
- A deterministic 18 Hz tone was evaluated on a deliberately coarse four-node
  C validation grid. Its maximum PSD was in the returned C frequency region at
  14.286609468 Hz. This validates regional recovery and C output behavior, not
  a precise-frequency or production-grid sensitivity claim.
- Seeded white noise produced finite nonnegative PSD values and positive C
  average counts. This is a basic stochastic sanity check, not a full ensemble
  normalization study.
- A four-node one-job C result exactly matched two consecutive two-node C jobs;
  frequencies were ordered and unique.
- Parser checks accepted strict `fDNRb` rows and rejected malformed rows. ASD is
  defined only as `sqrt(PSD)`; no second normalization is applied.
- Separate `rho_x` tone and `rho_y` noise paths generated separate C outputs and
  separate PSD/ASD traces, with no X/Y combination.
- Temporary unsmoothed log-log PSD and ASD plots rendered successfully.

Remaining gates before real data: complete production-scale resource benchmark;
tone amplitude/normalization study on the approved dense grid; white-noise
ensemble confidence/tolerance study; persistent non-overwriting result/artifact
writer and provenance audit; full cache segment inventory; comparison review;
physics/statistics review; final audit; and separate explicit real-data approval.

Dense synthetic validation is run by `tools/run_dense_synthetic_benchmark.py`.
It records deterministic near-low/mid/high tone responses, four seeded noise
realizations, increasing small-grid C jobs, batch equivalence, configuration,
commands, hashes, logs, and non-overwriting run records. Noise acceptance is
provisional: finite nonnegative PSD and reproducible ensemble band summaries;
production tolerances require a larger approved ensemble.
