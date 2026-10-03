# Scientific Changelog

## 2026-09-24 — Project scaffold and proposed contract

- Created a separate project scaffold and governance roles.
- Recorded reference-project, cache-format, and C-source inspection findings.
- Added a proposed Analysis Contract and validation plan.
- Added no scientific implementation.
- Ran no synthetic or real-data spectrum.

## 2026-09-24 — Approved synthetic-only C-LPSD bridge

- Added configuration-driven cache inspection, global valid-run selection,
  float64 HDF5 transport, C command construction, output parser, and CLI guard.
- Approved C profile: Kaiser -2/200 dB, linear regression, 50% overlap,
  desired averages 100, minimum averages 1, scale 1, `fDNRb`, and Q=1e6.
- Retained `fmin=fs/N`; changed minimum averages from the earlier proposed 10 to
  approved 1 because the C program rejects the former combination.
- Added synthetic HDF5/order, no-gap-stitching, longest-across-block, and grid tests.
- No real cache spectrum was run.

## 2026-09-24 — Approved C-window feasibility bound

- Replaced the infeasible theoretical C lower limit with an automatic effective
  C-window bound. `1/T` remains a diagnostic frequency.
- The wrapper runs a temporary synthetic C probe, parses the actual configured
  `SBIN`, and uses the exact C lower-bound equation. No value such as 8.15 is
  embedded as a production constant.
- Synthetic C probe confirmed `SBIN=8.15` for Kaiser `-w -2`, PSLL 200 dB,
  50% overlap, and minimum averages 1.

## 2026-09-24 — Controlled compiled-C synthetic checks

- Added temporary tone/noise C executions, strict `fDNRb` parsing checks,
  single-batch versus split-batch equivalence, separate X/Y artifacts, and
  unsmoothed PSD/ASD rendering.

## 2026-09-25 — approved C-native 65–9000 Hz grid planner

- Added a no-data planner that mirrors unchanged `lpsd.c` method-0 node and DFT-length equations, selects the densest feasible `Jdes` for the approved validated segment, and writes immutable manifests.
- No signal data, HDF5 transport, C spectrum, calibration, cleaning, or estimator settings changed.

## 2026-09-25 — Future full-record FFT extension scaffold

- Added a separate, unexecuted full-record one-sided voltage FFT extension.
- It is explicitly distinct from C-LPSD; no C-window equivalence is claimed.
- Added configuration, conservative memory preflight design, future artifact
  design, and an unexecuted validation plan. No experimental signal was read,
  no HDF5 was written, and no spectrum was calculated.

## 2026-09-25 — Ingestion/cleaning-review/segment-selection scaffold

- Added manual and chronology-gated automatic ingestion architecture.
- Automatic mode calls the configured reference Modules 2A, 2B, and 3 rather
  than substituting cleaning/calibration equations.
- Added explicit reference-default YAML settings, source-manifest requirement,
  cross-file continuity rule, frozen segment-manifest gate, and future-only
  diagnostic/test scaffolds. No raw data were opened or analyzed.

## 2026-09-25 — Approved calibrated quadrature review

- Added a read-only review calculation that uses the frozen reference scalar
  denominator separately for `pem1_x` and `pem1_Y` with the exact 1000:1
  raw/reduced mapping.
- Added a synthetic, no-Dewesoft reference-equivalence test: the quadrature
  norm agrees with reference `rho` within `rtol=1e-14`, `atol=1e-15`.
- Added non-overwriting X/Y ellipticity and denominator review plots plus
  calibration-profile provenance. No phase/gain/leakage correction, X/Y
  combination, HDF5, FFT, C-LPSD, candidates, or limits was added.
