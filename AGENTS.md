# APE C-LPSD Comparison Project

## Core rule

The researcher is the final decision-maker.

For any scientific or substantial code change:

1. Inspect first.
2. Explain what currently happens.
3. Propose the change.
4. Report which files would be affected.
5. Wait for explicit researcher approval.
6. Implement only the approved change.
7. Validate it.
8. Report exactly what changed.

## Safety

- Never modify `/Users/user/Desktop/Code/Projects/APE-data-analysis` or `/Users/user/APE-main-code-inbox`.
- Never modify the C-LPSD source or executable unless explicitly instructed.
- Never modify, overwrite, or delete raw experimental data or calibrated cache files.
- Never silently change calibration constants, equations, statistical methods, units, conventions, channel definitions, or cleaning rules.
- Never silently smooth, filter, interpolate, resample, fill, reject, join, or remove data.
- Never install packages without approval.
- Never execute hardware-control code.
- Never perform unrelated refactoring.

## Scientific workflow

Scientific correctness is more important than code elegance. Report improvements
before implementing them, using:

```text
PROBLEM:
PROPOSED CHANGE:
REASON:
FILES AFFECTED:
```

Then wait for approval.

## Reporting

After implementation, report:

```text
REQUESTED:
CHANGED:
FILES MODIFIED:
VALIDATION:
NOT CHANGED:
WARNINGS / UNCERTAINTIES:
```

## Mandatory implementation rules

### A. Confirm before coding

No scientific analysis code may be written or modified until the researcher
explicitly approves an Analysis Contract based on
`docs/ANALYSIS_CONTRACT_TEMPLATE.md`.

### B. No method substitution

Implement the approved method exactly. The compiled C `lpsd-exec` must genuinely
be invoked; a Python LPSD, FFT, Welch estimate, sinusoidal fit, copied spectrum,
or relabeled result is not a substitute.

### C. No helpful extras

Do not add filters, smoothing, cuts, normalization, interpolation, resampling,
gap filling, refactoring, extra estimators, or interpretation without approval.

### D. Small modules

Each module must have defined inputs, outputs, units, scientific purpose, and a
validation method.

### E. Test before real data

New modules must pass controlled synthetic validation. Applying the estimator to
real APE data requires a later explicit researcher approval.

### F. Centralized settings

Scientific constants and configurable parameters must come from the approved
configuration system. Record every effective C command-line option.

### G. Scientific equations

Authoritative equations belong in `src/physics/`. Do not silently change them.

### H. Documentation synchronization

Keep method documentation and `docs/EQUATION_REFERENCE.tex` synchronized with
the approved implementation and document assumptions and limitations.

### I. Audit

Audit USER REQUEST -> APPROVED CONTRACT -> CODE -> VALIDATION -> DOCUMENTATION.
Report every mismatch.

### J. Scope and immutability

Implement only the requested scope. Treat all reference projects, notebooks, C
sources/binaries, raw data, cache files, and prior APE results as read-only.

