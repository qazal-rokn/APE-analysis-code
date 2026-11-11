# APE external-equivalent working copy (baseline)

This branch starts from the verified, unmodified external program:
- Upstream: https://github.com/AlexandreGoettel/Scalar-Dark-Matter-LPSD
- Pinned upstream commit: `f9e78cbe1031406e217c0bfb63a2f6f0043a0a5a`
- Exact upstream Git tree SHA: `237b3ba91539ca0ac188cc78c2906fb035b6c301`
- Unmodified reference: `vendor_direct_copy/original_program/`
- Independently editable working copy: `APE_equivalent_external/`

**Baseline status:** As created, the working copy's Git tree is exactly the
upstream tree, including all 68 tracked files, bytes and executable modes.
**No APE-related implementation change has been made yet.**

## Change boundary

1. Keep the immutable `vendor_direct_copy/original_program/` untouched.
2. Only edit `APE_equivalent_external/` for changes demonstrably required by
   APE detector response, detector channels and input/output adaptation.
3. Preserve the original background-model, skew-normal residual distribution,
   likelihood, signed parameter optimization, Fisher uncertainty, candidate
   threshold/trial correction, and upper-limit statistics exactly wherever
   physically applicable. Do not silently repair original bugs.
4. Do not reinterpret LIGO's H1/L1 consistency tests as equivalent to APE
   quadratures. Any necessary exception must be explicitly justified.
5. Before scientific use, document each difference against the pinned vendor
   source and run quantitative parity, injection and end-to-end tests.
6. Keep all original attribution, LICENSE and other license notices.

## Status

This is a **source-identical starting baseline**, not a certified APE-equivalent
pipeline and not a demonstrated end-to-end execution.
