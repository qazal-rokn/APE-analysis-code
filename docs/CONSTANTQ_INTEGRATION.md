# Constant-Q controller V1 (synthetic-only integration)

Source contract: docs/ANALYSIS_CONTRACT_APE_CONSTANTQ_INTEGRATION.md.

The original compiled external C estimator remains unchanged. The reference
scientific profile is frozen at 20kHz, 65–10000 Hz, 5,035,956 bins and its
original 16 uneven batches. One method-3 FFT command per quadrature feeds
16 method-2 batch commands. X and Y are never time-domain combined.

Modules:
- spectral/constantq.py: immutable profile, exact index plan and C argv.
- spectral/combine.py: six-column method-2 input validation, four-column export.
- backends/local.py: guarded synthetic execution with controlled LPSDCFN/logs.
- backends/slurm.py: pure preview generation with FFT dependency requirement.
- tests/test_constantq_integration.py: small fixtures, command/format checks.

**Not authorized**: real-data processing, cluster submission, modifying the
external executable, or overwriting archived spectra. No science results are
claimed by these source files. Actual C execution and historical command-byte
equivalence remain to be tested in an approved environment.

The compatibility profile intentionally preserves the historical 15435.0s
C endpoint rather than silently using the full 15435.85s calibrated input.
No method-0 flags/outputs are accepted as modern method-2 products.

Install/use: place the repository `src` directory on PYTHONPATH and import
`ape_c_lpsd.spectral.constantq`. Further CLI integration remains pending
until an exact synthetic execution contract is checked against the vendor binary.
