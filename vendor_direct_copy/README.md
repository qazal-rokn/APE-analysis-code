# APE vendor-direct-copy v1 (isolated third implementation)

This branch starts from `main`, NOT from the earlier experimental adapters.

## Immutable upstream boundary
The existing APE ingestion, optical calibration and original C-LPSD execution
remain unchanged. Consume only their native X and Y C-LPSD outputs.

## Vendor source policy
Copy the complete original `scalar_dark_matter_lpsd/lvk/` Python source
tree and `LPSDIO.py` into `vendor_original/` without editing, retaining
original attribution/license. Record SHA-256 hashes. No AST extraction, no
reimplemented candidate or upper-limit statistics, and no monkeypatches.

## Authorized APE-specific modifications
1. Add a separate post-C-LPSD entrypoint that sums X/Y **linear power**
   frequency-bin by frequency-bin, with exact native-grid checks.
2. Replace the LIGO instrument response with APE optical response, only at
   the response boundary.
3. Exclude H1/L1-specific dual-interferometer consistency/veto and plots.
4. Retain all remaining original vendor background, likelihood, candidate,
   upper-limit algorithms, parameters and sequencing.

## Status
BRANCH CREATED; ORIGINAL SOURCE COPY AND INTEGRATION NOT YET COMPLETED.
Do not run this branch for scientific limits or treat it as validated.
