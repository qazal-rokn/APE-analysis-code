# Minimal-difference audit: APE-to-vendor post-LPSD input (Step 1)

## Scope and pinned sources

- External source: `AlexandreGoettel/Scalar-Dark-Matter-LPSD` at commit `f9e78cbe1031406e217c0bfb63a2f6f0043a0a5a`.
- APE source: `qazal-rokn/APE-analysis-code` main at commit `03722edb975e91a3653d3da673e5a521bb3a2854`.
- Adapted copy: `APE_equivalent_external/ape_post_lpsd/prepare_vendor_lpsd.py`.
- Compared interfaces: `scalar_dark_matter_lpsd/LPSDIO.py`, especially `LPSDOutput.get_lpsd_kwargs`, `read`, `read_hdf5`, and `LPSDDataGroup`; `scalar_dark_matter_lpsd/lvk/fit_background.py`.
- Original vendor statistical code and C estimator were **not modified**.

## Direct evidence

1. APE's verified method-2 combiner writes **exactly four columns**:
   frequency, one-sided PSD, N_avg, N_kernel. Its output contains no vendor
   command-line header. The original `LPSDOutput` text reader rejects a
   synthetic representative file with `OSError: Invalid output file.`.
2. Adding a compatible command-line header is **not enough to guarantee native
   frequency identity**: the original text reader parses PSD from column 2,
   then derives a fresh `np.logspace(fmin,fmax,Jdes)` grid. A controlled
   three-bin spectrum with native frequencies [1,2,3] was read back as
   [1,sqrt(3),3]. This is a deliberately non-logarithmic counterexample,
   demonstrating loss of native-grid fidelity. Whether a specific APE
   production grid happens to agree with the reconstruction still needs a
   full-grid comparison.
3. The original HDF5 reader accepts `frequency` and `logPSD` datasets,
   reading the frequency vector **directly** rather than regenerating it.
   HDF5 metadata attributes `fmin`, `fmax`, `fs`, `Jdes` are required
   for constructing `LPSDVars`. The unchanged `LPSDDataGroup` recognizes
   filenames ending in `_<t0>_<t1>_<ifo>.h5`.
4. The APE adapter writes this native interface without rebinning, rescaling,
   applying a detector response, summing X/Y, or editing the vendor source.
   The input PSD remains available in an extra HDF5 dataset, with
   `logPSD=np.log(PSD)` for the unchanged vendor reader.

## Necessity classification

| Feature | Scientific/runtime necessity | Decision |
|---|---|---|
| Convert four-column APE spectrum into a vendor-readable format | **Necessary** for current APE files without changing vendor IO | Keep |
| Preserve the actual native C frequency vector | **Necessary** to prevent silent frequency changes | Keep HDF5 native-frequency interface |
| Store `frequency`, `logPSD`, and required HDF5 attributes | **Necessary** for original vendor reader | Keep |
| Keep X and Y distinct at the handoff | **Necessary at present**; any combination requires separate statistical justification | Keep |
| Preserve `PSD`, `N_avg`, `N_kernel` as auxiliary HDF5 datasets | Not required by vendor `LPSDOutput`; useful scientific metadata and reproducibility | Keep as non-invasive safeguards |
| SHA-256 provenance, strict input checks, non-overwrite, temp files | Not required by likelihood mathematics; engineering safety/audit safeguards | Keep outside vendor code |
| Filename convention `APE_<t0>_<t1>_X/Y.h5` | **Necessary** for existing vendor `sep_label` and background JSON labels | Keep |
| Refit background, alter likelihood, C-LPSD, or rescale PSD | **Not needed** for this interface step | Prohibited here |

## Actual tests executed in local environment

The imported Step-1 adapter was checked byte-for-byte against its GitHub source
(blob `5d3913278e00b3d6e8571154ba1af5b0973b21c9`).
The unchanged original LPSDIO implementation was also checked against the
upstream blob (SHA `3f0138e7f93e4a126c30598fadbe5e9c38df83a4`).
Nine synthetic tests passed:
- Independent native-grid/PSD/metadata round-trip
- Unchanged vendor `LPSDOutput` / `LPSDDataGroup` accept the adapter output
- Original vendor text-reader rejection of headerless APE data
- Counterexample proving that header-only text changes the native grid
- Reject mismatched X/Y frequencies
- Reject mismatched row count
- Reject zero/nonfinite PSD
- Reject malformed column count
- Refuse overwriting existing results

No full APE production-spectrum conversion or numerical background fitting was
run in this audit. An existing miniaturized helper package was used for the
same import signatures in the temporary local test environment, with the real,
SHA-matching unchanged `LPSDIO.py`. An end-to-end checkout test on the full
external package remains to be performed in a configured environment.

## Minimality conclusion

**Required new scientific algorithms in Step 1: zero.**
The only genuinely required integration is a format bridge that exposes the
unchanged C-LPSD data in the original vendor HDF5 interface, preserving native
frequency coordinates. The current implementation additionally includes
non-numerical guardrails and provenance, not new scientific algorithms. These
extra protections are not *mathematically necessary* and should never be
presented as such. The code has not been demonstrated to have the minimum
possible line count or the minimum storage footprint.

Step 2 (background fitting) is **not automatically certified** by this result.
It requires its own controlled invocation and tests using original vendor
fitting functions, without changing their mathematics.
